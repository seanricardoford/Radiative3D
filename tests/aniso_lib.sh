# tests/aniso_lib.sh
#
# Shared helpers for the anisotropic-scale-length regression tests. This
# file is meant to be sourced (". tests/aniso_lib.sh"), not executed
# directly -- it only defines functions, and relies on the caller having
# already done "set -eu".
#
# Functions:
#   normalize_stdout FILE
#   normalize_stdout_noargs FILE
#   compare_runs DIR_A DIR_B MODE      (MODE is "strict" or "noargs")
#   trace_has_signal FILE
#   require_sct_bands SCENARIO REPORT BAND...

##
## normalize_stdout FILE
##
##   Prints FILE's content with heap addresses masked (0x... -> 0xADDR) and
##   the BUILD STATS line removed. The address mask covers both the
##   scatterer dump's "address" column (printed on every run, with or
##   without --dump-grid) and any other 0x... token. The BUILD STATS line
##   embeds the git revision string, which differs between the baseline and
##   current builds, so it must not be part of a byte-for-byte comparison.
##
normalize_stdout() {
  sed -E -e 's/0x[0-9a-f]+/0xADDR/g' -e '/BUILD STATS/d' "$1"
}

##
## normalize_stdout_noargs FILE
##
##   As normalize_stdout(), but additionally masks out:
##     - the line right after 'Args provided for "compiled" model are:'
##       (the echoed model-args vector itself), and
##     - the "with N args" count printed while constructing the grid.
##   This lets a comparison tolerate a different number/value of model args
##   between two runs (e.g. once a later task adds anisotropic scale-length
##   args), while this task's own isotropic regression test still uses the
##   stricter normalize_stdout().
##
normalize_stdout_noargs() {
  normalize_stdout "$1" | sed -E \
    -e '/Args provided for "compiled" model are:/{n;d;}' \
    -e 's/with [0-9]+ args/with N args/'
}

##
## compare_runs DIR_A DIR_B MODE
##
##   Byte-for-byte compares the outputs of two run_aniso_case.sh runs found
##   in directories DIR_A and DIR_B.
##
##   Always compared:
##     - seis_traces_asc.dat, with "cell: 0x...." addresses masked the same
##       way test_lopnor_thread_equivalence.sh does.
##     - every seis_*.octv file, compared verbatim.
##
##   MODE "strict" additionally compares out_mparams.octv verbatim and
##   compares normalize_stdout() output of stdout.txt.
##
##   MODE "noargs" skips out_mparams.octv (which echoes the model-args
##   vector, so it cannot match once that vector changes) and instead
##   compares normalize_stdout_noargs() output of stdout.txt.
##
compare_runs() {
  _car_dir_a=$1
  _car_dir_b=$2
  _car_mode=$3

  sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$_car_dir_a/seis_traces_asc.dat" \
    > "$_car_dir_a/.seis_traces_asc.normalized"
  sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$_car_dir_b/seis_traces_asc.dat" \
    > "$_car_dir_b/.seis_traces_asc.normalized"
  # Every cmp below is checked explicitly and returns on the first mismatch,
  # rather than relying on the caller's "set -e" to stop us: a caller that
  # invokes compare_runs inside "if"/"&&"/"||" (a natural way to use it)
  # would otherwise suspend errexit for this whole function, silently
  # masking all but the last comparison.
  cmp "$_car_dir_a/.seis_traces_asc.normalized" "$_car_dir_b/.seis_traces_asc.normalized" || return 1

  for _car_f in "$_car_dir_a"/seis_*.octv; do
    _car_base=$(basename "$_car_f")
    cmp "$_car_f" "$_car_dir_b/$_car_base" || return 1
  done

  case "$_car_mode" in
    strict)
      cmp "$_car_dir_a/out_mparams.octv" "$_car_dir_b/out_mparams.octv" || return 1
      normalize_stdout "$_car_dir_a/stdout.txt" > "$_car_dir_a/.stdout.normalized"
      normalize_stdout "$_car_dir_b/stdout.txt" > "$_car_dir_b/.stdout.normalized"
      cmp "$_car_dir_a/.stdout.normalized" "$_car_dir_b/.stdout.normalized" || return 1
      ;;
    noargs)
      normalize_stdout_noargs "$_car_dir_a/stdout.txt" > "$_car_dir_a/.stdout.normalized"
      normalize_stdout_noargs "$_car_dir_b/stdout.txt" > "$_car_dir_b/.stdout.normalized"
      cmp "$_car_dir_a/.stdout.normalized" "$_car_dir_b/.stdout.normalized" || return 1
      ;;
    *)
      echo "compare_runs: unknown mode '$_car_mode' (want strict|noargs)" >&2
      return 1
      ;;
  esac
}

##
## trace_has_signal FILE
##
##   True (exit 0) if FILE contains at least one nonzero value inside a
##   trace block. Identical all-zero traces would make compare_runs pass
##   vacuously, so callers should require this on at least one side of
##   every comparison. Copied verbatim from
##   tests/test_lopnor_thread_equivalence.sh (that file is owned by another
##   worker and is not touched here).
##
trace_has_signal() {
  awk '
    /#### BEGIN TRACE ####/ { in_trace = 1; next }
    /#### END TRACE ####/ { in_trace = 0; next }
    in_trace { for (i = 1; i <= NF; i++) if ($i + 0 != 0) found = 1 }
    END { exit !found }
  ' "$1"
}

##
## require_sct_bands SCENARIO REPORT BAND...
##
##   Counts "SCT:" scatter-report lines in REPORT whose location falls in
##   each named BAND, and fails -- naming the offending band -- if any
##   count is zero. REPORT must come from a run with REPORTS=SCT (so that
##   reports.dat contains only SCT: lines).
##
##   SCENARIO selects the coordinate interpretation of each SCT line's
##   "xyz:( x y z )" triple:
##
##     lopnor:      the band test applies directly to internal z (the 3rd
##                   number). The Lop Nor grid uses a flat (ENU) mapping
##                   with Earth-flattening, and the simulator's default
##                   output transform un-flattens z but leaves x,y alone,
##                   so z here is already the depth-equivalent coordinate.
##                   The "sediments" band additionally caps range =
##                   sqrt(x^2+y^2) at 700: each layer boundary is a plane
##                   fit through 3 widely-spaced control points (Lop Nor,
##                   MAK, WUS -- see user_LopNorCyl_inc.cpp / RCUCylinder),
##                   and far from those points (beyond ~800, within the
##                   model's 1200 range cylinder) that plane's extrapolation
##                   drifts enough to put z > 1.0 below the true local
##                   sediment/crust interface; the "crust" and "mantle"
##                   bands were checked the same way and keep a comfortable
##                   (>8 km) margin at every range, so they need no cap.
##
##     crustpinch:  the band test applies to depth = -z and
##                   range = sqrt(x^2+y^2). The crust-pinch grid uses a
##                   curved (RAE_CURVED) mapping, and the simulator's
##                   default output transform (ecs.cpp's OutConvert(),
##                   case RAE_CURVED/OUT_ENU_ORTHO) already converts the
##                   true internal xyz into a curvature-corrected
##                   range/azimuth projection for x,y and true elevation
##                   for z -- confirmed by reproducing Convert()+
##                   OutConvert() by hand and matching --dump-grid's
##                   printed node coordinates to 5 decimal places. So depth
##                   and range are already un-distorted here and need no
##                   further spherical correction.
##                     Applying a *second* spherical correction on top of
##                   this -- i.e. treating the reported x,y,z as if they
##                   were still raw/uncorrected internal Cartesian
##                   coordinates, via depth = 6371 - |xyz - (0,0,-6371)| --
##                   double-counts the curvature and produces wrong,
##                   sometimes negative, depths away from range zero. E.g.
##                   a real scatter event at true depth 5.1 km, range 504 km
##                   (reported xyz = 425.8, 269.4, -5.12) evaluates to a
##                   sane depth=5.1 under the simple formula used here, but
##                   to a nonsensical depth=-14.8 under that double
##                   correction.
##
require_sct_bands() {
  _rsb_scenario=$1
  _rsb_report=$2
  shift 2

  for _rsb_band in "$@"; do
    _rsb_count=$(awk -v scenario="$_rsb_scenario" -v band="$_rsb_band" '
      /^SCT:/ {
        found_xyz = 0
        for (i = 1; i <= NF; i++) {
          if ($i == "xyz:(") { x = $(i+1) + 0; y = $(i+2) + 0; z = $(i+3) + 0; found_xyz = 1; break }
        }
        if (!found_xyz) next

        hit = 0
        if (scenario == "lopnor") {
          range = sqrt(x * x + y * y)
          if (band == "sediments")    { if (z > 1.0 && range < 700) hit = 1 }
          else if (band == "crust")   { if (z > -30 && z < -15)     hit = 1 }
          else if (band == "mantle")  { if (z < -60)                hit = 1 }
        } else if (scenario == "crustpinch") {
          depth = -z
          range = sqrt(x * x + y * y)
          if (band == "sediments")    { if (depth < 1.5) hit = 1 }
          else if (band == "crust")   { if (depth > 12 && depth < 20 && (range < 300 || range > 550)) hit = 1 }
          else if (band == "pinched") { if (depth > 12 && depth < 20 && range > 380 && range < 460) hit = 1 }
          else if (band == "moho")    { if (depth > 33 && depth < 36 && (range < 300 || range > 550)) hit = 1 }
          else if (band == "mantle")  { if (depth > 50) hit = 1 }
        }
        if (hit) count++
      }
      END { print count + 0 }
    ' "$_rsb_report")

    if [ "$_rsb_count" -eq 0 ]; then
      echo "require_sct_bands: no SCT hits in band '$_rsb_band' for scenario '$_rsb_scenario'" >&2
      return 1
    fi
  done
}
