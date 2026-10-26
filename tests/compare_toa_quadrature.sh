#!/bin/sh
set -eu
#
# tests/compare_toa_quadrature.sh
#
# Manual scenario-level validation of --toa-quadrature=area against the
# default --toa-quadrature=equal: how much do whole-simulation outputs
# (station energies) move, how much do per-scatterer mean free paths
# shift, and what does area mode cost to build at a high take-off-angle
# degree. Not named test_*, so the full test-suite check skips it; run it
# directly. See docs/superpowers/plans/2026-10-25-toa-area-quadrature.md
# (Task 5) and docs/superpowers/specs/2026-10-25-toa-area-
# quadrature-design.md.
#
# Environment, with defaults:
#   TQ_DEGREE               --toa-degree for the energy/MFP-table runs.
#                                                             Default: 7
#   TQ_REPLICATES            K, replicates per mode per config. Must be
#                             >= 2 (need a sample variance per group).
#                                                             Default: 12
#   TQ_THREADS               --threads for every run, or "auto" for the
#                             logical CPU count.              Default: auto
#   TQ_SEED_BASE              Base seed added to (i*2K + k)*2^32 for each
#                             configuration i / replicate k (see "Seeds"
#                             below).                         Default: 596
#   TQ_LOPNOR_PHONONS         --num-phonons for lopnor runs.   Default: 200000
#   TQ_CRUSTPINCH_PHONONS     --num-phonons for crustpinch runs. Default: 50000
#   TQ_CONFIGS                Which (scenario, variant) configs to run for
#                             the energy and MFP tables; see syntax below.
#            Default: "lopnor:iso lopnor:axi4 crustpinch:iso crustpinch:axi4"
#   TQ_BUILD_CHECK            1 to also run the degree-9 build check below,
#                             0 to skip it.                    Default: 1
#
# TQ_CONFIGS syntax: a space-separated list of "SCENARIO:VARIANT" tokens,
# e.g.
#   TQ_CONFIGS="lopnor:iso lopnor:axi10"
# SCENARIO is "lopnor" or "crustpinch" (tests/run_aniso_case.sh's
# scenarios). VARIANT is "iso" (scenario_model SCENARIO 0.01 iso) or
# "axiR" for an integer R (scenario_model SCENARIO 0.01 axi R) -- e.g.
# "axi4" is tests/aniso_lib.sh's mode "axi" with R=4. A single config,
# e.g. TQ_CONFIGS="lopnor:iso", is how the smoke test below exercises the
# script cheaply.
#
# Seeds: for the configuration at 0-based position i in TQ_CONFIGS,
# replicate k = 1..K of that config runs in equal mode with seed
# TQ_SEED_BASE + (i*2K + k)*2^32; replicate k = K+1..2K runs in area mode
# (same formula). The i*2K term gives every configuration its own block
# of seeds, disjoint from every other configuration's, so no two
# configurations -- and, within one configuration, no two replicates or
# modes -- ever draw from the same seed; this matters because the energy
# table below treats different configurations as statistically
# independent samples. 2^32 is far larger than any phonon count used
# here, so MixRandomSeed(seed ^ phonon_id) (model.cpp) never lets two
# replicates' per-phonon streams collide (doc/MANUAL.md's --seed entry).
# Computed with awk ('%.0f' on an IEEE double), exact as long as the
# result stays under 2^53 -- enforced below by rejecting any
# TQ_SEED_BASE/TQ_REPLICATES/config-count combination whose largest seed
# would reach or exceed it, so (unlike tests/compare_axi_general.sh's
# seed_for(), which must also tolerate a TQ_SEED_BASE near 2^64) no
# arbitrary-precision arithmetic is needed here. The degree-9 build check
# (below) is not one of TQ_CONFIGS's own configurations; it reuses this
# same seed_for() with i=0, which collapses the formula back to plain
# TQ_SEED_BASE + k*2^32, since its three replicates per cell never need
# to be independent of the energy/MFP configurations above.
#
# Energy table, per config: K equal-mode and K area-mode runs are each
# reduced, via extract_total_energy() below, to one "station energy"
# number per station per run -- the sum of that station's X+Y+Z
# axis-energy trace columns over every trace bin of the whole run (the
# same per-bin quantity tests/compare_axi_general.sh's extract_energy()
# sums, but totaled over the whole run instead of split into 10-second
# windows: a whole-run total is never legitimately zero the way an
# individual window can be, so no positivity filter is needed). For each
# station, let meanE/meanA be the mean over the K equal/area replicates
# and seE/seA their standard errors (sample standard deviation / sqrt(K)).
# r_s = meanA/meanE, and SE(r_s) is obtained by first-order error
# propagation treating the two means as independent (they come from
# disjoint seed blocks):
#   SE(r_s) = r_s * sqrt((seA/meanA)^2 + (seE/meanE)^2)
# The printed row covers every station of the config's scenario: station
# count, min/median/max r_s across stations, and how many stations have
# |r_s - 1| > 3*SE(r_s).
#
# Energy table, per-station detail (spec decision 11): the summary row
# alone can't be interpreted -- a median r_s a long way from 1 could be
# one wildly-off station or a genuine broad shift, and the two look the
# same in a single summary number. So, per config, an appendix table
# gives one row per station: station id, mean(E, equal), mean(E, area),
# r_s, SE(r_s), and z = (r_s - 1)/SE(r_s) (how many of its own standard
# errors a station's ratio sits from 1; the summary row's "flagged"
# count is just the stations with |z| > 3 here). A final "all stations"
# row gives the scenario-wide ratio R: for each replicate, sum that
# run's energy over every station first, then treat those K per-replicate
# whole-scenario totals exactly like a single station's K values --
# R = mean(sum_s E_s, area) / mean(sum_s E_s, equal), with SE(R) by the
# same error-propagation formula as SE(r_s) above. R answers "does area
# mode move the scenario's total energy", which is a different (and
# usually much less noisy) question than "does it move any one station's
# share of it" -- the per-station rows answer that one.
#
# MFP table, per config: the first equal-mode run's and first area-mode
# run's stdout.txt each contain a scatterer dump (the lines between
# "BEGIN SCATTERER DUMP" and "END SCATTERERS" that don't start with "#");
# both dumps list the same scatterers in the same order (the dump order
# depends only on construction order, not on --toa-quadrature), so row i
# of one matches row i of the other. Each printed row holds the matched
# pair's a (really ax; column 3), ay, az and mode (columns 12-14; "-" for
# an isotropic row, which has none of those columns) and, for both MFP(P) and
# MFP(S) (columns 8 and 9), the equal value, the area value, and the
# percent change (area/equal - 1)*100.
#
# Build check (TQ_BUILD_CHECK=1): at take-off-angle degree 9 and
# --num-phonons=0 (so only model/grid/TOA-array construction runs, not
# phonon casting), Lop Nor "iso" and "axi10" are each run 3 times in
# equal mode and 3 times in area mode, single data point per run: the
# wall-clock seconds and maximum resident set size (RSS) of the "main"
# process itself. R3D_EXE is set to a generated wrapper script
# ($TIME_WRAPPER, below) that execs
#   /usr/bin/time -l "$ROOT/main" "$@"            (Darwin; this machine)
#   /usr/bin/time -v "$ROOT/main" "$@"            (elsewhere)
# with time's own stderr (which also carries anything main itself writes
# to stderr, since the child inherits the wrapper's redirected fd 2)
# redirected to time_output.txt in the run directory -- tests/
# run_aniso_case.sh always cd's into OUTDIR before running R3D_EXE, so a
# relative filename lands in the right place without the wrapper needing
# to know OUTDIR's path. parse_time_output() (below) extracts the wall
# seconds ("real", macOS; "Elapsed (wall clock) time", Linux, converted
# from [h:]mm:ss) and the RSS in bytes ("maximum resident set size",
# macOS, already bytes; "Maximum resident set size", Linux, kbytes*1024).
# The printed table gives, per variant, the median (of 3) wall seconds
# and median RSS for each mode and the area-minus-equal delta of each.
#
# Streams: progress goes to stderr; the two (or three, with the build
# check) Markdown tables go to stdout. The script exits nonzero if any
# run fails (tests/run_aniso_case.sh itself reports a nonzero exit) or if
# a run's stdout.txt lacks __SIMULATION_COMPLETE__ -- including a
# --num-phonons=0 build-check run, which this task confirmed still
# prints that marker (Model::RunSimulation() prints it unconditionally
# once the phonon loop -- zero iterations or not -- and all worker
# threads have joined).

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
. "$ROOT/tests/aniso_lib.sh"

TQ_DEGREE=${TQ_DEGREE:-7}
TQ_REPLICATES=${TQ_REPLICATES:-12}
TQ_THREADS=${TQ_THREADS:-auto}
TQ_SEED_BASE=${TQ_SEED_BASE:-596}
TQ_LOPNOR_PHONONS=${TQ_LOPNOR_PHONONS:-200000}
TQ_CRUSTPINCH_PHONONS=${TQ_CRUSTPINCH_PHONONS:-50000}
TQ_CONFIGS_DEFAULT="lopnor:iso lopnor:axi4 crustpinch:iso crustpinch:axi4"
TQ_CONFIGS=${TQ_CONFIGS:-$TQ_CONFIGS_DEFAULT}
TQ_BUILD_CHECK=${TQ_BUILD_CHECK:-1}

TQ_EPS=0.01   # Fixed eps for every config and for the build check, per
              # the brief's "axiR means scenario_model SCENARIO 0.01 axi
              # R" / "iso means scenario_model SCENARIO 0.01 iso".

if [ "$TQ_REPLICATES" -lt 2 ]; then
  echo "ERROR: TQ_REPLICATES must be >= 2 (need a sample variance per group)" >&2
  exit 1
fi

case "$TQ_SEED_BASE" in
  ''|*[!0-9]*)
    echo "ERROR: TQ_SEED_BASE must be a non-negative integer (got \"$TQ_SEED_BASE\")" >&2
    exit 1
    ;;
esac

# Reject a TQ_SEED_BASE/TQ_REPLICATES/TQ_CONFIGS combination whose largest
# seed (configuration index TQ_NUM_CONFIGS-1, k=2K; see the header
# comment's "Seeds" section) would reach 2^53, where seed_for()'s
# IEEE-double arithmetic (below) stops being exact.
TQ_NUM_CONFIGS=$(set -- $TQ_CONFIGS; echo $#)
awk -v base="$TQ_SEED_BASE" -v n="$TQ_NUM_CONFIGS" -v k="$TQ_REPLICATES" '
  BEGIN {
    max_seed = base + n*2*k*4294967296
    if (max_seed >= 9007199254740992) {
      printf "ERROR: TQ_SEED_BASE=%s, TQ_REPLICATES=%s and %d TQ_CONFIGS " \
             "entries would use seeds up to %.0f, at or above 2^53 " \
             "(9007199254740992), where seed arithmetic in this script " \
             "stops being exact; reduce TQ_SEED_BASE, TQ_REPLICATES, or " \
             "the number of TQ_CONFIGS entries\n", \
             base, k, n, max_seed > "/dev/stderr"
      exit 1
    }
  }
' || exit 1

if [ "$TQ_THREADS" = "auto" ]; then
  TQ_THREADS=$(auto_threads)
fi

echo "Building current executable..." >&2
make >/dev/null

R3D_EXE_CURRENT="$ROOT/main"
if [ ! -x "$R3D_EXE_CURRENT" ]; then
  echo "ERROR: build did not produce an executable at $R3D_EXE_CURRENT" >&2
  exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

K=$TQ_REPLICATES

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

# seed_for CFG_INDEX K_OFFSET -> prints TQ_SEED_BASE +
# (CFG_INDEX*2K + K_OFFSET)*2^32 as an exact unsigned integer (see the
# header comment's "Seeds" section). K here is $K (TQ_REPLICATES, set
# above). The degree-9 build check (below) is not one of TQ_CONFIGS's own
# configurations; it calls this with CFG_INDEX=0, collapsing the formula
# back to plain TQ_SEED_BASE + K_OFFSET*2^32.
seed_for() {
  awk -v base="$TQ_SEED_BASE" -v cfgidx="$1" -v koff="$2" -v repk="$K" \
      'BEGIN{printf "%.0f", base + (cfgidx*2*repk + koff)*4294967296}'
}

# phonons_for_scenario SCENARIO -> prints the configured phonon count.
phonons_for_scenario() {
  case "$1" in
    lopnor) echo "$TQ_LOPNOR_PHONONS" ;;
    crustpinch) echo "$TQ_CRUSTPINCH_PHONONS" ;;
    *) echo "phonons_for_scenario: unknown scenario '$1'" >&2; exit 1 ;;
  esac
}

# tq_config_fields TOKEN -> prints "SCENARIO VARIANT" (space-separated)
# on stdout. Exits the whole script nonzero, with an error to stderr, if
# TOKEN is not exactly two non-empty, colon-separated fields.
tq_config_fields() {
  awk -v tok="$1" '
    BEGIN {
      n = split(tok, f, ":")
      if (n != 2 || f[1] == "" || f[2] == "") {
        printf "ERROR: TQ_CONFIGS token \"%s\" is not SCENARIO:VARIANT (two non-empty colon-separated fields)\n", tok > "/dev/stderr"
        exit 1
      }
      print f[1], f[2]
    }
  '
}

# model_for_variant SCENARIO VARIANT -> prints the --model-args value.
# VARIANT is "iso" or "axiR" (R a non-negative integer); see the header
# comment's "TQ_CONFIGS syntax".
model_for_variant() {
  _mfv_scenario=$1; _mfv_variant=$2
  case "$_mfv_variant" in
    iso)
      scenario_model "$_mfv_scenario" "$TQ_EPS" iso
      ;;
    axi*)
      _mfv_r=${_mfv_variant#axi}
      case "$_mfv_r" in
        ''|*[!0-9]*)
          echo "ERROR: TQ_CONFIGS variant \"$_mfv_variant\" is not 'iso' or 'axi' followed by an integer" >&2
          exit 1
          ;;
      esac
      scenario_model "$_mfv_scenario" "$TQ_EPS" axi "$_mfv_r"
      ;;
    *)
      echo "ERROR: TQ_CONFIGS variant \"$_mfv_variant\" is not 'iso' or 'axiR'" >&2
      exit 1
      ;;
  esac
}

# run_case SCENARIO OUTDIR PHONONS MODEL_ARGS SEED QUADMODE DEGREE [R3D_EXE]
#
#   Runs one tests/run_aniso_case.sh case (QUADMODE "equal" or "area"),
#   using R3D_EXE if given, else "$ROOT/main". Exits the whole script
#   nonzero -- after dumping whatever logs exist -- if the run itself
#   fails, or if it succeeds but stdout.txt lacks __SIMULATION_COMPLETE__.
run_case() {
  _rc_scenario=$1; _rc_outdir=$2; _rc_phonons=$3; _rc_model=$4
  _rc_seed=$5; _rc_quadmode=$6; _rc_degree=$7; _rc_exe=${8:-$R3D_EXE_CURRENT}

  _rc_extra=""
  if [ "$_rc_quadmode" = "area" ]; then
    _rc_extra="--toa-quadrature=area"
  fi

  if ! R3D_EXE="$_rc_exe" TOA_DEGREE="$_rc_degree" SEED="$_rc_seed" EXTRA_ARGS="$_rc_extra" \
        "$ROOT/tests/run_aniso_case.sh" "$_rc_scenario" "$TQ_THREADS" "$_rc_outdir" \
        "$_rc_phonons" "$_rc_model" \
        >"$_rc_outdir/harness-stdout.txt" 2>"$_rc_outdir/harness-stderr.txt"; then
    echo "ERROR: run failed (scenario=$_rc_scenario quadrature=$_rc_quadmode" \
         "seed=$_rc_seed degree=$_rc_degree outdir=$_rc_outdir)" >&2
    dump_failure_logs "$_rc_outdir"
    if [ -s "$_rc_outdir/time_output.txt" ]; then
      echo "----- last 20 lines of time_output.txt -----" >&2
      tail -n 20 "$_rc_outdir/time_output.txt" >&2
    fi
    exit 1
  fi

  if ! grep -q '__SIMULATION_COMPLETE__' "$_rc_outdir/stdout.txt"; then
    echo "ERROR: run completed but stdout.txt lacks __SIMULATION_COMPLETE__" \
         "(scenario=$_rc_scenario quadrature=$_rc_quadmode seed=$_rc_seed" \
         "degree=$_rc_degree outdir=$_rc_outdir)" >&2
    dump_failure_logs "$_rc_outdir"
    exit 1
  fi
}

# extract_total_energy SEIS_TRACES_FILE OUT_TSV
#
#   Writes "station<TAB>energy" lines to OUT_TSV: one per station, energy
#   being the sum of that station's X+Y+Z axis-energy trace columns over
#   every trace bin of the whole run (see the header comment's "Energy
#   table" section). Errors out (and so, under this script's "set -e",
#   aborts the whole run) if any station's total comes out non-positive,
#   since that would signal a station that never saw a phonon -- expected
#   behavior for an individual 10-second window (as in
#   tests/compare_axi_general.sh) but not for a whole 600-second run with
#   tens of thousands of phonons.
extract_total_energy() {
  awk '
    /Seismometer Number:/ {
      for (i=1;i<=NF;i++) if ($i=="Number:") { station=$(i+1)+0; break }
      next
    }
    /#### BEGIN TRACE ####/ { intrace=1; next }
    /#### END TRACE ####/   { intrace=0; next }
    intrace {
      sum[station] += ($1+$2+$3)
      next
    }
    END {
      for (s in sum) {
        if (sum[s] <= 0) {
          printf "extract_total_energy: non-positive total energy %.6e for station %s in %s\n", \
                 sum[s], s, FILENAME > "/dev/stderr"
          exit 1
        }
        printf "%s\t%.10e\n", s, sum[s]
      }
    }
  ' "$1" > "$2"
}

# energy_table_row EQLIST ARLIST K CFG DETAILOUT -> prints "stations
# min_r median_r max_r flagged" on stdout, for the existing summary row.
# EQLIST/ARLIST are files with K energy.tsv paths each, one per line (as
# written by extract_total_energy, one row per station per run). See the
# header comment's "Energy table" section for the summary statistic.
#
# As a side effect, writes CFG's per-station detail appendix rows (see
# the header comment's "Energy table, per-station detail" section) to
# DETAILOUT, as ready-to-print Markdown table rows: one per station
# (sorted by station number) plus a final "all stations" row.
energy_table_row() {
  _etr_eq=$1; _etr_ar=$2; _etr_k=$3; _etr_cfg=$4; _etr_detailout=$5
  awk -v K="$_etr_k" -v EQLIST="$_etr_eq" -v ARLIST="$_etr_ar" \
      -v CFG="$_etr_cfg" -v DETAILOUT="$_etr_detailout" '
    function stderr_of(n, sumv, sumsq,    m, var) {
      m = sumv / n
      var = (sumsq - n*m*m) / (n-1)
      if (var < 0) var = 0   # guards tiny negative from fp roundoff
      return sqrt(var / n)
    }
    # ratio_and_se(meanA, meanE, seA, seE) -> prints "r se_r" via global
    # out_r/out_se_r (awk has no multi-return; see call sites). r =
    # meanA/meanE; SE(r) by the same first-order error-propagation
    # formula the header comment gives for both r_s and R.
    function ratio_and_se(meanA, meanE, seA, seE) {
      out_r = meanA / meanE
      out_se_r = out_r * sqrt((seA/meanA)*(seA/meanA) + (seE/meanE)*(seE/meanE))
    }
    BEGIN {
      neq = 0
      while ((getline line < EQLIST) > 0) { neq++; eqfiles[neq] = line }
      close(EQLIST)
      nar = 0
      while ((getline line < ARLIST) > 0) { nar++; arfiles[nar] = line }
      close(ARLIST)
      if (neq != K || nar != K) {
        printf "energy_table_row: expected %d files per group, got %d and %d\n", \
               K, neq, nar > "/dev/stderr"
        exit 1
      }

      for (k=1; k<=K; k++) {
        eqtotal[k] = 0
        while ((getline line < eqfiles[k]) > 0) {
          split(line, f, "\t")
          s = f[1]; v = f[2] + 0
          eqsum[s] += v; eqsumsq[s] += v*v; eqn[s]++
          seen[s] = 1
          eqtotal[k] += v
        }
        close(eqfiles[k])
      }
      for (k=1; k<=K; k++) {
        artotal[k] = 0
        while ((getline line < arfiles[k]) > 0) {
          split(line, f, "\t")
          s = f[1]; v = f[2] + 0
          arsum[s] += v; arsumsq[s] += v*v; arn[s]++
          seen[s] = 1
          artotal[k] += v
        }
        close(arfiles[k])
      }

      nstations = 0
      flagged = 0
      for (s in seen) {
        if (eqn[s] != K || arn[s] != K) {
          printf "energy_table_row: station %s present in %d equal and %d area replicates (want %d each)\n", \
                 s, eqn[s]+0, arn[s]+0, K > "/dev/stderr"
          exit 1
        }
        meq = eqsum[s] / K
        mar = arsum[s] / K
        seeq = stderr_of(K, eqsum[s], eqsumsq[s])
        separ = stderr_of(K, arsum[s], arsumsq[s])

        ratio_and_se(mar, meq, separ, seeq)
        r = out_r; se_r = out_se_r

        nstations++
        stations[nstations] = s + 0
        ratio[nstations] = r
        diff = (r >= 1) ? r - 1 : 1 - r
        if (diff > 3*se_r) flagged++

        z = (se_r > 0) ? (r - 1)/se_r : 0   # se_r==0 needs identical
                                             # floating-point sums across
                                             # K independent-seed runs --
                                             # not observed in practice
        detail[s] = sprintf("| %s | %s | %.6e | %.6e | %.6f | %.6f | %+.3f |", \
                             CFG, s, meq, mar, r, se_r, z)
      }

      # "All stations" row: treat each replicate'"'"'s whole-scenario total
      # (summed over every station of that one run) as that replicate'"'"'s
      # single data point, K of them per mode, same formula as a station.
      eqt_sum = 0; eqt_sumsq = 0
      art_sum = 0; art_sumsq = 0
      for (k=1; k<=K; k++) { eqt_sum += eqtotal[k]; eqt_sumsq += eqtotal[k]*eqtotal[k] }
      for (k=1; k<=K; k++) { art_sum += artotal[k]; art_sumsq += artotal[k]*artotal[k] }
      meq_all = eqt_sum / K
      mar_all = art_sum / K
      seeq_all = stderr_of(K, eqt_sum, eqt_sumsq)
      separ_all = stderr_of(K, art_sum, art_sumsq)
      ratio_and_se(mar_all, meq_all, separ_all, seeq_all)
      R = out_r; se_R = out_se_r
      z_all = (se_R > 0) ? (R - 1)/se_R : 0

      # Sort stations ascending by number for a stable, readable,
      # reproducible appendix (awk'"'"'s for-in order is unspecified).
      for (i=2; i<=nstations; i++) {    # insertion sort; nstations is small
        keys = stations[i]; j = i-1
        while (j>=1 && stations[j] > keys) { stations[j+1]=stations[j]; j-- }
        stations[j+1] = keys
      }
      for (i=1; i<=nstations; i++) print detail[stations[i]] > DETAILOUT
      printf "| %s | all | %.6e | %.6e | %.6f | %.6f | %+.3f |\n", \
             CFG, meq_all, mar_all, R, se_R, z_all > DETAILOUT
      close(DETAILOUT)

      for (i=2; i<=nstations; i++) {    # separate sort of ratio[], for
        key = ratio[i]; j = i-1         # min/median/max below
        while (j>=1 && ratio[j] > key) { ratio[j+1]=ratio[j]; j-- }
        ratio[j+1] = key
      }
      if (nstations % 2 == 1) median = ratio[(nstations+1)/2]
      else                    median = (ratio[nstations/2] + ratio[nstations/2+1]) / 2

      printf "%d %.6f %.6f %.6f %d\n", nstations, ratio[1], median, ratio[nstations], flagged
    }
  '
}

# scatterer_dump_rows FILE -> prints FILE's scatterer-dump data rows
# (between "BEGIN SCATTERER DUMP" and "END SCATTERERS", excluding lines
# starting with "#"), one per line, in order. Pattern copied from
# tests/test_aniso_inputs.sh (that file is owned by another worker and is
# not touched here).
scatterer_dump_rows() {
  awk '
    /^#  BEGIN SCATTERER DUMP/ { f=1; next }
    /^#  END SCATTERERS/       { f=0; next }
    f && $0 !~ /^#/ { print }
  ' "$1"
}

# mfp_table_rows EQ_STDOUT AR_STDOUT CONFIG -> prints one Markdown row
# per matched scatterer-dump data-row pair (see the header comment's
# "MFP table" section) to stdout. Exits the whole script nonzero if the
# two dumps don't have the same number of data rows.
mfp_table_rows() {
  _mtr_eq=$1; _mtr_ar=$2; _mtr_cfg=$3
  _mtr_eqrows="$work/.mfp_eq_rows"
  _mtr_arrows="$work/.mfp_ar_rows"
  scatterer_dump_rows "$_mtr_eq" > "$_mtr_eqrows"
  scatterer_dump_rows "$_mtr_ar" > "$_mtr_arrows"

  _mtr_neq=$(wc -l < "$_mtr_eqrows" | tr -d ' ')
  _mtr_nar=$(wc -l < "$_mtr_arrows" | tr -d ' ')
  if [ "$_mtr_neq" -ne "$_mtr_nar" ]; then
    echo "ERROR: scatterer dump row count differs between equal ($_mtr_neq) and" \
         "area ($_mtr_nar) runs for config $_mtr_cfg" >&2
    exit 1
  fi

  paste "$_mtr_eqrows" "$_mtr_arrows" | awk -v cfg="$_mtr_cfg" -F'\t' '
    {
      neq = split($1, e, " ")
      split($2, a, " ")
      if (neq >= 14) { ay = e[12]; az = e[13]; mode = e[14] }
      else           { ay = "-";   az = "-";   mode = "-" }

      mfpP_eq = e[8] + 0; mfpP_ar = a[8] + 0
      mfpS_eq = e[9] + 0; mfpS_ar = a[9] + 0
      pchgP = (mfpP_eq != 0) ? (mfpP_ar/mfpP_eq - 1)*100 : 0
      pchgS = (mfpS_eq != 0) ? (mfpS_ar/mfpS_eq - 1)*100 : 0

      printf "| %s | %s | %s | %s | %s | %s | %s | %+.3f%% | %s | %s | %+.3f%% |\n", \
             cfg, e[3], ay, az, mode, e[8], a[8], pchgP, e[9], a[9], pchgS
    }
  '
  rm -f "$_mtr_eqrows" "$_mtr_arrows"
}

# median_of FILE -> prints the median of the numbers in FILE (one per line).
median_of() {
  sort -n "$1" | awk '{a[NR]=$1} END{
    n=NR
    if (n%2==1) printf "%.6f\n", a[(n+1)/2]
    else        printf "%.6f\n", (a[n/2]+a[n/2+1])/2
  }'
}

# parse_time_output FILE -> prints "elapsed_seconds rss_bytes" on stdout,
# from a /usr/bin/time -l (macOS) or -v (Linux) report possibly mixed
# with the measured command's own stderr (see the header comment's
# "Build check" section).
parse_time_output() {
  _pto_file=$1
  if grep -q "maximum resident set size" "$_pto_file"; then
    _pto_elapsed=$(awk '$2=="real"{print $1; exit}' "$_pto_file")
    _pto_rss=$(awk '/maximum resident set size/{print $1; exit}' "$_pto_file")
  elif grep -q "Maximum resident set size" "$_pto_file"; then
    _pto_elapsed=$(awk -F': ' '/Elapsed \(wall clock\) time/{print $2; exit}' "$_pto_file" | \
      awk -F: '{
        if (NF==3)      printf "%.3f", $1*3600+$2*60+$3
        else if (NF==2) printf "%.3f", $1*60+$2
        else            printf "%.3f", $1
      }')
    _pto_rss_kb=$(awk -F': ' '/Maximum resident set size/{print $2; exit}' "$_pto_file")
    _pto_rss=$(awk -v kb="$_pto_rss_kb" 'BEGIN{printf "%.0f", kb*1024}')
  else
    echo "parse_time_output: no recognizable RSS line in $_pto_file" >&2
    exit 1
  fi
  printf '%s %s\n' "$_pto_elapsed" "$_pto_rss"
}

# ---------------------------------------------------------------------
# Energy table + MFP table
# ---------------------------------------------------------------------

echo "" >&2
echo "Degree: $TQ_DEGREE, Threads: $TQ_THREADS, Replicates (K): $K," \
     "Seed base: $TQ_SEED_BASE, Configs: $TQ_CONFIGS" >&2

echo ""
echo "# TOA Quadrature Scenario Comparison"
echo ""
echo "Degree: $TQ_DEGREE, Threads: $TQ_THREADS, Replicates (K): $K, Seed base: $TQ_SEED_BASE"
echo ""
echo "## Energy table"
echo ""
echo "E_s = a station's total trace energy (X+Y+Z, every bin, whole run)." \
     "r_s = mean(E_s, area) / mean(E_s, equal)."
echo ""
echo "| config | stations | min r | median r | max r | stations with \\|r-1\\| > 3 SE |"
echo "|---|---:|---:|---:|---:|---:|"

: > "$work/mfp_table_body"
: > "$work/energy_detail_body"

cfg_index=0
for cfg in $TQ_CONFIGS; do
  fields=$(tq_config_fields "$cfg")
  set -- $fields
  scenario=$1; variant=$2

  model=$(model_for_variant "$scenario" "$variant")
  phonons=$(phonons_for_scenario "$scenario")

  echo "Config: $cfg (scenario=$scenario variant=$variant, K=$K, index=$cfg_index)" >&2

  : > "$work/eq.list"
  : > "$work/ar.list"
  first_eq_stdout=""
  first_ar_stdout=""

  k=1
  while [ "$k" -le "$K" ]; do
    seed=$(seed_for "$cfg_index" "$k")
    outdir=$(mktemp -d "$work/run.XXXXXX")
    echo "  Running equal k=$k seed=$seed..." >&2
    run_case "$scenario" "$outdir" "$phonons" "$model" "$seed" equal "$TQ_DEGREE"
    tsv="$outdir/energy.tsv"
    extract_total_energy "$outdir/seis_traces_asc.dat" "$tsv"
    echo "$tsv" >> "$work/eq.list"
    if [ "$k" -eq 1 ]; then first_eq_stdout="$outdir/stdout.txt"; fi
    k=$((k + 1))
  done

  kend=$((2 * K))
  k=$((K + 1))
  while [ "$k" -le "$kend" ]; do
    seed=$(seed_for "$cfg_index" "$k")
    outdir=$(mktemp -d "$work/run.XXXXXX")
    echo "  Running area k=$k seed=$seed..." >&2
    run_case "$scenario" "$outdir" "$phonons" "$model" "$seed" area "$TQ_DEGREE"
    tsv="$outdir/energy.tsv"
    extract_total_energy "$outdir/seis_traces_asc.dat" "$tsv"
    echo "$tsv" >> "$work/ar.list"
    if [ "$k" -eq $((K + 1)) ]; then first_ar_stdout="$outdir/stdout.txt"; fi
    k=$((k + 1))
  done

  echo "  Computing energy-table statistic..." >&2
  result=$(energy_table_row "$work/eq.list" "$work/ar.list" "$K" "$cfg" "$work/.energy_detail_rows")
  set -- $result
  n_stations=$1; min_r=$2; med_r=$3; max_r=$4; flagged=$5

  printf '| %s | %s | %s | %s | %s | %s |\n' \
    "$cfg" "$n_stations" "$min_r" "$med_r" "$max_r" "$flagged"
  cat "$work/.energy_detail_rows" >> "$work/energy_detail_body"
  rm -f "$work/.energy_detail_rows"

  echo "  Computing MFP-table rows..." >&2
  mfp_table_rows "$first_eq_stdout" "$first_ar_stdout" "$cfg" >> "$work/mfp_table_body"

  # Free this config's run directories before the next one.
  rm -rf "$work"/run.* "$work/eq.list" "$work/ar.list"
  cfg_index=$((cfg_index + 1))
done

echo ""
echo "## Energy table: per-station detail"
echo ""
echo "One row per station (sorted by station number), plus a final \"all\"" \
     "row: R = mean(sum-over-stations E, area) / mean(sum-over-stations E, equal)," \
     "the same formula applied to each replicate's whole-scenario total" \
     "instead of one station. z = (r_s - 1)/SE(r_s) (or (R-1)/SE(R) for \"all\")."
echo ""
echo "| config | station | mean E equal | mean E area | r_s | SE(r_s) | z |"
echo "|---|---|---:|---:|---:|---:|---:|"
cat "$work/energy_detail_body"

echo ""
echo "## MFP table"
echo ""
echo "From the first run of each mode. Isotropic rows (no ay/az/mode) show \`-\`." \
     "%change = (area/equal - 1) * 100."
echo ""
echo "| config | a | ay | az | mode | MFP(P) equal | MFP(P) area | MFP(P) %change | MFP(S) equal | MFP(S) area | MFP(S) %change |"
echo "|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|"
cat "$work/mfp_table_body"

# ---------------------------------------------------------------------
# Degree-9 build check
# ---------------------------------------------------------------------

if [ "$TQ_BUILD_CHECK" = "1" ]; then
  BUILD_DEGREE=9
  BUILD_REPS=3

  echo "" >&2
  echo "Build check: Lop Nor iso + axi10, degree $BUILD_DEGREE, threads $TQ_THREADS," \
       "$BUILD_REPS runs per cell" >&2

  TIME_WRAPPER="$work/time_wrapper.sh"
  cat > "$TIME_WRAPPER" <<EOF
#!/bin/sh
set -eu
case "\$(uname -s)" in
  Darwin) exec /usr/bin/time -l "$ROOT/main" "\$@" 2>time_output.txt ;;
  *)      exec /usr/bin/time -v "$ROOT/main" "\$@" 2>time_output.txt ;;
esac
EOF
  chmod +x "$TIME_WRAPPER"

  model_iso=$(scenario_model lopnor "$TQ_EPS" iso)
  model_axi10=$(scenario_model lopnor "$TQ_EPS" axi 10)

  echo ""
  echo "## Degree-9 build check (--num-phonons=0)"
  echo ""
  echo "Lop Nor, degree $BUILD_DEGREE, threads $TQ_THREADS, median of $BUILD_REPS runs per cell." \
      "RSS is the main process's own maximum resident set size."
  echo ""
  echo "| variant | median s equal | median s area | delta s (area-equal) | median RSS MB equal | median RSS MB area | delta RSS MB (area-equal) |"
  echo "|---|---:|---:|---:|---:|---:|---:|"

  for variant_label in iso axi10; do
    if [ "$variant_label" = "iso" ]; then model="$model_iso"; else model="$model_axi10"; fi

    for mode in equal area; do
      : > "$work/bc.times"
      : > "$work/bc.rss"
      k=1
      while [ "$k" -le "$BUILD_REPS" ]; do
        seed=$(seed_for 0 "$k")
        outdir=$(mktemp -d "$work/bcrun.XXXXXX")
        echo "  Build check: variant=$variant_label mode=$mode k=$k seed=$seed..." >&2
        run_case lopnor "$outdir" 0 "$model" "$seed" "$mode" "$BUILD_DEGREE" "$TIME_WRAPPER"
        parsed=$(parse_time_output "$outdir/time_output.txt")
        set -- $parsed
        echo "$1" >> "$work/bc.times"
        echo "$2" >> "$work/bc.rss"
        k=$((k + 1))
      done
      med_s=$(median_of "$work/bc.times")
      med_rss=$(median_of "$work/bc.rss")
      if [ "$mode" = "equal" ]; then
        eq_s=$med_s; eq_rss=$med_rss
      else
        ar_s=$med_s; ar_rss=$med_rss
      fi
      rm -rf "$work"/bcrun.* "$work/bc.times" "$work/bc.rss"
    done

    eq_s_fmt=$(awk -v v="$eq_s" 'BEGIN{printf "%.3f", v}')
    ar_s_fmt=$(awk -v v="$ar_s" 'BEGIN{printf "%.3f", v}')
    ds_fmt=$(awk -v a="$ar_s" -v e="$eq_s" 'BEGIN{printf "%+.3f", a-e}')
    eq_rss_mb=$(awk -v b="$eq_rss" 'BEGIN{printf "%.1f", b/1e6}')
    ar_rss_mb=$(awk -v b="$ar_rss" 'BEGIN{printf "%.1f", b/1e6}')
    drss_mb=$(awk -v a="$ar_rss" -v e="$eq_rss" 'BEGIN{printf "%+.1f", (a-e)/1e6}')

    printf '| %s | %s | %s | %s | %s | %s | %s |\n' \
      "$variant_label" "$eq_s_fmt" "$ar_s_fmt" "$ds_fmt" "$eq_rss_mb" "$ar_rss_mb" "$drss_mb"
  done
fi

echo "" >&2
echo "Comparison complete" >&2

exit 0
