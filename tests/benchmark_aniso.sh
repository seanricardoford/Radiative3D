#!/bin/sh
set -eu
#
# tests/benchmark_aniso.sh
#
# Manual performance benchmark for anisotropic scale-length heterogeneity.
# Not named test_*, so the full test-suite check skips it; run it directly.
#
# For each of the "lopnor" and "crustpinch" scenarios, and for epsilon in
# {do-script value (0.01), 0.05}, times:
#   - the isotropic run (do-script "a" values, five-value layout), and
#   - axisymmetric (ax=ay=r*a, az=a), near-axisymmetric-general (ax=r*a,
#     ay=ax*(1+1e-5), az=a -- the same medium as the axisymmetric row, to
#     within 1e-5) and triaxial (ax=r*a, ay=sqrt(r)*a, az=a) runs at
#     aspect ratios r in {2, 4, 10}.
#
# Anisotropy is applied to the crust and mantle regions:
#   - lopnor: crust, mantle (sedi stays isotropic).
#   - crustpinch: crust, pinched crust, mantle (sedi and moho stay
#     isotropic). This mirrors the natural reading of the model's five
#     regions used elsewhere in this task's validation (see
#     docs/anisotropy-performance.md).
#
# Radiative3D requires every region of a given invocation to use the same
# layout width (all five-value, or all seven-value -- see doc/MANUAL.md,
# "Seven-value (anisotropic) model-args layouts"). So whenever a run
# contains any anisotropic region, *every* region in that run is written in
# seven-value form; regions that are meant to stay isotropic simply get
# ax=ay=az=a (byte-identical to the five-value form, per MANUAL.md).
#
# Each (scenario, eps, mode, r) row also gets a "build s" column: the same
# run's wall-clock time with --num-phonons=0, i.e. model construction
# (table build) alone, with the simulation loop doing no work. This
# isolates the construction-cost gap between samplers (see
# docs/superpowers/specs/2026-10-17-axisymmetric-sampler-performance-
# design.md's "Measured starting point") from phonon-tracing time. The
# build column is a single measurement, not a repeat median (only the
# "total s" column is -- see BENCH_REPEATS below).
#
# Environment:
#   BASELINE_EXE              Required. Pre-anisotropy executable for the
#                             isotropic no-slowdown check.
#   BENCH_DEGREE              --toa-degree for all runs.      Default: 9
#   BENCH_LOPNOR_PHONONS       --num-phonons for lopnor runs.  Default: 1000000
#   BENCH_CRUSTPINCH_PHONONS   --num-phonons for crustpinch.   Default: 200000
#   BENCH_THREADS              --threads for all runs, or "auto" for the
#                             logical CPU count (aniso_lib.sh's
#                             auto_threads).                   Default: auto
#   BENCH_REPEATS              Repeats of each (scenario, eps, mode, r)
#                             run; the table's "total s" column is the
#                             median of these. The isotropic current/
#                             baseline check below is unaffected -- it
#                             always uses a fixed 3 repeats per side.
#                                                               Default: 3
#
# Output: a Markdown table (scenario, eps, mode, r, build s, total s,
# slowdown vs. isotropic, gate vs. isotropic, axi/gen same-medium ratio,
# gate) and the isotropic current/baseline ratio lines go to stdout.
# Progress messages go to stderr.
#
# Exit status: nonzero if BASELINE_EXE is missing/unusable, if any
# simulation run fails, if either scenario's isotropic current/baseline
# ratio exceeds 1.05, or if any axi row's same-medium axi/gen ratio
# (design decision 5) exceeds 1.05 (reported as "SLOW"). A vs-isotropic
# slowdown over 2x ("GATE") is reported in its own column but does NOT
# affect exit status: per task-7-brief.md Step 2, a GATE row is the
# controller's call (e.g. a fallback such as AUTO choosing GENERAL for
# axisymmetric specs below some M), not something this script should fail
# the build on by itself.

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
. "$ROOT/tests/aniso_lib.sh"   # region5, region7, mode_mults,
                               # lopnor_model, crustpinch_model,
                               # scenario_model, auto_threads,
                               # dump_failure_logs (and their mul/sqrt_of
                               # helpers) -- shared with
                               # tests/compare_axi_general.sh.

BENCH_DEGREE=${BENCH_DEGREE:-9}
BENCH_LOPNOR_PHONONS=${BENCH_LOPNOR_PHONONS:-1000000}
BENCH_CRUSTPINCH_PHONONS=${BENCH_CRUSTPINCH_PHONONS:-200000}
BENCH_THREADS=${BENCH_THREADS:-auto}
BENCH_REPEATS=${BENCH_REPEATS:-3}

case "$BENCH_REPEATS" in
  ''|*[!0-9]*)
    echo "ERROR: BENCH_REPEATS must be a positive integer (got '$BENCH_REPEATS')" >&2
    exit 1 ;;
esac
if [ "$BENCH_REPEATS" -lt 1 ]; then
  echo "ERROR: BENCH_REPEATS must be >= 1 (got '$BENCH_REPEATS')" >&2
  exit 1
fi

if [ "$BENCH_THREADS" = "auto" ]; then
  BENCH_THREADS=$(auto_threads)
fi

if [ -z "${BASELINE_EXE:-}" ] || [ ! -x "$BASELINE_EXE" ]; then
  echo "ERROR: BASELINE_EXE must be set to an executable pre-anisotropy build" >&2
  exit 1
fi

echo "Building current executable..." >&2
make >/dev/null

CURRENT_EXE="$ROOT/main"
if [ ! -x "$CURRENT_EXE" ]; then
  echo "ERROR: build did not produce an executable at $CURRENT_EXE" >&2
  exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

DOSCRIPT_EPS=0.01   # The epsilon value used by do-lopnor.sh and
                    # do-crustpinch.sh for every region, by default.
EPS_LIST="$DOSCRIPT_EPS 0.05"
R_LIST="2 4 10"

FAIL=0   # Set on: a failed isotropic current/baseline ratio check, or a
         # same-medium axi/gen "SLOW" row. NOT set by a vs-isotropic
         # "GATE" row -- see the Exit status note above. A failed
         # simulation run exits immediately from time_run() instead of
         # going through this accumulator.

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
#
# region5, region7, mode_mults, lopnor_model, crustpinch_model,
# scenario_model (the --model-args layout builders), auto_threads, and
# dump_failure_logs now live in tests/aniso_lib.sh, sourced above --
# shared with tests/compare_axi_general.sh, which also needs mode "gen".

# time_run exe scenario phonons model_args  -> seconds (wall clock) on stdout
# Exits the whole script nonzero if the run fails.
time_run() {
  exe=$1; scenario=$2; phonons=$3; model_args=$4
  outdir=$(mktemp -d "$work/run.XXXXXX")
  start=$(date +%s.%N)
  if ! R3D_EXE="$exe" TOA_DEGREE="$BENCH_DEGREE" "$ROOT/tests/run_aniso_case.sh" \
        "$scenario" "$BENCH_THREADS" "$outdir" "$phonons" "$model_args" \
        >"$outdir/harness-stdout.txt" 2>"$outdir/harness-stderr.txt"; then
    echo "ERROR: run failed (exe=$exe scenario=$scenario phonons=$phonons" \
         "model=$model_args)" >&2
    dump_failure_logs "$outdir"
    exit 1
  fi
  end=$(date +%s.%N)
  awk "BEGIN{printf \"%.3f\", $end - $start}"
}

# median3 a b c -> median
median3() {
  printf '%s\n%s\n%s\n' "$1" "$2" "$3" | sort -n | sed -n '2p'
}

# median_n v1 v2 ... -> median of N values (any count >= 1; BENCH_REPEATS
# need not be 3, unlike median3's fixed-3 isotropic check above).
median_n() {
  _mn_n=$#
  _mn_sorted=$(printf '%s\n' "$@" | sort -n)
  _mn_mid=$(( (_mn_n + 1) / 2 ))
  if [ $((_mn_n % 2)) -eq 1 ]; then
    printf '%s\n' "$_mn_sorted" | sed -n "${_mn_mid}p"
  else
    _mn_lo=$(printf '%s\n' "$_mn_sorted" | sed -n "${_mn_mid}p")
    _mn_hi=$(printf '%s\n' "$_mn_sorted" | sed -n "$((_mn_mid + 1))p")
    awk "BEGIN{printf \"%.3f\", ($_mn_lo+$_mn_hi)/2}"
  fi
}

# bench_total exe scenario phonons model_args  -> median of BENCH_REPEATS
# time_run totals (wall-clock seconds, phonons included).
bench_total() {
  _bt_exe=$1; _bt_scenario=$2; _bt_phonons=$3; _bt_model=$4
  _bt_vals=""
  _bt_i=1
  while [ "$_bt_i" -le "$BENCH_REPEATS" ]; do
    _bt_t=$(time_run "$_bt_exe" "$_bt_scenario" "$_bt_phonons" "$_bt_model")
    echo "    repeat $_bt_i/$BENCH_REPEATS: ${_bt_t}s" >&2
    _bt_vals="$_bt_vals $_bt_t"
    _bt_i=$((_bt_i + 1))
  done
  median_n $_bt_vals
}

# ---------------------------------------------------------------------
# Isotropic no-slowdown check (BASELINE_EXE vs. current build)
# ---------------------------------------------------------------------

check_isotropic_ratio() {
  scenario=$1; phonons=$2
  echo "Isotropic regression check: $scenario (phonons=$phonons)" >&2

  iso_model=$(scenario_model "$scenario" "$DOSCRIPT_EPS" iso)

  b1=$(time_run "$BASELINE_EXE" "$scenario" "$phonons" "$iso_model")
  echo "  baseline run 1: ${b1}s" >&2
  b2=$(time_run "$BASELINE_EXE" "$scenario" "$phonons" "$iso_model")
  echo "  baseline run 2: ${b2}s" >&2
  b3=$(time_run "$BASELINE_EXE" "$scenario" "$phonons" "$iso_model")
  echo "  baseline run 3: ${b3}s" >&2

  c1=$(time_run "$CURRENT_EXE" "$scenario" "$phonons" "$iso_model")
  echo "  current run 1: ${c1}s" >&2
  c2=$(time_run "$CURRENT_EXE" "$scenario" "$phonons" "$iso_model")
  echo "  current run 2: ${c2}s" >&2
  c3=$(time_run "$CURRENT_EXE" "$scenario" "$phonons" "$iso_model")
  echo "  current run 3: ${c3}s" >&2

  med_base=$(median3 "$b1" "$b2" "$b3")
  med_curr=$(median3 "$c1" "$c2" "$c3")
  ratio=$(awk "BEGIN{printf \"%.4f\", $med_curr/$med_base}")

  echo "Isotropic ratio ($scenario, current/baseline, median of 3): $ratio"

  if awk "BEGIN{exit !($ratio > 1.05)}"; then
    echo "ERROR: $scenario isotropic ratio $ratio exceeds 1.05 threshold" >&2
    FAIL=1
  fi
}

# ---------------------------------------------------------------------
# Main benchmark table
# ---------------------------------------------------------------------

# print_row scenario eps mode r build total t_iso same_ratio same_gate
#
#   Prints one Markdown table row. Computes the vs-isotropic slowdown
#   (total/t_iso) and its gate here (ok <= 2x, else "GATE" -- report-only,
#   see the Exit status note above). same_ratio/same_gate are passed
#   through as-is: "-"/"-" for every row except axi, whose caller supplies
#   the axi/gen same-medium ratio and its gate (ok <= 1.05, else "SLOW",
#   which the caller has already folded into FAIL -- see decision 5).
print_row() {
  _pr_scenario=$1; _pr_eps=$2; _pr_mode=$3; _pr_r=$4
  _pr_build=$5; _pr_total=$6; _pr_tiso=$7
  _pr_sameratio=$8; _pr_samegate=$9

  _pr_slowdown=$(awk "BEGIN{printf \"%.3f\", $_pr_total/$_pr_tiso}")
  if awk "BEGIN{exit !($_pr_total/$_pr_tiso > 2.0)}"; then
    _pr_gate="GATE"
  else
    _pr_gate="ok"
  fi

  printf '| %s | %s | %s | %s | %.2f | %.2f | %sx | %s | %s | %s |\n' \
    "$_pr_scenario" "$_pr_eps" "$_pr_mode" "$_pr_r" \
    "$_pr_build" "$_pr_total" "$_pr_slowdown" "$_pr_gate" \
    "$_pr_sameratio" "$_pr_samegate"
}

benchmark_scenario() {
  scenario=$1; phonons=$2

  for eps in $EPS_LIST; do
    echo "Measuring $scenario eps=$eps isotropic (build)..." >&2
    iso_model=$(scenario_model "$scenario" "$eps" iso)
    b_iso=$(time_run "$CURRENT_EXE" "$scenario" 0 "$iso_model")
    echo "Measuring $scenario eps=$eps isotropic (total, x$BENCH_REPEATS)..." >&2
    t_iso=$(bench_total "$CURRENT_EXE" "$scenario" "$phonons" "$iso_model")
    print_row "$scenario" "$eps" iso - "$b_iso" "$t_iso" "$t_iso" - -

    for r in $R_LIST; do
      echo "Measuring $scenario eps=$eps mode=axi r=$r (build)..." >&2
      axi_model=$(scenario_model "$scenario" "$eps" axi "$r")
      b_axi=$(time_run "$CURRENT_EXE" "$scenario" 0 "$axi_model")
      echo "Measuring $scenario eps=$eps mode=axi r=$r (total, x$BENCH_REPEATS)..." >&2
      t_axi=$(bench_total "$CURRENT_EXE" "$scenario" "$phonons" "$axi_model")

      echo "Measuring $scenario eps=$eps mode=gen r=$r (build, same medium as axi)..." >&2
      gen_model=$(scenario_model "$scenario" "$eps" gen "$r")
      b_gen=$(time_run "$CURRENT_EXE" "$scenario" 0 "$gen_model")
      echo "Measuring $scenario eps=$eps mode=gen r=$r (total, x$BENCH_REPEATS)..." >&2
      t_gen=$(bench_total "$CURRENT_EXE" "$scenario" "$phonons" "$gen_model")

      same_ratio=$(awk "BEGIN{printf \"%.3f\", $t_axi/$t_gen}")
      if awk "BEGIN{exit !($t_axi/$t_gen > 1.05)}"; then
        same_gate="SLOW"
        echo "ERROR: $scenario eps=$eps r=$r axi/gen ratio $same_ratio exceeds 1.05 threshold" >&2
        FAIL=1
      else
        same_gate="ok"
      fi

      print_row "$scenario" "$eps" axi "$r" "$b_axi" "$t_axi" "$t_iso" \
        "$same_ratio" "$same_gate"
      print_row "$scenario" "$eps" gen "$r" "$b_gen" "$t_gen" "$t_iso" - -

      echo "Measuring $scenario eps=$eps mode=tri r=$r (build)..." >&2
      tri_model=$(scenario_model "$scenario" "$eps" tri "$r")
      b_tri=$(time_run "$CURRENT_EXE" "$scenario" 0 "$tri_model")
      echo "Measuring $scenario eps=$eps mode=tri r=$r (total, x$BENCH_REPEATS)..." >&2
      t_tri=$(bench_total "$CURRENT_EXE" "$scenario" "$phonons" "$tri_model")
      print_row "$scenario" "$eps" tri "$r" "$b_tri" "$t_tri" "$t_iso" - -
    done
  done
}

# ---------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------

check_isotropic_ratio lopnor "$BENCH_LOPNOR_PHONONS"
check_isotropic_ratio crustpinch "$BENCH_CRUSTPINCH_PHONONS"

echo ""
echo "# Anisotropic Scale-Length Performance Benchmark"
echo ""
echo "Degree: $BENCH_DEGREE, Threads: $BENCH_THREADS, Repeats: $BENCH_REPEATS"
echo ""
echo "| Scenario | eps | Mode | r | Build s | Total s | Slowdown vs isotropic | Gate vs isotropic | axi/gen ratio | Gate same-medium |"
echo "|---|---|---|---|---:|---:|---:|---|---:|---|"

benchmark_scenario lopnor "$BENCH_LOPNOR_PHONONS"
benchmark_scenario crustpinch "$BENCH_CRUSTPINCH_PHONONS"

echo "" >&2
echo "Benchmark complete" >&2

exit "$FAIL"
