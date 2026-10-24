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
#   - axisymmetric (ax=ay=r*a, az=a) and triaxial (ax=r*a, ay=sqrt(r)*a,
#     az=a) runs at aspect ratios r in {2, 4, 10}.
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
# Environment:
#   BASELINE_EXE             Required. Pre-anisotropy executable for the
#                             isotropic no-slowdown check.
#   BENCH_DEGREE              --toa-degree for all runs.      Default: 9
#   BENCH_LOPNOR_PHONONS       --num-phonons for lopnor runs.  Default: 1000000
#   BENCH_CRUSTPINCH_PHONONS   --num-phonons for crustpinch.   Default: 200000
#   BENCH_THREADS              --threads for all runs, or "auto" to use the
#                             logical CPU count (sysctl -n hw.ncpu on macOS,
#                             nproc elsewhere).                Default: auto
#
# Output: a Markdown table (scenario, eps, mode, r, seconds, slowdown vs.
# isotropic, gate) and the isotropic current/baseline ratio lines go to
# stdout. Progress messages go to stderr.
#
# Exit status: nonzero if BASELINE_EXE is missing/unusable, if any
# simulation run fails, or if either scenario's isotropic current/baseline
# ratio exceeds 1.05.

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
. "$ROOT/tests/aniso_lib.sh"   # region5, region7, mode_mults,
                               # lopnor_model, crustpinch_model,
                               # scenario_model (and their mul/sqrt_of
                               # helpers) -- shared with
                               # tests/compare_axi_general.sh.

BENCH_DEGREE=${BENCH_DEGREE:-9}
BENCH_LOPNOR_PHONONS=${BENCH_LOPNOR_PHONONS:-1000000}
BENCH_CRUSTPINCH_PHONONS=${BENCH_CRUSTPINCH_PHONONS:-200000}
BENCH_THREADS=${BENCH_THREADS:-auto}

if [ "$BENCH_THREADS" = "auto" ]; then
  if [ "$(uname -s)" = "Darwin" ]; then
    BENCH_THREADS=$(sysctl -n hw.ncpu 2>/dev/null || echo 1)
  else
    BENCH_THREADS=$(nproc 2>/dev/null || echo 1)
  fi
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
MODE_LIST="axi tri"
R_LIST="2 4 10"

GATE_FAIL=0

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
#
# region5, region7, mode_mults, lopnor_model, crustpinch_model, and
# scenario_model (the --model-args layout builders) now live in
# tests/aniso_lib.sh, sourced above -- shared with
# tests/compare_axi_general.sh, which also needs mode "gen".

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
         "model=$model_args); see $outdir/stderr.txt and" \
         "$outdir/harness-stderr.txt" >&2
    exit 1
  fi
  end=$(date +%s.%N)
  awk "BEGIN{printf \"%.3f\", $end - $start}"
}

# median3 a b c -> median
median3() {
  printf '%s\n%s\n%s\n' "$1" "$2" "$3" | sort -n | sed -n '2p'
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
    GATE_FAIL=1
  fi
}

# ---------------------------------------------------------------------
# Main benchmark table
# ---------------------------------------------------------------------

benchmark_scenario() {
  scenario=$1; phonons=$2

  for eps in $EPS_LIST; do
    echo "Measuring $scenario eps=$eps isotropic..." >&2
    iso_model=$(scenario_model "$scenario" "$eps" iso)
    t_iso=$(time_run "$CURRENT_EXE" "$scenario" "$phonons" "$iso_model")
    printf '| %s | %s | iso | - | %.2f | 1.00x | ok |\n' \
      "$scenario" "$eps" "$t_iso"

    for mode in $MODE_LIST; do
      for r in $R_LIST; do
        echo "Measuring $scenario eps=$eps mode=$mode r=$r..." >&2
        model=$(scenario_model "$scenario" "$eps" "$mode" "$r")
        t=$(time_run "$CURRENT_EXE" "$scenario" "$phonons" "$model")
        slowdown=$(awk "BEGIN{printf \"%.3f\", $t/$t_iso}")
        if awk "BEGIN{exit !($slowdown > 2.0)}"; then
          gate="GATE"
          GATE_FAIL=1
        else
          gate="ok"
        fi
        printf '| %s | %s | %s | %s | %.2f | %sx | %s |\n' \
          "$scenario" "$eps" "$mode" "$r" "$t" "$slowdown" "$gate"
      done
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
echo "Degree: $BENCH_DEGREE, Threads: $BENCH_THREADS"
echo ""
echo "| Scenario | eps | Mode | r | Seconds | Slowdown vs isotropic | Gate |"
echo "|---|---|---|---|---:|---:|---|"

benchmark_scenario lopnor "$BENCH_LOPNOR_PHONONS"
benchmark_scenario crustpinch "$BENCH_CRUSTPINCH_PHONONS"

echo "" >&2
echo "Benchmark complete" >&2

exit "$GATE_FAIL"
