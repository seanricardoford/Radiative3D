#!/bin/sh
set -eu

# tests/test_toa_quadrature.sh
#
# Integration tests for the --toa-quadrature command-line option and its
# area-mode provenance output. Validates that:
#  - --toa-quadrature=area|equal (both "=value" and "value" forms) are
#    accepted; a bad or missing value is rejected with the standard
#    command-line error message
#  - Explicit --toa-quadrature=equal is byte-identical to the default
#    (no flag)
#  - Area mode prints exactly the three provenance lines the brief
#    specifies, in the right places; default (equal) mode prints no
#    "quadrature" text at all
#  - Area mode actually changes the traces relative to equal mode
#  - Area mode is deterministic: reproducible at THREADS=1, and
#    THREADS=12 agrees with THREADS=1 after address-masking -- for an
#    isotropic layout and for the anisotropic lopnor_21/crustpinch_42
#    layouts used by tests/test_aniso_scenarios.sh
#  - Area mode preserves the five-value/seven-value isotropic-layout
#    equivalence tests/test_aniso_inputs.sh checks for equal mode
#  - --event-test with --toa-quadrature=area neither crashes nor changes
#    the displayed radiation pattern (only the parameter echo above it
#    gains the area line)

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

# shellcheck source=tests/aniso_lib.sh
. "$ROOT/tests/aniso_lib.sh"

make >/dev/null

work=$(mktemp -d "${TMPDIR:-/tmp}/radiative3d-toa-quadrature.XXXXXX")
cleanup() {
  rm -rf "$work"
}
trap cleanup EXIT HUP INT TERM

##
## Options
##

echo "Testing --toa-quadrature option parsing..."

if ! ./main --toa-quadrature=area --help >"$work/.opt-area-eq.out" 2>&1; then
  echo "FAIL: --toa-quadrature=area --help should exit 0" >&2
  cat "$work/.opt-area-eq.out" >&2
  exit 1
fi

if ! ./main --toa-quadrature equal --help >"$work/.opt-equal-sp.out" 2>&1; then
  echo "FAIL: --toa-quadrature equal --help should exit 0" >&2
  cat "$work/.opt-equal-sp.out" >&2
  exit 1
fi

if ./main --toa-quadrature=bogus --help >"$work/.opt-bogus.out" 2>&1; then
  echo "FAIL: --toa-quadrature=bogus --help should exit nonzero" >&2
  exit 1
fi
if ! grep -q "Error processing command-line option" "$work/.opt-bogus.out"; then
  echo "FAIL: --toa-quadrature=bogus --help output missing 'Error processing command-line option'" >&2
  cat "$work/.opt-bogus.out" >&2
  exit 1
fi

if ./main --toa-quadrature --help >"$work/.opt-noval.out" 2>&1; then
  echo "FAIL: --toa-quadrature --help (no value) should exit nonzero" >&2
  exit 1
fi
if ! grep -q "Error processing command-line option" "$work/.opt-noval.out"; then
  echo "FAIL: --toa-quadrature --help (no value) output missing 'Error processing command-line option'" >&2
  cat "$work/.opt-noval.out" >&2
  exit 1
fi

echo "  option parsing: PASS"

##
## Explicit equal equals the default; provenance; effect
##

echo "Testing explicit equal vs default, provenance, and effect..."

eq_default="$work/eq-default"
TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$eq_default" 20000

eq_equal="$work/eq-equal"
EXTRA_ARGS="--toa-quadrature=equal" TOA_DEGREE=5 \
  "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$eq_equal" 20000

if ! compare_runs "$eq_default" "$eq_equal" strict; then
  echo "FAIL: explicit --toa-quadrature=equal differs from the default (no flag)" >&2
  exit 1
fi
echo "  explicit equal equals default: PASS"

eq_area="$work/eq-area"
EXTRA_ARGS="--toa-quadrature=area" TOA_DEGREE=5 \
  "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$eq_area" 20000

# Provenance: the three area-only lines, verbatim, each on their own line.
if ! grep -qxF "TOA Quadrature: area" "$eq_area/stdout.txt"; then
  echo "FAIL: area run stdout.txt missing 'TOA Quadrature: area'" >&2
  exit 1
fi
if ! grep -qxF "|          (Area-weighted take-off-angle quadrature.)" "$eq_area/stdout.txt"; then
  echo "FAIL: area run stdout.txt missing '|          (Area-weighted take-off-angle quadrature.)'" >&2
  exit 1
fi
if ! grep -qxF "#  TOA quadrature: area" "$eq_area/stdout.txt"; then
  echo "FAIL: area run stdout.txt missing '#  TOA quadrature: area'" >&2
  exit 1
fi
echo "  area-mode provenance lines present: PASS"

# Default (no-flag) run must not mention "quadrature" anywhere.
qcount=$(grep -ci quadrature "$eq_default/stdout.txt" || true)
if [ "$qcount" -ne 0 ]; then
  echo "FAIL: default run stdout.txt contains $qcount 'quadrature' mention(s), expected 0" >&2
  exit 1
fi
echo "  default run has no quadrature mentions: PASS"

# Effect: area-mode traces differ from equal/default-mode traces.
sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$eq_default/seis_traces_asc.dat" \
  > "$eq_default/.seis_traces_asc.normalized"
sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$eq_area/seis_traces_asc.dat" \
  > "$eq_area/.seis_traces_asc.normalized"
if cmp -s "$eq_default/.seis_traces_asc.normalized" "$eq_area/.seis_traces_asc.normalized"; then
  echo "FAIL: area-mode and default traces are identical (expected a difference)" >&2
  exit 1
fi
echo "  area-mode traces differ from default: PASS"

##
## Determinism in area mode
##

# Anisotropic layouts, copied verbatim from tests/test_aniso_scenarios.sh.
lopnor_21="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,2.0,2.0,0.5,0.3,1000,0.8,0.01,1.4,0.7,0.35,0.5,300"
crustpinch_42="0.8,0.01,0.2,0.2,0.2,0.2,200,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.6,0.3,0.15,0.4,1500,0.8,0.01,0.2,0.2,0.2,0.5,900,2.0,30.0,5.0,.3666667,.4736842,1,1"

##
## check_area_determinism LABEL SCENARIO MODEL_ARGS
##
## Two THREADS=1 area-mode runs must match strictly, and a THREADS=12
## area-mode run must match THREADS=1 after address-masking.
##
check_area_determinism() {
  _label=$1
  _scenario=$2
  _model_args=$3

  echo "Testing $_label area-mode determinism..."

  _d1="$work/det-$_label-d1"
  EXTRA_ARGS="--toa-quadrature=area" TOA_DEGREE=5 \
    "$ROOT/tests/run_aniso_case.sh" "$_scenario" 1 "$_d1" 20000 "$_model_args"

  _d2="$work/det-$_label-d2"
  EXTRA_ARGS="--toa-quadrature=area" TOA_DEGREE=5 \
    "$ROOT/tests/run_aniso_case.sh" "$_scenario" 1 "$_d2" 20000 "$_model_args"

  if ! compare_runs "$_d1" "$_d2" strict; then
    echo "FAIL: $_label area-mode THREADS=1 runs differ (not reproducible)" >&2
    exit 1
  fi
  echo "  $_label area-mode THREADS=1 reproducible: PASS"

  _d12="$work/det-$_label-d12"
  EXTRA_ARGS="--toa-quadrature=area" TOA_DEGREE=5 \
    "$ROOT/tests/run_aniso_case.sh" "$_scenario" 12 "$_d12" 20000 "$_model_args"

  sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$_d1/seis_traces_asc.dat" \
    > "$_d1/.seis_traces_asc.normalized"
  sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$_d12/seis_traces_asc.dat" \
    > "$_d12/.seis_traces_asc.normalized"
  if ! cmp "$_d1/.seis_traces_asc.normalized" "$_d12/.seis_traces_asc.normalized"; then
    echo "FAIL: $_label area-mode THREADS=12 traces differ from THREADS=1 after address masking" >&2
    exit 1
  fi
  echo "  $_label area-mode THREADS=12 matches THREADS=1 (address-masked): PASS"
}

check_area_determinism lopnor-default lopnor ""
check_area_determinism lopnor-21 lopnor "$lopnor_21"
check_area_determinism crustpinch-42 crustpinch "$crustpinch_42"

##
## Seven- vs five-value equivalence in area mode
##

echo "Testing seven- vs five-value equivalence in area mode..."

# Lop Nor 15-value vs 21-value isotropic layouts, copied verbatim from
# tests/test_aniso_inputs.sh.
args_5_lopnor_15="0.8,0.01,0.5,0.2,50,0.8,0.01,0.5,0.3,1000,0.8,0.01,0.7,0.5,300"
args_7_lopnor_21="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,0.5,0.5,0.5,0.3,1000,0.8,0.01,0.7,0.7,0.7,0.5,300"

sf_dir5="$work/sf-lopnor-15"
sf_dir7="$work/sf-lopnor-21"
EXTRA_ARGS="--toa-quadrature=area" TOA_DEGREE=4 \
  "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$sf_dir5" 20000 "$args_5_lopnor_15"
EXTRA_ARGS="--toa-quadrature=area" TOA_DEGREE=4 \
  "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$sf_dir7" 20000 "$args_7_lopnor_21"

if ! compare_runs "$sf_dir5" "$sf_dir7" noargs; then
  echo "FAIL: area-mode 5-value and 7-value (isotropic) lopnor 15->21 output differ" >&2
  exit 1
fi
echo "  area-mode lopnor 15-value vs 21-value: PASS"

##
## Event test
##

echo "Testing --event-test with area quadrature..."

if ! ./main --event-test --source=SDR,125,40,90 \
    >"$work/.event-plain.out" 2>&1; then
  echo "FAIL: --event-test --source=SDR,125,40,90 should exit 0" >&2
  cat "$work/.event-plain.out" >&2
  exit 1
fi

if ! ./main --event-test --source=SDR,125,40,90 --toa-quadrature=area \
    >"$work/.event-area.out" 2>&1; then
  echo "FAIL: --event-test --source=SDR,125,40,90 --toa-quadrature=area should exit 0" >&2
  cat "$work/.event-area.out" >&2
  exit 1
fi
echo "  --event-test exits 0 in both modes: PASS"

# Stdout must be identical from the @@ __EVENT_SOURCE_TEST__ line onward.
awk '/@@ __EVENT_SOURCE_TEST__/{f=1} f' "$work/.event-plain.out" \
  > "$work/.event-plain.tail"
awk '/@@ __EVENT_SOURCE_TEST__/{f=1} f' "$work/.event-area.out" \
  > "$work/.event-area.tail"

if ! cmp "$work/.event-plain.tail" "$work/.event-area.tail"; then
  echo "FAIL: --event-test output from @@ __EVENT_SOURCE_TEST__ onward differs between modes" >&2
  exit 1
fi
echo "  --event-test output identical from @@ __EVENT_SOURCE_TEST__ onward: PASS"

# The parameter echo above that line must gain exactly the area line --
# removing it from the area-mode head must leave the plain-mode head.
awk '/@@ __EVENT_SOURCE_TEST__/{exit} {print}' "$work/.event-plain.out" \
  > "$work/.event-plain.head"
awk '/@@ __EVENT_SOURCE_TEST__/{exit} {print}' "$work/.event-area.out" \
  > "$work/.event-area.head"
grep -vxF "TOA Quadrature: area" "$work/.event-area.head" \
  > "$work/.event-area.head.stripped"

if ! cmp "$work/.event-plain.head" "$work/.event-area.head.stripped"; then
  echo "FAIL: area --event-test parameter echo differs from plain by more than the area line" >&2
  exit 1
fi
echo "  --event-test parameter echo gains exactly the area line: PASS"

echo ""
echo "All test_toa_quadrature checks passed"
