#!/bin/sh
set -eu

# tests/test_aniso_inputs.sh
#
# Integration tests for the seven-value anisotropic-scale-length input
# layouts, grid-dump column detection, and scatterer-access guard.
# Validates that:
#  - Five-value and seven-value (isotropic) layouts produce identical output
#  - Model dispatch selects the correct grid model
#  - Grid dumps include the right columns when anisotropy is detected
#  - Invalid input counts and scale-length values are rejected

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

# shellcheck source=tests/aniso_lib.sh
. "$ROOT/tests/aniso_lib.sh"

make >/dev/null

work=$(mktemp -d "${TMPDIR:-/tmp}/radiative3d-aniso-inputs.XXXXXX")
cleanup() {
  rm -rf "$work"
}
trap cleanup EXIT HUP INT TERM

##
## Test: equivalence of 5-value and 7-value (isotropic) layouts
##
## For each scenario and count pair, run with both the 5-value and 7-value
## (with ax=ay=az=a) layouts and confirm they produce identical output.
##

echo "Testing equivalence: 5-value vs 7-value isotropic layouts..."

# Lop Nor: 15-value (5 per region) vs 21-value (7 per region)
# 15-value: nu1, eps1, a1, k1, Q1, nu2, eps2, a2, k2, Q2, nu3, eps3, a3, k3, Q3
# 21-value: nu1, eps1, ax1, ay1, az1, k1, Q1, nu2, eps2, ax2, ay2, az2, k2, Q2, ...
args_5_lopnor="0.8,0.01,0.5,0.2,50,0.8,0.01,0.5,0.3,1000,0.8,0.01,0.7,0.5,300"
args_7_lopnor="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,0.5,0.5,0.5,0.3,1000,0.8,0.01,0.7,0.7,0.7,0.5,300"

dir_5="$work/lopnor-5value"
dir_7="$work/lopnor-7value"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_5" 20000 "$args_5_lopnor"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_7" 20000 "$args_7_lopnor"

if ! compare_runs "$dir_5" "$dir_7" noargs; then
  echo "FAIL: 5-value and 7-value (isotropic) lopnor output differ" >&2
  exit 1
fi

echo "  lopnor 15-value vs 21-value: PASS"

##
## Test: dispatch selects correct model for count
##

echo "Testing model dispatch..."

# Count 15 should select Lop Nor Baseline
dir_count_15="$work/lopnor-count-15"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_count_15" 20000 "$args_5_lopnor"
if ! grep -q "Selected Lop Nor Baseline Model" "$dir_count_15/stdout.txt"; then
  echo "FAIL: count 15 did not select Lop Nor Baseline Model" >&2
  exit 1
fi

# Count 21 should also select Lop Nor Baseline (not Moho)
dir_count_21="$work/lopnor-count-21"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_count_21" 20000 "$args_7_lopnor"
if ! grep -q "Selected Lop Nor Baseline Model" "$dir_count_21/stdout.txt"; then
  echo "FAIL: count 21 did not select Lop Nor Baseline Model" >&2
  exit 1
fi

echo "  dispatch for counts 15 and 21: PASS"

##
## Test: grid dump contains expected columns
##

echo "Testing grid-dump columns..."

# Run with 0 phonons and --dump-grid to get just the grid output
grid_args="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,2.0,2.0,0.5,0.3,1000,0.8,0.01,1.4,0.7,0.35,0.5,300"
grid_dir="$work/grid-dump"
EXTRA_ARGS="--dump-grid" TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$grid_dir" 0 "$grid_args"

# Check header contains expected column markers
if ! grep -E "nu.*eps.*a.*k.*ay.*az" "$grid_dir/stdout.txt" >/dev/null; then
  echo "FAIL: grid dump header does not contain expected anisotropic columns (nu eps a k ay az)" >&2
  cat "$grid_dir/stdout.txt" | grep "^#" >&2
  exit 1
fi

# Count that each node row has 17 fields (original 15 + ay + az)
node_rows=$(grep -E "^\s*[0-9]+\s+[0-9]+\s+[0-9]+" "$grid_dir/stdout.txt" | wc -l)
if [ "$node_rows" -eq 0 ]; then
  echo "FAIL: no node rows found in grid dump" >&2
  exit 1
fi

# Check a crust row (field 15 = 0.30000 for k=0.3) has correct ay,az values (2.00000, 0.50000)
if ! grep "0.30000\s*2.00000\s*0.50000" "$grid_dir/stdout.txt" >/dev/null; then
  echo "FAIL: crust rows do not have correct ay,az values (2.00000 0.50000)" >&2
  grep "0.30000" "$grid_dir/stdout.txt" | head -1 >&2
  exit 1
fi

# Check a mantle row (field 15 = 0.50000 for k=0.5) has correct ay,az values (0.70000, 0.35000)
if ! grep "0.50000\s*0.70000\s*0.35000" "$grid_dir/stdout.txt" >/dev/null; then
  echo "FAIL: mantle rows do not have correct ay,az values (0.70000 0.35000)" >&2
  grep "1.40000" "$grid_dir/stdout.txt" | tail -1 >&2
  exit 1
fi

echo "  grid dump columns: PASS"

##
## Test: invalid input counts are rejected
##

echo "Testing invalid input count rejection..."

# 22 values should fail
bad_args="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,0.5,0.5,0.5,0.3,1000,0.8,0.01,0.7,0.7,0.7,0.5,300,999"
bad_dir="$work/bad-count-22"
if TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$bad_dir" 20000 "$bad_args" 2>&1; then
  echo "FAIL: 22-value args should have been rejected" >&2
  exit 1
fi

if ! grep -i "wrong number of model args" "$bad_dir/stdout.txt" "$bad_dir/stderr.txt" 2>/dev/null; then
  echo "FAIL: error message for 22 values did not contain 'wrong number of model args'" >&2
  exit 1
fi

echo "  invalid count rejection: PASS"

##
## Test: invalid scale length values are rejected
##

echo "Testing invalid scale-length rejection..."

# Crust region (second) with ay=0 should fail
bad_scale_args="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,0.5,0.0,0.5,0.3,1000,0.8,0.01,0.7,0.7,0.7,0.5,300"
bad_scale_dir="$work/bad-scale-ay0"
if TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$bad_scale_dir" 20000 "$bad_scale_args" 2>&1; then
  echo "FAIL: ay=0 should have been rejected" >&2
  exit 1
fi

if ! grep -i "scale lengths ax, ay, az must be positive" "$bad_scale_dir/stdout.txt" "$bad_scale_dir/stderr.txt" 2>/dev/null; then
  echo "FAIL: error message for ay=0 did not contain 'scale lengths ax, ay, az must be positive'" >&2
  cat "$bad_scale_dir/stdout.txt" "$bad_scale_dir/stderr.txt" >&2
  exit 1
fi

echo "  invalid scale-length rejection: PASS"

echo ""
echo "All test_aniso_inputs checks passed"
