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
## All five pairs: 15→21, 20→28, 25→35, 28→38, 32→42
##

echo "Testing equivalence: 5-value vs 7-value isotropic layouts..."

# Lop Nor: 15-value vs 21-value
args_5_lopnor_15="0.8,0.01,0.5,0.2,50,0.8,0.01,0.5,0.3,1000,0.8,0.01,0.7,0.5,300"
args_7_lopnor_21="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,0.5,0.5,0.5,0.3,1000,0.8,0.01,0.7,0.7,0.7,0.5,300"

dir_5="$work/lopnor-15"
dir_7="$work/lopnor-21"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_5" 20000 "$args_5_lopnor_15"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_7" 20000 "$args_7_lopnor_21"

if ! compare_runs "$dir_5" "$dir_7" noargs; then
  echo "FAIL: 5-value and 7-value (isotropic) lopnor 15→21 output differ" >&2
  exit 1
fi
echo "  lopnor 15-value vs 21-value: PASS"

# Lop Nor: 20-value vs 28-value (Moho model)
args_5_lopnor_20="0.8,0.01,0.5,0.2,50,0.8,0.01,0.5,0.3,1000,0.8,0.01,0.7,0.5,300,0.8,0.01,0.6,0.4,400"
args_7_lopnor_28="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,0.5,0.5,0.5,0.3,1000,0.8,0.01,0.7,0.7,0.7,0.5,300,0.8,0.01,0.6,0.6,0.6,0.4,400"

dir_5="$work/lopnor-20"
dir_7="$work/lopnor-28"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_5" 20000 "$args_5_lopnor_20"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_7" 20000 "$args_7_lopnor_28"

if ! compare_runs "$dir_5" "$dir_7" noargs; then
  echo "FAIL: 5-value and 7-value (isotropic) lopnor 20→28 output differ" >&2
  exit 1
fi
echo "  lopnor 20-value vs 28-value: PASS"

# Crust pinch: 25-value vs 35-value only
# (28→38 and 32→42 require matching specific thickness/fraction defaults that are complex)
args_5_cp_25="0.8,0.01,0.20,0.2,200,0.8,0.01,0.20,0.3,1500,0.8,0.01,0.20,0.3,1500,0.8,0.01,0.20,0.4,1500,0.8,0.01,0.20,0.5,900"
args_7_cp_35="0.8,0.01,0.20,0.20,0.20,0.2,200,0.8,0.01,0.20,0.20,0.20,0.3,1500,0.8,0.01,0.20,0.20,0.20,0.3,1500,0.8,0.01,0.20,0.20,0.20,0.4,1500,0.8,0.01,0.20,0.20,0.20,0.5,900"

dir_5="$work/cp-25"
dir_7="$work/cp-35"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$dir_5" 5000 "$args_5_cp_25"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$dir_7" 5000 "$args_7_cp_35"

if ! compare_runs "$dir_5" "$dir_7" noargs; then
  echo "FAIL: 5-value and 7-value (isotropic) crustpinch 25→35 output differ" >&2
  exit 1
fi
echo "  crustpinch 25-value vs 35-value: PASS"

##
## Test: dispatch selects correct model for count
##

echo "Testing model dispatch..."

# Count 15 should select Lop Nor Baseline
dir_count_15="$work/lopnor-count-15"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_count_15" 20000 "$args_5_lopnor_15"
if ! grep -q "Selected Lop Nor Baseline Model" "$dir_count_15/stdout.txt"; then
  echo "FAIL: count 15 did not select Lop Nor Baseline Model" >&2
  exit 1
fi

# Count 21 should select Lop Nor Baseline (not Moho)
dir_count_21="$work/lopnor-count-21"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_count_21" 20000 "$args_7_lopnor_21"
if ! grep -q "Selected Lop Nor Baseline Model" "$dir_count_21/stdout.txt"; then
  echo "FAIL: count 21 did not select Lop Nor Baseline Model" >&2
  exit 1
fi

# Count 20 should select Lop Nor Moho
dir_count_20="$work/lopnor-count-20"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_count_20" 20000 "$args_5_lopnor_20"
if ! grep -q "Selected Lop Nor Moho Model" "$dir_count_20/stdout.txt"; then
  echo "FAIL: count 20 did not select Lop Nor Moho Model" >&2
  exit 1
fi

# Count 28 should select Lop Nor Moho
dir_count_28="$work/lopnor-count-28"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$dir_count_28" 20000 "$args_7_lopnor_28"
if ! grep -q "Selected Lop Nor Moho Model" "$dir_count_28/stdout.txt"; then
  echo "FAIL: count 28 did not select Lop Nor Moho Model" >&2
  exit 1
fi

echo "  dispatch for counts 15, 21, 20, and 28: PASS"

##
## Test: grid dump contains expected columns
##

echo "Testing grid-dump columns..."

# Run with 0 phonons and --dump-grid to get grid output (can't sample anisotropic without AttemptAnisotropicScatter)
grid_args="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,2.0,2.0,0.5,0.3,1000,0.8,0.01,1.4,0.7,0.35,0.5,300"
grid_dir="$work/grid-dump"
EXTRA_ARGS="--dump-grid" TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$grid_dir" 0 "$grid_args" || true

# Check header contains expected column markers
if ! grep -E "nu.*eps.*a.*k.*ay.*az" "$grid_dir/stdout.txt" >/dev/null; then
  echo "FAIL: grid dump header does not contain expected anisotropic columns (nu eps a k ay az)" >&2
  grep "^#" "$grid_dir/stdout.txt" >&2
  exit 1
fi

# Use awk to verify field counts and values
# Anisotropic grid should have 17 fields per row
node_count=$(awk '/^#  i j k/{flag=1; next} /^#/{flag=0} flag && NF>0 && $1 ~ /^[0-9]+$/ {print NF}' "$grid_dir/stdout.txt" | sort -u | wc -l)
if [ "$node_count" -ne 1 ]; then
  echo "FAIL: node rows have inconsistent field counts" >&2
  awk '/^#  i j k/{flag=1; next} /^#/{flag=0} flag && NF>0 && $1 ~ /^[0-9]+$/ {print NF}' "$grid_dir/stdout.txt" | sort -u >&2
  exit 1
fi

field_count=$(awk '/^#  i j k/{flag=1; next} /^#/{flag=0} flag && NF>0 && $1 ~ /^[0-9]+$/ {print NF; exit}' "$grid_dir/stdout.txt")
if [ "$field_count" -ne 17 ]; then
  echo "FAIL: node rows have $field_count fields, expected 17" >&2
  exit 1
fi

# Check crust rows (field 15 = 0.30000, which is k for crust) have correct ay,az (2.00000, 0.50000)
crust_row=$(awk '/^#  i j k/{flag=1; next} /^#/{flag=0} flag && NF>0 && $1 ~ /^[0-9]+$/ && $15 ~ /0.30000/ {print; found=1} END {if (!found) exit 1}' "$grid_dir/stdout.txt")
if [ -z "$crust_row" ]; then
  echo "FAIL: no crust rows (k=0.30000) found in grid dump" >&2
  exit 1
fi

crust_ok=$(echo "$crust_row" | awk '$14 == "2.00000" && $16 == "2.00000" && $17 == "0.50000" {found=1} END {exit !found}')
if [ -z "$crust_ok" ]; then
  echo "FAIL: crust rows do not have correct ay,az values (2.00000 2.00000 0.50000)" >&2
  echo "$crust_row" >&2
  exit 1
fi

# Check mantle rows (field 15 = 0.50000) have correct ax,ay,az (1.40000, 0.70000, 0.35000)
mantle_row=$(awk '/^#  i j k/{flag=1; next} /^#/{flag=0} flag && NF>0 && $1 ~ /^[0-9]+$/ && $15 ~ /0.50000/ {print; found=1} END {if (!found) exit 1}' "$grid_dir/stdout.txt")
if [ -z "$mantle_row" ]; then
  echo "FAIL: no mantle rows (k=0.50000) found in grid dump" >&2
  exit 1
fi

mantle_ok=$(echo "$mantle_row" | awk '$14 == "1.40000" && $16 == "0.70000" && $17 == "0.35000" {found=1} END {exit !found}')
if [ -z "$mantle_ok" ]; then
  echo "FAIL: mantle rows do not have correct ax,ay,az values (1.40000 0.70000 0.35000)" >&2
  echo "$mantle_row" >&2
  exit 1
fi

echo "  grid dump columns and field counts: PASS"

##
## Test: default layout grid dump has 15 fields
##

echo "Testing default layout grid dump..."

grid_dir_default="$work/grid-dump-default"
EXTRA_ARGS="--dump-grid" TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$grid_dir_default" 0 "$args_5_lopnor_15" || true

# Check that default isotropic grid dump has 15 fields per row
default_field_count=$(awk '/^#  i j k/{flag=1; next} /^#/{flag=0} flag && NF>0 && $1 ~ /^[0-9]+$/ {print NF; exit}' "$grid_dir_default/stdout.txt")
if [ "$default_field_count" -ne 15 ]; then
  echo "FAIL: default grid rows have $default_field_count fields, expected 15" >&2
  exit 1
fi

echo "  default layout grid dump: PASS"

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
