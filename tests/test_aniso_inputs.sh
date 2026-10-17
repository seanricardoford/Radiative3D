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

# Crust pinch: 25-value vs 35-value
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

# Crust pinch: 28-value vs 38-value (adds SediThick, CrustThick, MohoThick)
args_5_cp_28="${args_5_cp_25},2.0,30.0,10.0"
args_7_cp_38="${args_7_cp_35},2.0,30.0,10.0"

dir_5="$work/cp-28"
dir_7="$work/cp-38"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$dir_5" 5000 "$args_5_cp_28"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$dir_7" 5000 "$args_7_cp_38"

if ! compare_runs "$dir_5" "$dir_7" noargs; then
  echo "FAIL: 5-value and 7-value (isotropic) crustpinch 28→38 output differ" >&2
  exit 1
fi
echo "  crustpinch 28-value vs 38-value: PASS"

# Crust pinch: 32-value vs 42-value (adds PinchFrac, DepthFrac, SediFrac, MohoFrac)
args_5_cp_32="${args_5_cp_28},0.70,0.50,1.0,1.0"
args_7_cp_42="${args_7_cp_38},0.70,0.50,1.0,1.0"

dir_5="$work/cp-32"
dir_7="$work/cp-42"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$dir_5" 5000 "$args_5_cp_32"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$dir_7" 5000 "$args_7_cp_42"

if ! compare_runs "$dir_5" "$dir_7" noargs; then
  echo "FAIL: 5-value and 7-value (isotropic) crustpinch 32→42 output differ" >&2
  exit 1
fi
echo "  crustpinch 32-value vs 42-value: PASS"

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
# Anisotropic grid should have 17 fields per row (excluding *** rows which have 18 due to leading space)
# Data starts after the header line and ends at "#  END"
node_count=$(awk '/^#  i j k/{flag=1; next} /^#  END/{flag=0} flag && NF==17 {print NF}' "$grid_dir/stdout.txt" | sort -u | wc -l)
if [ "$node_count" -ne 1 ]; then
  echo "FAIL: node rows have inconsistent field counts" >&2
  # Debug output
  echo "Field count distribution for rows with 15-18 fields:" >&2
  awk '/^#  i j k/{flag=1; next} /^#  END/{flag=0} flag && NF>=15 && NF<=18 {print NF}' "$grid_dir/stdout.txt" | sort | uniq -c >&2
  exit 1
fi

field_count=$(awk '/^#  i j k/{flag=1; next} /^#  END/{flag=0} flag && NF==17 {print NF; exit}' "$grid_dir/stdout.txt")
if [ "$field_count" -ne 17 ]; then
  echo "FAIL: node rows have $field_count fields, expected 17" >&2
  exit 1
fi

# Check crust rows (field 15 = 0.30000, which is k for crust) have correct ay,az (2.00000, 0.50000)
crust_row=$(awk '/^#  i j k/{flag=1; next} /^#  END/{flag=0} flag && NF==17 && $15 ~ /0.30000/ {print; exit}' "$grid_dir/stdout.txt")
if [ -z "$crust_row" ]; then
  echo "FAIL: no crust rows (k=0.30000) found in grid dump" >&2
  exit 1
fi

# Extract and verify the specific fields
crust_ay=$(echo "$crust_row" | awk '{print $16}')
crust_az=$(echo "$crust_row" | awk '{print $17}')
if [ "$crust_ay" != "2.00000" ] || [ "$crust_az" != "0.50000" ]; then
  echo "FAIL: crust row does not have correct ay,az values (expected 2.00000 0.50000, got $crust_ay $crust_az)" >&2
  exit 1
fi

# Check mantle rows (field 15 = 0.50000) have correct ax,ay,az (1.40000, 0.70000, 0.35000)
mantle_row=$(awk '/^#  i j k/{flag=1; next} /^#  END/{flag=0} flag && NF==17 && $15 ~ /0.50000/ {print; exit}' "$grid_dir/stdout.txt")
if [ -z "$mantle_row" ]; then
  echo "FAIL: no mantle rows (k=0.50000) found in grid dump" >&2
  exit 1
fi

# Extract and verify the specific fields
mantle_ax=$(echo "$mantle_row" | awk '{print $14}')
mantle_ay=$(echo "$mantle_row" | awk '{print $16}')
mantle_az=$(echo "$mantle_row" | awk '{print $17}')
if [ "$mantle_ax" != "1.40000" ] || [ "$mantle_ay" != "0.70000" ] || [ "$mantle_az" != "0.35000" ]; then
  echo "FAIL: mantle row does not have correct ax,ay,az values (expected 1.40000 0.70000 0.35000, got $mantle_ax $mantle_ay $mantle_az)" >&2
  exit 1
fi

echo "  grid dump columns and field counts: PASS"

##
## Test: default layout grid dump has 15 fields
##

echo "Testing default layout grid dump..."

grid_dir_default="$work/grid-dump-default"
EXTRA_ARGS="--dump-grid" TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$grid_dir_default" 0 "$args_5_lopnor_15" || true

# Check that default isotropic grid dump has 15 fields per row (skip rows with leading space = ***)
default_field_count=$(awk '/^#  i j k/{flag=1; next} /^#  END/{flag=0} flag && NF==15 {print NF; exit}' "$grid_dir_default/stdout.txt")
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

##
## Test: anisotropic scatterer dump (Task 6 reporting)
##

echo "Testing anisotropic scatterer dump..."

# Zero-phonon run (PrintAllScatteringStats() always runs on model build,
# regardless of phonon count): one isotropic, one axisymmetric (AXI,
# ax=ay=2 != az=0.5), and one triaxial (ANI, 1.4/0.7/0.35) region.
scatdump_args="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,2.0,2.0,0.5,0.3,1000,0.8,0.01,1.4,0.7,0.35,0.5,300"
scatdump_dir="$work/scatterer-dump"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$scatdump_dir" 0 "$scatdump_args"
dump="$scatdump_dir/stdout.txt"

# Isotropic rows have 11 fields; AXI rows have 22; ANI rows have 23.
# (Scatterer-dump data rows are the non-"#"-prefixed lines between the
# BEGIN/END markers.)
iso_count=$(awk '/^#  BEGIN SCATTERER DUMP/{f=1;next} /^#  END SCATTERERS/{f=0} f && $0 !~ /^#/ && NF==11 {c++} END{print c+0}' "$dump")
if [ "$iso_count" -lt 1 ]; then
  echo "FAIL: expected at least one isotropic (11-field) scatterer row, got $iso_count" >&2
  exit 1
fi

axi_row=$(awk '/^#  BEGIN SCATTERER DUMP/{f=1;next} /^#  END SCATTERERS/{f=0} f && $0 !~ /^#/ && NF==22 {print; exit}' "$dump")
if [ -z "$axi_row" ]; then
  echo "FAIL: no AXI (22-field) scatterer row found" >&2
  exit 1
fi
axi_mode=$(echo "$axi_row" | awk '{print $14}')
if [ "$axi_mode" != "AXI" ]; then
  echo "FAIL: AXI row field 14 is '$axi_mode', expected AXI" >&2
  exit 1
fi
axi_ax=$(echo "$axi_row" | awk '{print $3}')
if [ "$axi_ax" != "2" ]; then
  echo "FAIL: AXI row field 3 (ax) is '$axi_ax', expected 2 (crust)" >&2
  exit 1
fi

ani_row=$(awk '/^#  BEGIN SCATTERER DUMP/{f=1;next} /^#  END SCATTERERS/{f=0} f && $0 !~ /^#/ && NF==23 {print; exit}' "$dump")
if [ -z "$ani_row" ]; then
  echo "FAIL: no ANI (23-field) scatterer row found" >&2
  exit 1
fi
ani_mode=$(echo "$ani_row" | awk '{print $14}')
if [ "$ani_mode" != "ANI" ]; then
  echo "FAIL: ANI row field 14 is '$ani_mode', expected ANI" >&2
  exit 1
fi

# el * MFP(S) (field 5 * field 9) must be positive on every row -- this
# is the same quantity do-halfspace-nearsrc50.sh computes from the
# scatterer dump, so an anisotropic dump must keep it meaningful.
neg_count=$(awk '/^#  BEGIN SCATTERER DUMP/{f=1;next} /^#  END SCATTERERS/{f=0} f && $0 !~ /^#/ && NF>=11 {if ($5*$9<=0) c++} END{print c+0}' "$dump")
if [ "$neg_count" -ne 0 ]; then
  echo "FAIL: found $neg_count scatterer row(s) with non-positive el*MFP(S)" >&2
  exit 1
fi

echo "  anisotropic scatterer dump field counts and values: PASS"

# The default (all-isotropic) layout's scatterer-dump header must not
# gain the appended 'ay' column text.
default_scat_dir="$work/scatterer-dump-default"
TOA_DEGREE=4 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$default_scat_dir" 0 "$args_5_lopnor_15"
header=$(awk '/^#  .*nu.*eps.*a.*kappa/{print; exit}' "$default_scat_dir/stdout.txt")
case "$header" in
  *ay*)
    echo "FAIL: default (isotropic) scatterer-dump header unexpectedly contains 'ay': $header" >&2
    exit 1
    ;;
esac

echo "  default layout scatterer-dump header has no ay: PASS"

echo ""
echo "All test_aniso_inputs checks passed"
