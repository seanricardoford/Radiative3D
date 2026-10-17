#!/bin/sh
set -eu

# tests/test_aniso_scenarios.sh
#
# End-to-end integration tests for anisotropic scale-length scattering.
# Validates that:
#  - Anisotropic runs complete with __SIMULATION_COMPLETE__
#  - Two THREADS=1 runs with anisotropic args produce identical output (strict)
#  - THREADS=12 traces match THREADS=1 traces after address masking
#  - Anisotropic traces contain nonzero signal (trace_has_signal)
#  - Anisotropic traces differ from isotropic default-args traces (! cmp -s)
#  - REPORTS=SCT with 20000 phonons hits the expected anisotropic-region bands

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

# shellcheck source=tests/aniso_lib.sh
. "$ROOT/tests/aniso_lib.sh"

make >/dev/null

work=$(mktemp -d "${TMPDIR:-/tmp}/radiative3d-aniso-scenarios.XXXXXX}")
cleanup() {
  rm -rf "$work"
}
trap cleanup EXIT HUP INT TERM

# Anisotropic layouts from the brief
lopnor_21="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,2.0,2.0,0.5,0.3,1000,0.8,0.01,1.4,0.7,0.35,0.5,300"
crustpinch_42="0.8,0.01,0.2,0.2,0.2,0.2,200,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.6,0.3,0.15,0.4,1500,0.8,0.01,0.2,0.2,0.2,0.5,900,2.0,30.0,5.0,.3666667,.4736842,1,1"

# Default isotropic layouts for comparison
lopnor_isotropic="0.8,0.01,0.5,0.2,50,0.8,0.01,0.5,0.3,1000,0.8,0.01,0.7,0.5,300"
crustpinch_isotropic="0.8,0.01,0.20,0.2,200,0.8,0.01,0.20,0.3,1500,0.8,0.01,0.20,0.3,1500,0.8,0.01,0.20,0.4,1500,0.8,0.01,0.20,0.5,900,2.0,30.0,5.0,.3666667,.4736842,1,1"

##
## lopnor scenario
##

echo "Testing lopnor anisotropic scenario..."

# Run 1: THREADS=1, anisotropic, first time
lopnor_a1="$work/lopnor-a1"
TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$lopnor_a1" "" "$lopnor_21"

if ! grep -q "__SIMULATION_COMPLETE__" "$lopnor_a1/stdout.txt"; then
  echo "FAIL: lopnor run 1 did not complete" >&2
  exit 1
fi

# Run 2: THREADS=1, anisotropic, second time (for strict compare)
lopnor_a2="$work/lopnor-a2"
TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$lopnor_a2" "" "$lopnor_21"

if ! compare_runs "$lopnor_a1" "$lopnor_a2" strict; then
  echo "FAIL: lopnor THREADS=1 anisotropic runs differ (not reproducible)" >&2
  exit 1
fi
echo "  lopnor THREADS=1 reproducible: PASS"

# Check signal
if ! trace_has_signal "$lopnor_a1/seis_traces_asc.dat"; then
  echo "FAIL: lopnor anisotropic traces are all zero; no phonons reached seismometer" >&2
  exit 1
fi
echo "  lopnor anisotropic trace has signal: PASS"

# Run 3: THREADS=12, anisotropic (for threading check)
lopnor_a12="$work/lopnor-a12"
TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" lopnor 12 "$lopnor_a12" "" "$lopnor_21"

# Compare THREADS=12 traces to THREADS=1 after address masking
sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$lopnor_a1/seis_traces_asc.dat" \
  > "$lopnor_a1/.seis_traces_asc.normalized"
sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$lopnor_a12/seis_traces_asc.dat" \
  > "$lopnor_a12/.seis_traces_asc.normalized"

if ! cmp "$lopnor_a1/.seis_traces_asc.normalized" "$lopnor_a12/.seis_traces_asc.normalized"; then
  echo "FAIL: lopnor THREADS=12 traces differ from THREADS=1 after address masking" >&2
  exit 1
fi
echo "  lopnor THREADS=12 matches THREADS=1 (address-masked): PASS"

# Run 4: default isotropic args (for diff comparison)
lopnor_iso="$work/lopnor-iso"
TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$lopnor_iso" "" "$lopnor_isotropic"

if cmp "$lopnor_a1/seis_traces_asc.dat" "$lopnor_iso/seis_traces_asc.dat"; then
  echo "FAIL: lopnor anisotropic and isotropic traces are identical (expected difference)" >&2
  exit 1
fi
echo "  lopnor anisotropic traces differ from isotropic: PASS"

# Run 5: REPORTS=SCT with 20000 phonons for band checks
lopnor_sct="$work/lopnor-sct"
REPORTS=SCT TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" lopnor 1 "$lopnor_sct" 20000 "$lopnor_21"

require_sct_bands lopnor "$lopnor_sct/reports.dat" crust mantle
echo "  lopnor SCT bands (crust, mantle) present: PASS"

##
## crustpinch scenario
##

echo "Testing crustpinch anisotropic scenario..."

# Run 1: THREADS=1, anisotropic, first time
cp_a1="$work/cp-a1"
TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$cp_a1" "" "$crustpinch_42"

if ! grep -q "__SIMULATION_COMPLETE__" "$cp_a1/stdout.txt"; then
  echo "FAIL: crustpinch run 1 did not complete" >&2
  exit 1
fi

# Run 2: THREADS=1, anisotropic, second time (for strict compare)
cp_a2="$work/cp-a2"
TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$cp_a2" "" "$crustpinch_42"

if ! compare_runs "$cp_a1" "$cp_a2" strict; then
  echo "FAIL: crustpinch THREADS=1 anisotropic runs differ (not reproducible)" >&2
  exit 1
fi
echo "  crustpinch THREADS=1 reproducible: PASS"

# Check signal
if ! trace_has_signal "$cp_a1/seis_traces_asc.dat"; then
  echo "FAIL: crustpinch anisotropic traces are all zero; no phonons reached seismometer" >&2
  exit 1
fi
echo "  crustpinch anisotropic trace has signal: PASS"

# Run 3: THREADS=12, anisotropic (for threading check)
cp_a12="$work/cp-a12"
TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" crustpinch 12 "$cp_a12" "" "$crustpinch_42"

# Compare THREADS=12 traces to THREADS=1 after address masking
sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$cp_a1/seis_traces_asc.dat" \
  > "$cp_a1/.seis_traces_asc.normalized"
sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$cp_a12/seis_traces_asc.dat" \
  > "$cp_a12/.seis_traces_asc.normalized"

if ! cmp "$cp_a1/.seis_traces_asc.normalized" "$cp_a12/.seis_traces_asc.normalized"; then
  echo "FAIL: crustpinch THREADS=12 traces differ from THREADS=1 after address masking" >&2
  exit 1
fi
echo "  crustpinch THREADS=12 matches THREADS=1 (address-masked): PASS"

# Run 4: default isotropic args (for diff comparison)
cp_iso="$work/cp-iso"
TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$cp_iso" "" "$crustpinch_isotropic"

if cmp "$cp_a1/seis_traces_asc.dat" "$cp_iso/seis_traces_asc.dat"; then
  echo "FAIL: crustpinch anisotropic and isotropic traces are identical (expected difference)" >&2
  exit 1
fi
echo "  crustpinch anisotropic traces differ from isotropic: PASS"

# Run 5: REPORTS=SCT with 20000 phonons for band checks
cp_sct="$work/cp-sct"
REPORTS=SCT TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" crustpinch 1 "$cp_sct" 20000 "$crustpinch_42"

require_sct_bands crustpinch "$cp_sct/reports.dat" crust pinched moho
echo "  crustpinch SCT bands (crust, pinched, moho) present: PASS"

echo "anisotropic scenario checks passed"
