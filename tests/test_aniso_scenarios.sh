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

##
## check_scenario SCENARIO ANISO_ARGS BAND...
##
## Runs comprehensive checks for one anisotropic scenario:
##  - Two THREADS=1 runs for strict reproducibility
##  - One THREADS=12 run for threading equivalence (address-masked)
##  - Signal validation (trace_has_signal)
##  - Difference from default-args isotropic run
##  - SCT band presence for REPORTS=SCT
##
check_scenario() {
  _scenario=$1
  _aniso_args=$2
  shift 2
  # "$@" now contains the band names

  echo "Testing $_scenario anisotropic scenario..."

  # Run 1: THREADS=1, anisotropic, first time
  _a1="$work/$_scenario-a1"
  TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" "$_scenario" 1 "$_a1" "" "$_aniso_args"

  if ! grep -q "__SIMULATION_COMPLETE__" "$_a1/stdout.txt"; then
    echo "FAIL: $_scenario run 1 did not complete" >&2
    exit 1
  fi

  # Run 2: THREADS=1, anisotropic, second time (for strict compare)
  _a2="$work/$_scenario-a2"
  TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" "$_scenario" 1 "$_a2" "" "$_aniso_args"

  if ! compare_runs "$_a1" "$_a2" strict; then
    echo "FAIL: $_scenario THREADS=1 anisotropic runs differ (not reproducible)" >&2
    exit 1
  fi
  echo "  $_scenario THREADS=1 reproducible: PASS"

  # Check signal
  if ! trace_has_signal "$_a1/seis_traces_asc.dat"; then
    echo "FAIL: $_scenario anisotropic traces are all zero; no phonons reached seismometer" >&2
    exit 1
  fi
  echo "  $_scenario anisotropic trace has signal: PASS"

  # Run 3: THREADS=12, anisotropic (for threading check)
  _a12="$work/$_scenario-a12"
  TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" "$_scenario" 12 "$_a12" "" "$_aniso_args"

  # Compare THREADS=12 traces to THREADS=1 after address masking
  sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$_a1/seis_traces_asc.dat" \
    > "$_a1/.seis_traces_asc.normalized"
  sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$_a12/seis_traces_asc.dat" \
    > "$_a12/.seis_traces_asc.normalized"

  if ! cmp "$_a1/.seis_traces_asc.normalized" "$_a12/.seis_traces_asc.normalized"; then
    echo "FAIL: $_scenario THREADS=12 traces differ from THREADS=1 after address masking" >&2
    exit 1
  fi
  echo "  $_scenario THREADS=12 matches THREADS=1 (address-masked): PASS"

  # Run 4: default isotropic args (for diff comparison)
  _iso="$work/$_scenario-iso"
  TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" "$_scenario" 1 "$_iso" "" ""

  if cmp -s "$_a1/seis_traces_asc.dat" "$_iso/seis_traces_asc.dat"; then
    echo "FAIL: $_scenario anisotropic and isotropic traces are identical (expected difference)" >&2
    exit 1
  fi
  echo "  $_scenario anisotropic traces differ from isotropic: PASS"

  # Run 5: REPORTS=SCT with 20000 phonons for band checks
  _sct="$work/$_scenario-sct"
  REPORTS=SCT TOA_DEGREE=5 "$ROOT/tests/run_aniso_case.sh" "$_scenario" 1 "$_sct" 20000 "$_aniso_args"

  require_sct_bands "$_scenario" "$_sct/reports.dat" "$@"
  echo "  $_scenario SCT bands ($*) present: PASS"
}

# Anisotropic layouts from the brief
lopnor_21="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,2.0,2.0,0.5,0.3,1000,0.8,0.01,1.4,0.7,0.35,0.5,300"
crustpinch_42="0.8,0.01,0.2,0.2,0.2,0.2,200,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.6,0.3,0.15,0.4,1500,0.8,0.01,0.2,0.2,0.2,0.5,900,2.0,30.0,5.0,.3666667,.4736842,1,1"

check_scenario lopnor "$lopnor_21" crust mantle
check_scenario crustpinch "$crustpinch_42" crust pinched moho

echo "anisotropic scenario checks passed"
