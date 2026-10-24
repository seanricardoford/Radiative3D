#!/bin/sh
set -eu

# tests/compare_aniso_reference.sh
#
# Manual comparison harness (not part of the automated test_*.sh suite,
# and not auto-discovered by any runner): runs the anisotropic lopnor
# and crustpinch scenarios from tests/test_aniso_scenarios.sh against a
# prior ("reference") build of 'main', across thread counts and
# take-off-angle degrees, and requires every (current, reference) run
# pair to match under aniso_lib.sh's compare_runs in "strict" mode.
#
# Written for Task 2 of docs/superpowers/specs/2026-10-17-axisymmetric-
# sampler-performance-design.md, to confirm the AxiProposal bound-kernel
# restructure (PrepareBoundRep()/BoundRatio()/the geometry cache) didn't
# change any simulation output; reused as-is by later tasks (3-5) in
# that same plan for their own restructures.
#
# REQUIRED environment variable:
#   REF_EXE       Path to the reference 'main' executable to compare
#                 against (e.g. a build of the commit a restructure
#                 started from).
#
# Optional environment variable:
#   CMP_THREADS   The non-1 thread count to test.          Default: 12
#
# Builds the current tree's 'main' (into $ROOT/main) and, for every
# (scenario, THREADS, TOA_DEGREE) combination in
#   scenario    in {lopnor, lopnor-prolate, crustpinch}
#   THREADS     in {1, CMP_THREADS}
#   TOA_DEGREE  in {5, 7}            (7 exercises representatives with
#                                     nonzero radii -- see aniso_axi.cpp)
# runs the scenario once with R3D_EXE=$ROOT/main and once with
# R3D_EXE=$REF_EXE (both via tests/run_aniso_case.sh, at that script's
# default phonon counts), then compares the two run directories with
# compare_runs in "strict" mode. Prints one "PASS" or "DIFF" line per
# case (all 12 cases always run, so a human gets full visibility) and
# exits nonzero if any case differed.

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

# shellcheck source=tests/aniso_lib.sh
. "$ROOT/tests/aniso_lib.sh"

: "${REF_EXE:?REF_EXE must be set to a reference 'main' executable to compare against}"
if [ ! -x "$REF_EXE" ]; then
  echo "compare_aniso_reference.sh: REF_EXE '$REF_EXE' is not an executable file" >&2
  exit 1
fi

CMP_THREADS=${CMP_THREADS:-12}

make >/dev/null

work=$(mktemp -d "${TMPDIR:-/tmp}/radiative3d-aniso-reference.XXXXXX")
cleanup() {
  rm -rf "$work"
}
trap cleanup EXIT HUP INT TERM

# Anisotropic layouts, copied verbatim from tests/test_aniso_scenarios.sh
# (not sourced from it, since that script runs its own checks
# immediately on load) -- keep these two lines in sync with that file's
# lopnor_21/crustpinch_42 if either changes.
lopnor_21="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,2.0,2.0,0.5,0.3,1000,0.8,0.01,1.4,0.7,0.35,0.5,300"
crustpinch_42="0.8,0.01,0.2,0.2,0.2,0.2,200,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.6,0.3,0.15,0.4,1500,0.8,0.01,0.2,0.2,0.2,0.5,900,2.0,30.0,5.0,.3666667,.4736842,1,1"

# Prolate case (task-5b-brief.md, Step 4): lopnor_21's own layout and
# sediments region, with the crust and mantle regions replaced by the
# brief's ax==ay < az pair -- exercises BoundRatio()/RatioRun()'s
# "ah < az" branch (lopnor_21's own crust/mantle regions are both
# "ah >= az", crust trivially so since there ax==ay==az).
lopnor_prolate_21="0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,0.5,0.5,1.0,0.3,1000,0.8,0.01,0.7,0.7,1.4,0.5,300"

any_diff=0

##
## run_case LABEL SCENARIO ANISO_ARGS THREADS DEGREE
##
## Runs SCENARIO (the run_aniso_case.sh scenario name: "lopnor" or
## "crustpinch", used to pick the right source/receiver/grid command
## line) once against the current build and once against REF_EXE, both
## with ANISO_ARGS as --model-args and the same THREADS/DEGREE/seed/
## phonon-count, then compares the two run directories strictly. LABEL
## names the case in PASS/DIFF output and the two runs' work directories
## -- distinct from SCENARIO so two cases sharing a SCENARIO (e.g.
## "lopnor" and the prolate variant, both still the "lopnor"
## run_aniso_case.sh scenario) get distinguishable output and don't
## overwrite each other's work directory. Records any DIFF in
## 'any_diff'.
##
run_case() {
  _rc_label_name=$1
  _rc_scenario=$2
  _rc_args=$3
  _rc_threads=$4
  _rc_degree=$5
  _rc_label="$_rc_label_name threads=$_rc_threads degree=$_rc_degree"

  _rc_cur="$work/$_rc_label_name-t$_rc_threads-d$_rc_degree-cur"
  _rc_ref="$work/$_rc_label_name-t$_rc_threads-d$_rc_degree-ref"

  R3D_EXE="$ROOT/main" TOA_DEGREE="$_rc_degree" \
    "$ROOT/tests/run_aniso_case.sh" "$_rc_scenario" "$_rc_threads" "$_rc_cur" "" "$_rc_args"
  R3D_EXE="$REF_EXE" TOA_DEGREE="$_rc_degree" \
    "$ROOT/tests/run_aniso_case.sh" "$_rc_scenario" "$_rc_threads" "$_rc_ref" "" "$_rc_args"

  if compare_runs "$_rc_cur" "$_rc_ref" strict; then
    echo "PASS $_rc_label"
  else
    echo "DIFF $_rc_label" >&2
    any_diff=1
  fi
}

for _degree in 5 7; do
  for _threads in 1 "$CMP_THREADS"; do
    run_case lopnor          lopnor     "$lopnor_21"         "$_threads" "$_degree"
    run_case lopnor-prolate  lopnor     "$lopnor_prolate_21" "$_threads" "$_degree"
    run_case crustpinch      crustpinch "$crustpinch_42"     "$_threads" "$_degree"
  done
done

if [ "$any_diff" -ne 0 ]; then
  echo "compare_aniso_reference.sh: at least one case differed from REF_EXE" >&2
  exit 1
fi

echo "compare_aniso_reference.sh: all cases match REF_EXE"
