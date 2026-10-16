#!/bin/sh
set -eu

# tests/test_aniso_isotropic_regression.sh
#
# Safety net for the anisotropic scale-length work: with no C++ changes yet
# made, this proves the test harness itself (run_aniso_case.sh,
# aniso_lib.sh) reproduces a baseline build of commit aaec540 byte-for-byte
# for both the lopnor and crustpinch scenarios, and that per-region SCT
# scatter reports land in the expected depth/range bands. Later tasks reuse
# this harness once real anisotropic behavior is added, at which point
# these same comparisons are expected to start showing differences.

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

# shellcheck source=tests/aniso_lib.sh
. "$ROOT/tests/aniso_lib.sh"

make >/dev/null

work=$(mktemp -d "${TMPDIR:-/tmp}/radiative3d-aniso.XXXXXX")
cleanup() {
  if [ -n "${BASELINE_WORKTREE:-}" ]; then
    git worktree remove --force "$BASELINE_WORKTREE" >/dev/null 2>&1 || true
  fi
  rm -rf "$work"
}
trap cleanup EXIT HUP INT TERM

if [ -z "${BASELINE_EXE:-}" ]; then
  BASELINE_WORKTREE="$work/baseline"
  git worktree add --detach "$BASELINE_WORKTREE" "${BASELINE_REV:-aaec540}" >/dev/null
  make -C "$BASELINE_WORKTREE" >/dev/null
  BASELINE_EXE="$BASELINE_WORKTREE/main"
fi

for scenario in lopnor crustpinch; do

  baseline_dir="$work/$scenario-baseline"
  current_dir="$work/$scenario-current"

  R3D_EXE="$BASELINE_EXE" "$ROOT/tests/run_aniso_case.sh" "$scenario" 8 "$baseline_dir"
  R3D_EXE="$ROOT/main"    "$ROOT/tests/run_aniso_case.sh" "$scenario" 8 "$current_dir"

  if ! trace_has_signal "$baseline_dir/seis_traces_asc.dat"; then
    echo "$scenario: baseline traces are all zero; no phonons reached a seismometer" >&2
    exit 1
  fi
  if ! trace_has_signal "$current_dir/seis_traces_asc.dat"; then
    echo "$scenario: current traces are all zero; no phonons reached a seismometer" >&2
    exit 1
  fi

  compare_runs "$baseline_dir" "$current_dir" strict

  sct_dir="$work/$scenario-sct"
  REPORTS=SCT "$ROOT/tests/run_aniso_case.sh" "$scenario" 1 "$sct_dir" 20000

  case "$scenario" in
    lopnor)
      require_sct_bands lopnor "$sct_dir/reports.dat" sediments crust mantle
      ;;
    crustpinch)
      require_sct_bands crustpinch "$sct_dir/reports.dat" \
        sediments crust pinched moho mantle
      ;;
  esac

done

echo "isotropic regression checks passed"
