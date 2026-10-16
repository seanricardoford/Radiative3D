#!/bin/sh
set -eu

# tests/test_anisotropy_unit.sh [FILTER]
#
# Runs the C++ unit-check harness for the anisotropic heterogeneity
# scale-length work (tests/anisotropy_tests.cpp), built via `make
# anisotropy-tests` into build/anisotropy_tests.
#
# With no argument: runs every test whose name doesn't start with
# "isolated_" in a single process, then runs each "isolated_*" test (if
# any) in its own fresh process -- those tests exercise overrides that
# are class-static and can't be reset once set, so each needs a clean
# process. No isolated_* tests exist yet, so that loop currently finds
# none and does nothing.
#
# With a FILTER argument: runs only `build/anisotropy_tests FILTER`.
# (anisotropy_tests itself only runs an isolated_* test when FILTER
# names it exactly; see tests/anisotropy_tests.cpp.)

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

make >/dev/null
make anisotropy-tests >/dev/null

FILTER="${1:-}"

if [ -n "$FILTER" ]; then
  build/anisotropy_tests "$FILTER"
else
  build/anisotropy_tests

  names=$(build/anisotropy_tests --list)
  for name in $names; do
    case "$name" in
      isolated_*)
        build/anisotropy_tests "$name"
        ;;
    esac
  done
fi

echo "anisotropy unit checks passed"
