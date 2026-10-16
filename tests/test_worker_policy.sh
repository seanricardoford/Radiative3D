#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
make >/dev/null

run_and_check() {
  expected=$1
  shift
  ./main --num-phonons 1 --threads "$@" --timetolive 0.001 2>&1 | grep "Workers: $expected" >/dev/null
}

run_and_check 1 1
run_and_check 7 7

echo "worker policy checks passed"
