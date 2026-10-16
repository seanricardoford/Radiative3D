#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
make >/dev/null

run_once() {
  output=$1
  ./main --num-phonons 2 --threads 1 --seed 123 --timetolive 0.001 \
    --report-file "$output.report" >"$output.log" 2>&1
  grep "__SIMULATION_COMPLETE__" "$output.log" >/dev/null
}

run_once tests/.seed-a
run_once tests/.seed-b
cmp tests/.seed-a.report tests/.seed-b.report

rm -f tests/.seed-a.* tests/.seed-b.*
echo "seed reproducibility checks passed"
