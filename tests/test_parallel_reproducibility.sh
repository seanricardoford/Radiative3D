#!/bin/bash

set -eu

test_root=$(mktemp -d /tmp/r3d-parallel-repro.XXXXXX)
trap 'rm -rf "$test_root"' EXIT

model_args="0.8,0.01,0.5,0.2,50,0.8,0.01,0.5,0.3,1000,0.8,0.01,0.7,0.5,300"

run_case() {
  local label=$1
  local workers=$2
  shift 2
  local output_dir="$test_root/$label-workers-$workers"
  mkdir -p "$output_dir"

  ./main \
    --reports=ALL_OFF \
    --output-dir="$output_dir" \
    --report-file=reports.dat \
    --num-phonons=512 \
    --toa-degree=3 \
    --source=EXPL \
    --source-loc=425.54,-169.53,-1.02 \
    --frequency=2.0 \
    --timetolive=10 \
    --binsize=1.0 \
    --grid-compiled=1 \
    --range=1200 \
    --flatten \
    --model-args="$model_args" \
    --seis-p2p=425.54,-169.53,0.98,-390.04,-167.18,1.457,1.0,2.0,2.0,4 \
    --workers="$workers" \
    --seed=0x5eedc0de12345678 \
    "$@" \
    > "$output_dir/stdout.txt"

  rg --fixed-strings --quiet '@@ __SIMULATION_COMPLETE__' \
    "$output_dir/stdout.txt"
  rg --regexp 'Invalidity:[[:space:]]+0([[:space:]]|$)' \
    "$output_dir/stdout.txt"

  rg 'Loss surfaces:|Timeout:|Invalidity:' "$output_dir/stdout.txt" \
    > "$output_dir/summary.txt"
}

run_case isotropic 1
run_case isotropic 2
diff -u "$test_root/isotropic-workers-1/summary.txt" \
        "$test_root/isotropic-workers-2/summary.txt"

run_case anisotropic 1 --scatter-horizontal=1.25 --scatter-vertical=0.625
run_case anisotropic 2 --scatter-horizontal=1.25 --scatter-vertical=0.625
diff -u "$test_root/anisotropic-workers-1/summary.txt" \
        "$test_root/anisotropic-workers-2/summary.txt"

echo "parallel reproducibility tests passed"
