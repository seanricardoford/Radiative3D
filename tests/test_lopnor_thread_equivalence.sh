#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

make >/dev/null

work=$(mktemp -d "${TMPDIR:-/tmp}/radiative3d-lopnor.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM

run_case() {
  "$ROOT/tests/run_lopnor_thread_case.sh" "$1" "$work/$2"
}

# Identical all-zero traces would make the comparisons below pass trivially,
# so require at least one nonzero energy or count value inside a trace block.
trace_has_signal() {
  awk '
    /#### BEGIN TRACE ####/ { in_trace = 1; next }
    /#### END TRACE ####/ { in_trace = 0; next }
    in_trace { for (i = 1; i <= NF; i++) if ($i + 0 != 0) found = 1 }
    END { exit !found }
  ' "$1"
}

run_case 1 one
run_case 12 twelve
run_case auto automatic

for case in one twelve automatic; do
  grep '__SIMULATION_COMPLETE__' "$work/$case/stdout.txt" >/dev/null
  grep '^Workers:' "$work/$case/stdout.txt" >/dev/null
  test -s "$work/$case/seis_traces_asc.dat"
  grep -E '^[[:space:]]*[0-9.eE+-]+[[:space:]]+' "$work/$case/seis_traces_asc.dat" >/dev/null
  if ! trace_has_signal "$work/$case/seis_traces_asc.dat"; then
    echo "Lop Nor '$case' traces are all zero; no phonons reached a seismometer" >&2
    exit 1
  fi
done

normalize_trace() {
  sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$1"
}

normalize_trace "$work/one/seis_traces_asc.dat" > "$work/one.trace"
normalize_trace "$work/twelve/seis_traces_asc.dat" > "$work/twelve.trace"
normalize_trace "$work/automatic/seis_traces_asc.dat" > "$work/automatic.trace"
cmp "$work/one.trace" "$work/twelve.trace"
cmp "$work/one.trace" "$work/automatic.trace"

echo "Lop Nor thread-equivalence checks passed without Octave"
