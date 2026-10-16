#!/bin/sh
set -eu

# tests/run_aniso_case.sh SCENARIO THREADS OUTDIR [PHONONS] [MODEL_ARGS]
#
# Runs the simulation-only part of do-lopnor.sh ("lopnor") or
# do-crustpinch.sh ("crustpinch") -- i.e. it stops after the executable
# returns, the same way run_lopnor_thread_case.sh does, without going on
# into the do-scripts' Octave figure generation.
#
# PHONONS defaults to 200000 (lopnor) or 50000 (crustpinch) when omitted or
# empty. MODEL_ARGS defaults to the do-script values for the scenario when
# omitted or empty.
#
# Environment variables:
#   R3D_EXE      Executable to run.            Default: "$ROOT/main"
#   TOA_DEGREE   --toa-degree value.            Default: 5
#   SEED         --seed value.                  Default: 596
#   REPORTS      --reports value.               Default: ALL_OFF
#   EXTRA_ARGS   Extra args appended verbatim (word-split), e.g.
#                EXTRA_ARGS=--dump-grid.        Default: (none)
#
# Writes stdout.txt, stderr.txt, reports.dat, out_mparams.octv,
# seis_traces_asc.dat, and seis_*.octv into OUTDIR. The command itself is
# run from inside OUTDIR, as run_lopnor_thread_case.sh does.

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

SCENARIO=$1
THREADS=$2
OUTDIR=$3
PHONONS=${4:-}
MODEL_ARGS=${5:-}

R3D_EXE=${R3D_EXE:-"$ROOT/main"}
TOA_DEGREE=${TOA_DEGREE:-5}
SEED=${SEED:-596}
REPORTS=${REPORTS:-ALL_OFF}
EXTRA_ARGS=${EXTRA_ARGS:-}

mkdir -p "$OUTDIR"
OUTDIR=$(CDPATH= cd -- "$OUTDIR" && pwd)

case "$SCENARIO" in
  lopnor)
    [ -n "$PHONONS" ] || PHONONS=200000
    [ -n "$MODEL_ARGS" ] || \
      MODEL_ARGS=0.8,0.01,0.5,0.2,50,0.8,0.01,0.5,0.3,1000,0.8,0.01,0.7,0.5,300

    XINXYZ=425.54,-169.53,0.98
    WUSXYZ=-390.04,-167.18,1.457
    MAKXYZ=-102.27,430.84,0.60
    SOURCETYP=SDR,125,40,90,0.0
    SOURCELOC=425.54,-169.53,-31.02

    set -- "$R3D_EXE" \
      --reports="$REPORTS" \
      --output-dir="$OUTDIR" \
      --report-file=reports.dat \
      --mparams-outfile=out_mparams.octv \
      --num-phonons="$PHONONS" \
      --toa-degree="$TOA_DEGREE" \
      --source="$SOURCETYP" \
      --source-loc="$SOURCELOC" \
      --frequency=2.0 \
      --timetolive=600 \
      --binsize=2.0 \
      --grid-compiled=1 \
      --range=1200 \
      --flatten \
      --model-args="$MODEL_ARGS" \
      --seed="$SEED" \
      --threads="$THREADS" \
      --seis-p2p="$XINXYZ,$WUSXYZ,1.0,2.0,40.0,8" \
      --seis-p2p="$XINXYZ,$MAKXYZ,1.0,2.0,40.0,8" \
      $EXTRA_ARGS
    ;;

  crustpinch)
    [ -n "$PHONONS" ] || PHONONS=50000
    [ -n "$MODEL_ARGS" ] || \
      MODEL_ARGS=0.8,0.01,0.20,0.2,200,0.8,0.01,0.20,0.3,1500,0.8,0.01,0.20,0.3,1500,0.8,0.01,0.20,0.4,1500,0.8,0.01,0.20,0.5,900,2.0,30.0,5.0,.3666667,.4736842,1,1

    SOURCETYP=SDR,22.5,90,0
    SOURCELOC=0,0,-10

    set -- "$R3D_EXE" \
      --reports="$REPORTS" \
      --output-dir="$OUTDIR" \
      --report-file=reports.dat \
      --mparams-outfile=out_mparams.octv \
      --num-phonons="$PHONONS" \
      --toa-degree="$TOA_DEGREE" \
      --source="$SOURCETYP" \
      --source-loc="$SOURCELOC" \
      --frequency=2.0 \
      --timetolive=600 \
      --binsize=2.0 \
      --grid-compiled=5 \
      --model-args="$MODEL_ARGS" \
      --seed="$SEED" \
      --threads="$THREADS" \
      --seis-p2p=0,67.5,0,950,67.5,0,1.0,2.0,40.0,8 \
      --seis-p2p=0,112.5,0,950,112.5,0,1.0,2.0,40.0,8 \
      --seis-p2p=0,90,0,950,90,0,1.0,2.0,40.0,8 \
      $EXTRA_ARGS
    ;;

  *)
    echo "run_aniso_case.sh: unknown scenario '$SCENARIO' (want lopnor|crustpinch)" >&2
    exit 1
    ;;
esac

(
  cd "$OUTDIR"
  "$@" >stdout.txt 2>stderr.txt
)
