#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

BENCH_DEGREE=${BENCH_DEGREE:-9}
BENCH_LOPNOR_PHONONS=${BENCH_LOPNOR_PHONONS:-1000000}
BENCH_CRUSTPINCH_PHONONS=${BENCH_CRUSTPINCH_PHONONS:-200000}
BENCH_THREADS=${BENCH_THREADS:-auto}

if [ "$BENCH_THREADS" = "auto" ]; then
  BENCH_THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)
fi

if [ -z "${BASELINE_EXE:-}" ] || [ ! -x "$BASELINE_EXE" ]; then
  echo "BASELINE_EXE not set or executable" >&2
  exit 1
fi

make >/dev/null

work=$(mktemp -d)
trap "rm -rf $work" EXIT

iso_lopnor="0.8,0.01,0.5,0.2,50,0.8,0.01,0.5,0.3,1000,0.8,0.01,0.7,0.5,300"
iso_crustpinch="0.8,0.01,0.2,0.2,200,0.8,0.01,0.2,0.3,1500,0.8,0.01,0.2,0.3,1500,0.8,0.01,0.2,0.4,1500,0.8,0.01,0.2,0.5,900,2.0,30.0,5.0,.3666667,.4736842,1,1"

echo "# Anisotropic Scale-Length Performance Benchmark"
echo ""
echo "Degree: $BENCH_DEGREE, Threads: $BENCH_THREADS"
echo ""
echo "| Scenario | ε | Mode | r | Seconds | Slowdown vs isotropic | Gate |"
echo "|---|---|---|---|---:|---:|---|"

time_run() {
  exe=$1
  scenario=$2
  phonons=$3
  model=$4
  
  outdir="$work/run-$RANDOM"
  mkdir -p "$outdir"
  
  grid_compiled=1
  flatten=""
  range=1200
  if [ "$scenario" = "crustpinch" ]; then
    grid_compiled=5
    range=1500
  else
    flatten="--flatten"
  fi
  
  start=$(date +%s)
  
  "$exe" \
    --reports=ALL_OFF \
    --output-dir="$outdir" \
    --report-file=reports.dat \
    --mparams-outfile=out_mparams.octv \
    --num-phonons="$phonons" \
    --toa-degree="$BENCH_DEGREE" \
    --grid-compiled="$grid_compiled" \
    $flatten \
    --range="$range" \
    --frequency=2.0 \
    --timetolive=600 \
    --binsize=2.0 \
    --model-args="$model" \
    --seed=596 \
    --threads="$BENCH_THREADS" \
    >/dev/null 2>&1
  
  end=$(date +%s)
  echo $((end - start))
}

benchmark_scenario() {
  scenario=$1
  phonons=$2
  iso_model=$3

  echo "Measuring $scenario..." >&2

  times_base=""
  times_curr=""
  for i in 1 2 3; do
    echo "  isotropic run $i..." >&2
    t=$(time_run "$BASELINE_EXE" "$scenario" "$phonons" "$iso_model")
    times_base="$times_base $t"
    t=$(time_run "$ROOT/main" "$scenario" "$phonons" "$iso_model")
    times_curr="$times_curr $t"
  done

  med_base=$(echo $times_base | tr ' ' '\n' | sort -n | sed -n '2p')
  med_curr=$(echo $times_curr | tr ' ' '\n' | sort -n | sed -n '2p')
  iso_ratio=$(awk "BEGIN {printf \"%.3f\", $med_curr / $med_base}")

  echo "  isotropic ratio (current/baseline): $iso_ratio" >&2

  if awk "BEGIN {exit !($iso_ratio > 1.05)}"; then
    echo "ERROR: Isotropic ratio exceeds 1.05" >&2
    exit 1
  fi

  echo "  isotropic baseline timing..." >&2
  t_iso=$(time_run "$ROOT/main" "$scenario" "$phonons" "$iso_model")

  for eps_val in 0.01 0.05; do
    for r in 2 4 10; do
      echo "  $scenario r=$r ε=$eps_val..." >&2
      
      if [ "$scenario" = "lopnor" ]; then
        r_sqrt=$(awk "BEGIN {printf \"%.3f\", sqrt($r)}")
        axi_crust_x=$(awk "BEGIN {printf \"%.1f\", $r * 0.5}")
        axi_crust_z="0.5"
        ani_mantle_x=$(awk "BEGIN {printf \"%.2f\", $r * 0.7}")
        ani_mantle_y=$(awk "BEGIN {printf \"%.2f\", $r_sqrt * 0.7}")
        ani_mantle_z="0.7"
        
        model="0.8,$eps_val,0.5,0.2,50,0.8,$eps_val,$axi_crust_x,$axi_crust_x,$axi_crust_z,0.3,1000,0.8,$eps_val,$ani_mantle_x,$ani_mantle_y,$ani_mantle_z,0.5,300"
      else
        r_sqrt=$(awk "BEGIN {printf \"%.3f\", sqrt($r)}")
        crust_x=$(awk "BEGIN {printf \"%.2f\", $r * 0.2}")
        crust_z="0.2"
        mantle_x=$(awk "BEGIN {printf \"%.2f\", $r * 0.2}")
        mantle_y=$(awk "BEGIN {printf \"%.2f\", $r_sqrt * 0.2}")
        mantle_z="0.2"
        
        model="0.8,$eps_val,0.2,0.2,200,0.8,$eps_val,$crust_x,$crust_x,$crust_z,0.3,1500,0.8,$eps_val,$crust_x,$crust_x,$crust_z,0.3,1500,0.8,$eps_val,0.2,0.2,0.2,0.4,1500,0.8,$eps_val,$mantle_x,$mantle_y,$mantle_z,0.5,900,2.0,30.0,5.0,.3666667,.4736842,1,1"
      fi
      
      t=$(time_run "$ROOT/main" "$scenario" "$phonons" "$model")
      slowdown=$(awk "BEGIN {printf \"%.2f\", $t / $t_iso}")
      
      if awk "BEGIN {exit !($slowdown > 2.0)}"; then
        gate="GATE"
      else
        gate="ok"
      fi
      
      printf "| %s | %.2f | axi | %d | %d | %.2fx | %s |\n" "$scenario" "$eps_val" "$r" "$t" "$slowdown" "$gate"
    done
  done

  echo "" >&2
  echo "Isotropic ratio (current/baseline): $iso_ratio" >&2
}

benchmark_scenario lopnor "$BENCH_LOPNOR_PHONONS" "$iso_lopnor"
benchmark_scenario crustpinch "$BENCH_CRUSTPINCH_PHONONS" "$iso_crustpinch"

echo ""
echo "Benchmark complete" >&2
