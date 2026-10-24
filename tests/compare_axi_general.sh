#!/bin/sh
set -eu
#
# tests/compare_axi_general.sh
#
# Manual scenario-level validation that the axisymmetric sampler agrees
# with the general sampler on the same (or a near-identical) medium --
# the theory-correctness companion to tests/benchmark_aniso.sh's
# same-medium performance gate. Not named test_*, so the full test-suite
# check skips it; run it directly. See docs/superpowers/specs/
# 2026-10-17-axisymmetric-sampler-performance-design.md, decision 6
# ("Validation of the theory with ax ~= ay"), and .superpowers/sdd/
# 2026-10-17-axisymmetric-sampler-performance/task-6-brief.md.
#
# Groups, per config (scenario, eps, r):
#   A   mode "axi": ax == ay == r*a (AUTO picks the axisymmetric sampler).
#   G   mode "gen": ax == r*a, ay == ax*(1+1e-5) (AUTO picks the general
#       sampler, on a medium indistinguishable from A's to 1e-5).
#   A'  mode "axi" again, with a different block of seeds: a null
#       control -- A and A' are statistically identical runs, so a
#       p-value test flagging a difference here would mean the test
#       itself (not the theory) is unreliable at this K/permutation
#       count.
#
# Seeds: replicate k of A uses 596 + k*2^32 (k = 1..K, K = CMP_REPLICATES);
# G uses k = K+1..2K; A' uses k = 2K+1..3K. 2^32 is far larger than any
# phonon count used here, so MixRandomSeed(seed ^ phonon_id)
# (model.cpp) never lets two replicates' per-phonon streams collide --
# see docs/superpowers/specs/2026-10-17-axisymmetric-sampler-performance-
# design.md, decision 6's seed rule. 596 + k*2^32 is always an exact
# integer in IEEE double precision for any k used here (doubles are
# exact up to 2^53), so it is computed with awk (seed_for(), below)
# rather than shell arithmetic, and passed through run_aniso_case.sh's
# SEED env var to --seed, which parses an unsigned 64-bit integer.
#
# Statistic: from each run's seis_traces_asc.dat, every station's
# energy (the sum of its X+Y+Z axis-energy columns, equal to the sum of
# its P+S by-type columns -- dataout.cpp's Seismometer::RecordHit()
# partitions the same per-phonon energy both ways) is summed into
# 10-second windows (5 of the simulation's 2-second bins each -- see
# run_aniso_case.sh's --binsize=2.0). Each (station, window) pair is one
# "window" below. For a given pair of groups being compared, a window is
# kept only if every one of the 2K runs being compared (all K of the
# first group and all K of the second) has positive energy there --
# i.e. "every replicate of both compared groups has positive energy", a
# per-comparison filter, so A-vs-G and A-vs-A' can (and in general will,
# by a handful of windows) keep slightly different window sets; the
# table's single "windows" column reports A-vs-G's count (A' is the
# null control, not the primary comparison). For each kept window, with
# K values per group, z = (mean_A - mean_G) / sqrt(var_A/K + var_G/K)
# (sample variance, K-1 denominator), and T = mean(z^2) over kept
# windows. A window whose two groups both have exactly zero variance
# (so the z denominator is itself zero) is excluded from T's average
# rather than producing a divide-by-zero; this has not been observed in
# practice (real, phonon-shot-noise energy sums) and would only happen
# for a pathologically degenerate run.
#
# p-value: the fraction of CMP_PERMUTATIONS random relabellings of the
# 2K runs (into two new groups of K) with T >= T_obs, using a fixed-seed
# (not time-seeded) linear congruential generator implemented in awk
# (lcg_next()/lcg_uniform(), in window_stat() below; the same seed,
# LCG_SEED, is reused for every comparison and every config, since nothing
# about this test should depend on wall-clock time). The observed
# labelling is NOT separately added to the N permutations tested (a
# label relabelling drawn at random could by chance reproduce it, but it
# is not forced in) -- with CMP_PERMUTATIONS=1000 (or the smoke test's
# 200) either convention changes p by at most 1/N, which does not matter
# at the 0.01 gate. Computed independently for A-vs-G and A-vs-A'.
#
# Exit status: nonzero if any simulation run fails, or if either
# p-value (A-vs-G or A-vs-A') is < 0.01 for any config.
#
# Environment, with defaults:
#   CMP_DEGREE               --toa-degree for all runs.      Default: 7
#   CMP_LOPNOR_PHONONS        --num-phonons for lopnor runs.  Default: 250000
#   CMP_CRUSTPINCH_PHONONS    --num-phonons for crustpinch.   Default: 50000
#   CMP_REPLICATES            K, replicates per group.        Default: 8
#   CMP_THREADS               --threads for all runs, or "auto" for the
#                             logical CPU count.              Default: auto
#   CMP_PERMUTATIONS          Random relabellings per comparison per
#                             config.                         Default: 1000
#   CMP_CONFIGS               Which (scenario, eps, r) configs to run;
#                             see syntax below.
#                             Default: all eight of {lopnor, crustpinch}
#                             x eps {do-script value 0.01, 0.05} x r {4, 10}
#
# CMP_CONFIGS syntax: a space-separated list of "SCENARIO:EPS:R" tokens,
# e.g.
#   CMP_CONFIGS="lopnor:0.05:4 crustpinch:0.01:10"
# SCENARIO is "lopnor" or "crustpinch" (tests/run_aniso_case.sh's
# scenarios); EPS and R are passed straight through to
# tests/aniso_lib.sh's scenario_model(). A single config, e.g.
#   CMP_CONFIGS="lopnor:0.05:4"
# is how the smoke test below exercises the script cheaply. A token that
# isn't exactly three non-empty, colon-separated fields (e.g. a missing
# "r", or an extra ":") is rejected with an error and a nonzero exit --
# see cmp_config_fields() below.
#
# Output: a Markdown table, one row per config, with columns scenario,
# eps, r, windows, p(A,G), p(A,A'), median wall-clock seconds for A,
# median wall-clock seconds for G, and the G/A time ratio (informational
# here; tests/benchmark_aniso.sh's median-of-3 same-medium gate is the
# authoritative performance check). Progress messages go to stderr.

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
. "$ROOT/tests/aniso_lib.sh"

CMP_DEGREE=${CMP_DEGREE:-7}
CMP_LOPNOR_PHONONS=${CMP_LOPNOR_PHONONS:-250000}
CMP_CRUSTPINCH_PHONONS=${CMP_CRUSTPINCH_PHONONS:-50000}
CMP_REPLICATES=${CMP_REPLICATES:-8}
CMP_THREADS=${CMP_THREADS:-auto}
CMP_PERMUTATIONS=${CMP_PERMUTATIONS:-1000}

DOSCRIPT_EPS=0.01   # The epsilon value used by do-lopnor.sh and
                    # do-crustpinch.sh for every region, by default --
                    # same constant as tests/benchmark_aniso.sh's.
CMP_CONFIGS_DEFAULT="lopnor:$DOSCRIPT_EPS:4 lopnor:$DOSCRIPT_EPS:10 lopnor:0.05:4 lopnor:0.05:10 crustpinch:$DOSCRIPT_EPS:4 crustpinch:$DOSCRIPT_EPS:10 crustpinch:0.05:4 crustpinch:0.05:10"
CMP_CONFIGS=${CMP_CONFIGS:-$CMP_CONFIGS_DEFAULT}

LCG_SEED=20261017   # Fixed (not time-based) seed for window_stat()'s
                    # permutation LCG; see the header comment above.

if [ "$CMP_REPLICATES" -lt 2 ]; then
  echo "ERROR: CMP_REPLICATES must be >= 2 (need a sample variance per group)" >&2
  exit 1
fi

if [ "$CMP_THREADS" = "auto" ]; then
  CMP_THREADS=$(auto_threads)
fi

echo "Building current executable..." >&2
make >/dev/null

R3D_EXE_CURRENT="$ROOT/main"
if [ ! -x "$R3D_EXE_CURRENT" ]; then
  echo "ERROR: build did not produce an executable at $R3D_EXE_CURRENT" >&2
  exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

GATE_FAIL=0
K=$CMP_REPLICATES

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

# seed_for K_OFFSET -> prints 596 + K_OFFSET*2^32 as an exact unsigned
# 64-bit integer (see the header comment's seed rule).
seed_for() {
  awk -v k="$1" 'BEGIN{printf "%.0f", 596 + k*4294967296}'
}

# cmp_config_fields TOKEN  -> prints "SCENARIO EPS R" (space-separated) on
# stdout. Exits the whole script nonzero, with an error to stderr, if
# TOKEN is not exactly three non-empty, colon-separated fields -- the old
# "${tok%%:*}"/"${tok#*:}" field-splitting alone would silently turn a
# malformed token (e.g. "lopnor", missing eps and r) into a scenario of
# "lopnor" and an eps/r of "lopnor" too, rather than failing loudly.
cmp_config_fields() {
  awk -v tok="$1" '
    BEGIN {
      n = split(tok, f, ":")
      if (n != 3 || f[1] == "" || f[2] == "" || f[3] == "") {
        printf "ERROR: CMP_CONFIGS token \"%s\" is not SCENARIO:EPS:R (three non-empty colon-separated fields)\n", tok > "/dev/stderr"
        exit 1
      }
      print f[1], f[2], f[3]
    }
  '
}

# phonons_for_scenario SCENARIO -> prints the configured phonon count.
phonons_for_scenario() {
  case "$1" in
    lopnor) echo "$CMP_LOPNOR_PHONONS" ;;
    crustpinch) echo "$CMP_CRUSTPINCH_PHONONS" ;;
    *) echo "phonons_for_scenario: unknown scenario '$1'" >&2; exit 1 ;;
  esac
}

# extract_energy SEIS_TRACES_FILE OUT_TSV
#
#   Writes "station<TAB>window<TAB>energy" lines to OUT_TSV: one per
#   (station, 10-second window) cell with a complete 5 bins of 2s each
#   (dataout.cpp's Seismometer::OutputTrace() always emits exactly
#   cmNumBins rows per station, so every cell is complete in practice;
#   the cnt[key]==5 check is defensive). energy is that cell's summed
#   X+Y+Z axis-energy columns -- see the header comment.
extract_energy() {
  awk '
    /Seismometer Number:/ {
      for (i=1;i<=NF;i++) if ($i=="Number:") { station=$(i+1)+0; break }
      next
    }
    /#### BEGIN TRACE ####/ { intrace=1; bin=0; next }
    /#### END TRACE ####/   { intrace=0; next }
    intrace {
      w = int(bin/5)
      key = station SUBSEP w
      sum[key] += ($1+$2+$3)
      cnt[key]++
      bin++
      next
    }
    END {
      for (key in sum) {
        if (cnt[key]==5) {
          split(key, p, SUBSEP)
          printf "%s\t%s\t%.10e\n", p[1], p[2], sum[key]
        }
      }
    }
  ' "$1" > "$2"
}

# median_of FILE -> prints the median of the numbers in FILE (one per line).
median_of() {
  sort -n "$1" | awk '{a[NR]=$1} END{
    n=NR
    if (n%2==1) printf "%.3f\n", a[(n+1)/2]
    else        printf "%.3f\n", (a[n/2]+a[n/2+1])/2
  }'
}

# run_group SCENARIO EPS R MODE KSTART KEND PREFIX
#
#   Runs replicates KSTART..KEND of MODE (axi|gen), via
#   tests/run_aniso_case.sh on the current build. Writes one
#   energy.tsv path per run to $work/PREFIX.list and one elapsed
#   wall-clock-seconds value per run to $work/PREFIX.times. Exits the
#   whole script nonzero (per run_aniso_case.sh's own stderr) if any run
#   fails.
run_group() {
  _rg_scenario=$1; _rg_eps=$2; _rg_r=$3; _rg_mode=$4
  _rg_kstart=$5; _rg_kend=$6; _rg_prefix=$7

  _rg_model=$(scenario_model "$_rg_scenario" "$_rg_eps" "$_rg_mode" "$_rg_r")
  _rg_phonons=$(phonons_for_scenario "$_rg_scenario")

  : > "$work/$_rg_prefix.list"
  : > "$work/$_rg_prefix.times"

  _rg_k=$_rg_kstart
  while [ "$_rg_k" -le "$_rg_kend" ]; do
    _rg_seed=$(seed_for "$_rg_k")
    _rg_outdir=$(mktemp -d "$work/run.XXXXXX")
    _rg_start=$(date +%s.%N)
    if ! R3D_EXE="$R3D_EXE_CURRENT" TOA_DEGREE="$CMP_DEGREE" SEED="$_rg_seed" \
          "$ROOT/tests/run_aniso_case.sh" "$_rg_scenario" "$CMP_THREADS" \
          "$_rg_outdir" "$_rg_phonons" "$_rg_model" \
          >"$_rg_outdir/harness-stdout.txt" 2>"$_rg_outdir/harness-stderr.txt"; then
      echo "ERROR: run failed (scenario=$_rg_scenario mode=$_rg_mode r=$_rg_r" \
           "eps=$_rg_eps k=$_rg_k seed=$_rg_seed)" >&2
      dump_failure_logs "$_rg_outdir"
      exit 1
    fi
    _rg_end=$(date +%s.%N)
    _rg_elapsed=$(awk "BEGIN{printf \"%.3f\", $_rg_end - $_rg_start}")
    echo "$_rg_elapsed" >> "$work/$_rg_prefix.times"

    _rg_tsv="$_rg_outdir/energy.tsv"
    extract_energy "$_rg_outdir/seis_traces_asc.dat" "$_rg_tsv"
    echo "$_rg_tsv" >> "$work/$_rg_prefix.list"

    _rg_k=$((_rg_k + 1))
  done
}

# window_stat GROUP1_LIST GROUP2_LIST K NPERM SEED
#
#   Prints "Tobs pvalue nkept" on stdout, computed as described in the
#   header comment above. GROUP1_LIST/GROUP2_LIST are files with K
#   energy.tsv paths each, one per line, in replicate order (as written
#   by run_group).
window_stat() {
  _ws_g1=$1; _ws_g2=$2; _ws_k=$3; _ws_nperm=$4; _ws_seed=$5
  awk -v K="$_ws_k" -v NPERM="$_ws_nperm" -v SEED="$_ws_seed" \
      -v G1LIST="$_ws_g1" -v G2LIST="$_ws_g2" '
    # Fixed-seed LCG (Numerical-Recipes constants): _state is the only
    # thing that persists between calls, seeded once from SEED below.
    function lcg_next() {
      _state = (1664525.0*_state + 1013904223.0)
      _state = _state - int(_state/4294967296.0)*4294967296.0
      return _state
    }
    function lcg_uniform() { return lcg_next()/4294967296.0 }

    # compute_T(grpOf): mean(z^2) over the shared "kept" window list,
    # grouping the 2K slots 1..total by grpOf[slot] (1 or 2). All
    # variables after grpOf are awk'"'"'s local-variable idiom (extra,
    # un-passed formal parameters), not globals.
    function compute_T(grpOf,    i,s,key,v,n1,n2,sum1,sum2,sumsq1,sumsq2,
                                  m1,m2,var1,var2,denom,z,sumz2,cnt) {
      sumz2 = 0; cnt = 0
      for (i=1; i<=nkept; i++) {
        key = kept[i]
        n1=0; n2=0; sum1=0; sum2=0; sumsq1=0; sumsq2=0
        for (s=1; s<=total; s++) {
          v = val[key,s]
          if (grpOf[s]==1) { n1++; sum1+=v; sumsq1+=v*v }
          else              { n2++; sum2+=v; sumsq2+=v*v }
        }
        m1 = sum1/n1; m2 = sum2/n2
        var1 = (sumsq1 - n1*m1*m1)/(n1-1)
        var2 = (sumsq2 - n2*m2*m2)/(n2-1)
        denom = sqrt(var1/n1 + var2/n2)
        if (denom > 0) {
          z = (m1-m2)/denom
          sumz2 += z*z
          cnt++
        }
      }
      if (cnt==0) {
        print "window_stat: no window had nonzero combined variance" > "/dev/stderr"
        exit 1
      }
      return sumz2/cnt
    }

    BEGIN {
      total = 2*K

      nfiles = 0
      while ((getline line < G1LIST) > 0) { nfiles++; files[nfiles]=line; slotgrp[nfiles]=1 }
      close(G1LIST)
      g1count = nfiles
      while ((getline line < G2LIST) > 0) { nfiles++; files[nfiles]=line; slotgrp[nfiles]=2 }
      close(G2LIST)
      if (g1count != K || (nfiles-g1count) != K) {
        printf "window_stat: expected %d files per group, got %d and %d\n", \
               K, g1count, nfiles-g1count > "/dev/stderr"
        exit 1
      }

      g1i = 0; g2i = 0
      for (fi=1; fi<=nfiles; fi++) {
        if (slotgrp[fi]==1) { g1i++; slot = g1i } else { g2i++; slot = K+g2i }
        while ((getline line < files[fi]) > 0) {
          split(line, f, "\t")
          key = f[1] SUBSEP f[2]
          val[key,slot] = f[3]+0
          seenkey[key] = 1
        }
        close(files[fi])
      }

      nkept = 0
      for (key in seenkey) {
        ok = 1
        for (s=1; s<=total; s++) { if (val[key,s] <= 0) { ok=0; break } }
        if (ok) { nkept++; kept[nkept] = key }
      }
      if (nkept == 0) {
        print "window_stat: no window had positive energy in every replicate of both groups" > "/dev/stderr"
        exit 1
      }

      for (s=1; s<=total; s++) baseGrp[s] = (s<=K) ? 1 : 2
      Tobs = compute_T(baseGrp)

      _state = SEED
      count_ge = 0
      for (perm=1; perm<=NPERM; perm++) {
        for (s=1; s<=total; s++) idx[s]=s
        for (s=total; s>=2; s--) {
          rpos = int(lcg_uniform()*s) + 1
          tmp=idx[s]; idx[s]=idx[rpos]; idx[rpos]=tmp
        }
        for (s=1; s<=total; s++) permGrp[idx[s]] = (s<=K) ? 1 : 2
        Tp = compute_T(permGrp)
        if (Tp >= Tobs) count_ge++
      }
      pval = count_ge / NPERM

      printf "%.6f %.6f %d\n", Tobs, pval, nkept
    }
  '
}

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

echo ""
echo "# Axisymmetric vs. Near-Axisymmetric-General Scenario Validation"
echo ""
echo "Degree: $CMP_DEGREE, Threads: $CMP_THREADS, Replicates (K): $CMP_REPLICATES, Permutations: $CMP_PERMUTATIONS"
echo ""
echo "| Scenario | eps | r | windows | p(A,G) | p(A,A') | median s (A) | median s (G) | G/A time ratio |"
echo "|---|---|---:|---:|---:|---:|---:|---:|---:|"

for _cag_cfg in $CMP_CONFIGS; do
  _cag_fields=$(cmp_config_fields "$_cag_cfg")
  set -- $_cag_fields
  _cag_scenario=$1; _cag_eps=$2; _cag_r=$3

  echo "Config: scenario=$_cag_scenario eps=$_cag_eps r=$_cag_r (K=$K)" >&2

  echo "  Running group A (axi, k=1..$K)..." >&2
  run_group "$_cag_scenario" "$_cag_eps" "$_cag_r" axi 1 "$K" A
  echo "  Running group G (gen, k=$((K+1))..$((2*K)))..." >&2
  run_group "$_cag_scenario" "$_cag_eps" "$_cag_r" gen "$((K+1))" "$((2*K))" G
  echo "  Running group A' (axi, k=$((2*K+1))..$((3*K)))..." >&2
  run_group "$_cag_scenario" "$_cag_eps" "$_cag_r" axi "$((2*K+1))" "$((3*K))" AP

  echo "  Computing A-vs-G statistic..." >&2
  _cag_result_ag=$(window_stat "$work/A.list" "$work/G.list" "$K" "$CMP_PERMUTATIONS" "$LCG_SEED")
                    # Captured into a plain assignment (not used directly
                    # as another command's argument) so that "set -e"
                    # actually aborts the script if window_stat's awk
                    # program hit one of its own `exit 1` error paths --
                    # "set -- $(window_stat ...)" would not: "set --"
                    # always itself "succeeds", so set -e would not
                    # trigger on it even if the substitution failed.
  set -- $_cag_result_ag
  p_AG=$2; windows_AG=$3

  echo "  Computing A-vs-A' statistic..." >&2
  _cag_result_aa=$(window_stat "$work/A.list" "$work/AP.list" "$K" "$CMP_PERMUTATIONS" "$LCG_SEED")
  set -- $_cag_result_aa
  p_AA=$2; windows_AA=$3

  if [ "$windows_AA" != "$windows_AG" ]; then
    echo "  NOTE: A-vs-A' kept $windows_AA windows, vs $windows_AG for A-vs-G" \
         "(each comparison filters positivity over its own 2K runs)" >&2
  fi

  med_A=$(median_of "$work/A.times")
  med_G=$(median_of "$work/G.times")
  ratio=$(awk -v a="$med_A" -v g="$med_G" 'BEGIN{printf "%.3f", g/a}')

  if awk -v p="$p_AG" 'BEGIN{exit !(p<0.01)}'; then
    echo "  GATE: p(A,G)=$p_AG < 0.01 for scenario=$_cag_scenario eps=$_cag_eps r=$_cag_r" >&2
    GATE_FAIL=1
  fi
  if awk -v p="$p_AA" 'BEGIN{exit !(p<0.01)}'; then
    echo "  GATE: p(A,A')=$p_AA < 0.01 for scenario=$_cag_scenario eps=$_cag_eps r=$_cag_r" >&2
    GATE_FAIL=1
  fi

  printf '| %s | %s | %s | %s | %.4f | %.4f | %s | %s | %sx |\n' \
    "$_cag_scenario" "$_cag_eps" "$_cag_r" "$windows_AG" "$p_AG" "$p_AA" \
    "$med_A" "$med_G" "$ratio"

  # Free this config's run directories before the next one, so a large
  # CMP_CONFIGS sweep doesn't accumulate every run's output on disk at
  # once (each run's own files are no longer needed once its energy.tsv
  # has been extracted and its elapsed time recorded).
  rm -rf "$work"/run.* "$work"/A.list "$work"/G.list "$work"/AP.list \
         "$work"/A.times "$work"/G.times "$work"/AP.times
done

echo "" >&2
echo "Comparison complete" >&2

exit "$GATE_FAIL"
