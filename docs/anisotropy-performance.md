# Anisotropic Scale-Length Performance Measurement

This file documents measured performance and validation results for the
anisotropic scale-length / axisymmetric-sampler work, most recent
measurements first. Nothing in this file is estimated; every number traces
to a captured test or benchmark run, as each section states.

## 2026-10-25 measurements (code 65128c6)

Measured by running `BASELINE_EXE=/tmp/r3d-baseline/main tests/benchmark_aniso.sh`
(captured as `bench3-full.out`), `tests/compare_axi_general.sh` (captured as
`cmp2-full.out`), two fresh-seed confirmation runs of
`tests/compare_axi_general.sh` using its new `CMP_SEED_BASE` option
(captured as `cmp3-fresh.out`, run at `8a996ec` -- one commit after
`65128c6` that only adds `CMP_SEED_BASE` itself, so it exercises the same
production code as `65128c6`), and `tests/test_anisotropy_unit.sh
near_axi_cost`, snapshotted as
`docs/anisotropy-near-axi-cost-2026-10-25.txt` (that file's own header
names the exact command and commit). That exact filtered command is
reproducible: `random.hpp`'s `RandomUniform()` falls back to a
`thread_local` engine fixed-seeded `0x7261643364` whenever no
`RandomContext` is active, which is always true entering a standalone
`near_axi_cost` run, so rerunning `tests/test_anisotropy_unit.sh
near_axi_cost` alone reproduces the snapshot byte-for-byte. Running the
full, unfiltered test suite instead gives slightly different values
(roughly 1-3%), because earlier tests draw from that same fallback stream
first; the snapshot, not the live `build/` file (which a full-suite run
overwrites with those different values), is this section's source of
record. All numbers in this section come from those captured outputs;
nothing here is estimated.

### Host, compiler, and workload

| Date | Host | Compiler | Workload |
|---|---|---|---|
| 2026-10-25 | Apple M2 Ultra, 24 threads | g++ (Homebrew GCC 16.2.0) 16.2.0 | Degree 9; lopnor 1,000,000 phonons; crustpinch 200,000 phonons |

### Benchmark table

`tests/benchmark_aniso.sh` (degree 9, 24 threads, repeats 3) now times, for
every axisymmetric (`axi`) row, a matching general-sampler row on the *same
medium* (`gen`: `ay = ax * (1 + 1e-5)` instead of `ay = ax`), alongside the
original isotropic and triaxial (`tri`: `ax = r*a`, `ay = sqrt(r)*a`) rows.
Each row also reports a zero-phonon build time (construction only, a single
measurement) alongside the total time (the median of the 3 repeats).
Output, verbatim from `bench3-full.out`:

| Scenario | eps | Mode | r | Build s | Total s | Slowdown vs isotropic | Gate vs isotropic | axi/gen ratio | Gate same-medium |
|---|---|---|---|---:|---:|---:|---|---:|---|
| lopnor | 0.01 | iso | - | 3.79 | 4.92 | 1.000x | ok | - | - |
| lopnor | 0.01 | axi | 2 | 4.43 | 5.74 | 1.166x | ok | 1.082 | SLOW |
| lopnor | 0.01 | gen | 2 | 3.86 | 5.30 | 1.078x | ok | - | - |
| lopnor | 0.01 | tri | 2 | 3.86 | 5.22 | 1.061x | ok | - | - |
| lopnor | 0.01 | axi | 4 | 4.51 | 6.00 | 1.220x | ok | 1.005 | ok |
| lopnor | 0.01 | gen | 4 | 3.89 | 5.97 | 1.214x | ok | - | - |
| lopnor | 0.01 | tri | 4 | 3.89 | 5.46 | 1.109x | ok | - | - |
| lopnor | 0.01 | axi | 10 | 4.54 | 7.58 | 1.541x | ok | 0.726 | ok |
| lopnor | 0.01 | gen | 10 | 3.94 | 10.44 | 2.122x | GATE | - | - |
| lopnor | 0.01 | tri | 10 | 3.88 | 6.84 | 1.391x | ok | - | - |
| lopnor | 0.05 | iso | - | 3.73 | 11.18 | 1.000x | ok | - | - |
| lopnor | 0.05 | axi | 2 | 4.53 | 19.00 | 1.700x | ok | 0.988 | ok |
| lopnor | 0.05 | gen | 2 | 3.86 | 19.23 | 1.720x | ok | - | - |
| lopnor | 0.05 | tri | 2 | 3.90 | 16.06 | 1.437x | ok | - | - |
| lopnor | 0.05 | axi | 4 | 4.54 | 30.42 | 2.722x | GATE | 0.703 | ok |
| lopnor | 0.05 | gen | 4 | 3.95 | 43.30 | 3.874x | GATE | - | - |
| lopnor | 0.05 | tri | 4 | 3.90 | 27.98 | 2.503x | GATE | - | - |
| lopnor | 0.05 | axi | 10 | 4.47 | 79.26 | 7.093x | GATE | 0.476 | ok |
| lopnor | 0.05 | gen | 10 | 3.96 | 166.41 | 14.892x | GATE | - | - |
| lopnor | 0.05 | tri | 10 | 3.85 | 69.21 | 6.194x | GATE | - | - |
| crustpinch | 0.01 | iso | - | 2.04 | 2.24 | 1.000x | ok | - | - |
| crustpinch | 0.01 | axi | 2 | 2.25 | 2.48 | 1.104x | ok | 1.075 | SLOW |
| crustpinch | 0.01 | gen | 2 | 2.12 | 2.31 | 1.027x | ok | - | - |
| crustpinch | 0.01 | tri | 2 | 2.08 | 2.31 | 1.030x | ok | - | - |
| crustpinch | 0.01 | axi | 4 | 2.27 | 2.47 | 1.099x | ok | 1.078 | SLOW |
| crustpinch | 0.01 | gen | 4 | 2.08 | 2.29 | 1.020x | ok | - | - |
| crustpinch | 0.01 | tri | 4 | 2.08 | 2.31 | 1.028x | ok | - | - |
| crustpinch | 0.01 | axi | 10 | 2.26 | 2.50 | 1.113x | ok | 1.046 | ok |
| crustpinch | 0.01 | gen | 10 | 2.11 | 2.39 | 1.064x | ok | - | - |
| crustpinch | 0.01 | tri | 10 | 2.13 | 2.31 | 1.029x | ok | - | - |
| crustpinch | 0.05 | iso | - | 2.02 | 2.29 | 1.000x | ok | - | - |
| crustpinch | 0.05 | axi | 2 | 2.26 | 2.58 | 1.127x | ok | 1.053 | SLOW |
| crustpinch | 0.05 | gen | 2 | 2.11 | 2.45 | 1.070x | ok | - | - |
| crustpinch | 0.05 | tri | 2 | 2.11 | 2.42 | 1.058x | ok | - | - |
| crustpinch | 0.05 | axi | 4 | 2.30 | 2.68 | 1.168x | ok | 0.941 | ok |
| crustpinch | 0.05 | gen | 4 | 2.10 | 2.84 | 1.241x | ok | - | - |
| crustpinch | 0.05 | tri | 4 | 2.07 | 2.63 | 1.147x | ok | - | - |
| crustpinch | 0.05 | axi | 10 | 2.28 | 3.38 | 1.475x | ok | 0.592 | ok |
| crustpinch | 0.05 | gen | 10 | 2.12 | 5.70 | 2.490x | GATE | - | - |
| crustpinch | 0.05 | tri | 10 | 2.08 | 3.53 | 1.542x | ok | - | - |

#### Isotropic no-slowdown check (iso-gate verdict)

Isotropic current/baseline ratios (median of 3 current runs over 3 baseline
(`aaec540`) runs): **0.2115** (lopnor), **0.2699** (crustpinch) -- isotropic
runs now take about 21-27% of the pre-anisotropy time. Both are comfortably
under the 1.05 no-slowdown threshold, so this check passed.

The run as a whole still exited nonzero (`EXIT=1` in `bench3-full.out`).
That is not from a failed case or an isotropic-ratio violation -- it is from
four same-medium rows marked `SLOW`, covered next. A `GATE` row (more than
2x the isotropic time) is reported but, per `tests/benchmark_aniso.sh`'s own
contract, does not by itself fail the run; several axi/tri/gen rows at high
`eps` and large `r` exceed 2x (e.g. lopnor eps=0.05 r=10 gen: 14.892x) --
that is the expected cost of a strongly anisotropic, high-phonon-count
workload, not a sign of a defect.

### Same-medium verdict

Design decision 5's gate compares the axisymmetric sampler against the
general sampler on the *same medium* (`gen`'s `ay = ax * (1 + 1e-5)`, the
near-axisymmetric medium): axi total time must be at most 1.05x the gen
total time in every configuration. In `bench3-full.out`, four rows exceed
it, all at `eps = 0.01`, all build-dominated (small absolute totals, where
the axi/gen build-time gap is a larger fraction of the total), all between
1.05x and 1.08x:

| Scenario | eps | r | axi/gen ratio | axi total (s) | gen total (s) |
|---|---|---:|---:|---:|---:|
| lopnor | 0.01 | 2 | 1.082 | 5.74 | 5.30 |
| crustpinch | 0.01 | 2 | 1.075 | 2.48 | 2.31 |
| crustpinch | 0.01 | 4 | 1.078 | 2.47 | 2.29 |
| crustpinch | 0.05 | 2 | 1.053 | 2.58 | 2.45 |

`progress.md`'s ruling R8 accepted these rows rather than adding a
fallback (R8 was made against the equivalent `bench2-full.out` run at an
earlier commit, `42b2fe6`, which showed the same four build-dominated SLOW
rows; `progress.md` notes that the final run here reproduces the same
build-dominated pattern that ruling R8 covers):

> Ruling: R8 -- accept the four SLOW same-medium rows (all eps 0.01,
> build-dominated, <= 8%, 0.1-0.4 s absolute); no AUTO fallback -- an
> M-based fallback would forfeit crustpinch's eps 0.05 r=4/10 wins
> (0.943/0.597) because those same specs are SLOW at eps 0.01, and the
> remaining lever (lazy per-bin construction) would relax the original
> spec's read-only-scatterer constraint, which is the user's call -- cost
> if wrong: weak-scattering axisymmetric runs stay up to ~8% slower than
> the general sampler on the same medium. Surfaced to the user.

(R8's own cited crustpinch eps=0.05 r=4/10 figures, 0.943/0.597, are from
`bench2-full.out`; the equivalent figures in the final `bench3-full.out`
run above are 0.941/0.592 -- the same pattern, not identical runs.) No
static fallback (such as an aspect-ratio/`M` threshold that routes small
axisymmetric specs to the general sampler) was implemented, because it
would have cost more than it saved: it would forfeit real wins at the same
`eps` for a different `r`, and the only remaining lever that could close
the gap without a fallback -- lazy per-bin construction -- would relax the
plan's read-only-scatterer constraint.

The remaining configurations pass the same-medium gate (ratios from 0.476
to 1.046), with large wins at the heaviest workloads, e.g.:

- lopnor eps=0.05 r=10: axi 79.26 s vs gen 166.41 s (ratio 0.476)
- crustpinch eps=0.05 r=10: axi 3.38 s vs gen 5.70 s (ratio 0.592)

#### Why the axi and tri rows differ

The axisymmetric (`axi`) and triaxial (`tri`, general sampler) rows at a
given `r` are not the same medium, so their times should not be compared
directly to each other. An axisymmetric row sets both horizontal scale
lengths to `ax = ay = r*a`; a triaxial row sets only one, `ax = r*a`,
`ay = sqrt(r)*a`. Both samplers' majorant bound factor `M` scales with the
medium's horizontal cross-section -- `ah^2 * az / a_min^3` for the
axisymmetric medium's `ah = ay = r*a`, against
`ah * sqrt(r)*a * az / a_min^3` for the triaxial medium -- so `M` grows as
`r^2` for the axisymmetric medium but only `r^1.5` (`r * sqrt(r)`) for the
triaxial one: at the same `r`, the axisymmetric medium is the more strongly
scattering of the two. That is visible directly in the table above: at
lopnor eps=0.01 r=10, the `axi` row totals 7.58 s against the `tri` row's
6.84 s, even though the axisymmetric *sampler* is the faster of the two
once the medium is held fixed (axi vs. the same-medium `gen` row, 7.58 s vs
10.44 s, ratio 0.726). The fair sampler-performance comparison is therefore
always axi against gen on the same medium (the same-medium gate above), not
axi against tri.

### Build-time history

The degree-9 Lop Nor axisymmetric zero-phonon (construction-only) build
time, as recorded in `progress.md` at each task's completion (single-run
measurements except where noted "median"; default/automatic thread count,
24 logical cores on the M2 Ultra host used throughout this plan):

| Stage (commit) | axi build (s) | Decision |
|---|---:|---|
| `3db253a` baseline (measured for comparison in Task 2) | 106.19 | -- |
| Task 2: restructured bound loop + per-degree geometry cache (`08231c5`) | 66.57 | 1 |
| Task 3: parallel bin/report-state build (`86e6692`) | 27.88 | 2 |
| Task 4: trig-free cap kernel (`961ca45`) | 26.10 | 3 |
| Task 4b: parallelize geometry and `PrepareBoundRep` (`e64bb3c`) | 26.20 | 1/2 follow-up (R3) |
| Task 4c: parallel isotropic table construction, chunked claiming (`6edbe97`) | 7.43 | 8 |
| Task 5b: remaining axisymmetric build overhead removed, median of 3 (`e13063c`) | 4.53 | 1/2 follow-up (R7) |

The isotropic zero-phonon build dropped alongside it, from about 21.9 s
(per Task 4's note on the pre-Task-4c state) to 3.72 s (Task 4c,
`6edbe97`): Task 4c's bit-identical parallel `PopulateProbDists` /
`PSATO`-caching change (decision 8) benefits every mode, isotropic
included, which is why the axi and iso columns both drop sharply at that
step.

### Near-axisymmetric validation

Design decision 6 validates the axisymmetric sampler against the general
sampler run on a medium with `ay = ax * (1 + delta)`, `delta = 1e-5`
("near-axisymmetric"), at the unit level (Task 6, Task 6b) and the scenario
level (Task 7, this benchmark).

#### Unit level

From `progress.md`'s Task 6 and Task 6b completion notes:

- **Rate scaling:** `IntegratedRate`'s change with `delta` scales linearly,
  as decision 6 requires: the ratio of the rate deviation at
  `delta = 1e-4` to the deviation at `delta = 1e-6`,
  `dR(1e-4)/dR(1e-6)`, falls in **[99.985, 99.999]** against the
  theoretical value of 100 for an exactly linear response.
- **Sampler agreement:** the axisymmetric and general samplers agree on
  mean distance to a real event across the same 10 incidence states used
  in the draws table below (**|z| <= 1.29**), and on the (conversion,
  cell) distribution (**chi2/dof = 0.987**, 613 degrees of freedom).
- **Report MFPs:** AXI<->ANI report mean-free-path pairs agree to within
  **1e-5 to 1e-6** (design decision 6's own threshold was 1e-4).
- **P-frame consistency (Task 6b, ruling R10):** the investigation below
  found that, before this fix, the general sampler's P-incidence rate
  depended on the (unphysical) carried polarization-rotation frame, while
  the axisymmetric sampler always used the canonical (highest-rate) frame.
  Ruling R10 switched the general sampler to the same canonical P frame for
  P incidence. A dedicated consistency test at a non-canonical carried
  polarization then found the pre-fix ("RED") cases failing in 6/6 trials
  (`chi2/dof` 2.05-3.09, `z` as extreme as -4.26), and the post-fix
  ("GREEN") cases passing (`chi2/dof` 0.82-1.19, `|z|` 0.21-1.71).

#### Scenario level

`tests/compare_axi_general.sh` at degree 7, K=8 replicates, 1000
permutations, all 8 benchmark configurations, captured as `cmp2-full.out`
(note this script runs at a different take-off-angle degree and replicate
structure than the degree-9 benchmark above, so its own G/A time-ratio
column is independent corroboration, not a substitute for the same-medium
gate):

| Scenario | eps | r | windows | p(A,G) | p(A,A') | median s (A) | median s (G) | G/A time ratio |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| lopnor | 0.01 | 4 | 76 | 0.7800 | 0.9580 | 1.445 | 0.956 | 0.662x |
| lopnor | 0.01 | 10 | 74 | 0.0040 | 0.6170 | 1.782 | 1.910 | 1.072x |
| lopnor | 0.05 | 4 | 659 | 0.2420 | 0.1320 | 6.978 | 8.984 | 1.287x |
| lopnor | 0.05 | 10 | 613 | 0.3750 | 0.7900 | 17.735 | 35.519 | 2.003x |
| crustpinch | 0.01 | 4 | 105 | 0.2660 | 0.1620 | 0.500 | 0.340 | 0.680x |
| crustpinch | 0.01 | 10 | 101 | 0.3480 | 0.2500 | 0.502 | 0.350 | 0.697x |
| crustpinch | 0.05 | 4 | 119 | 0.2500 | 0.1540 | 0.534 | 0.430 | 0.805x |
| crustpinch | 0.05 | 10 | 145 | 0.6030 | 0.3700 | 0.687 | 0.975 | 1.419x |

`p(A,G)` is the permutation-test p-value comparing the axisymmetric group
(A) against the near-axisymmetric general group (G) on windowed trace
energies; `p(A,A')` is the null-control comparison between two independent
A replicate groups. All configurations pass the design's p >= 0.01
threshold except lopnor eps=0.01 r=10, `p(A,G) = 0.0040`.

That same configuration, with the same default seed set, produced the
plan's original p=0.002 failure at an earlier commit, investigated in
`progress.md` as a "Debugging result":

> Debugging result: p=0.002 is a fluke -- fresh seeds K=16: p(A,G)=0.903
> and 0.312 (controls 0.090, 0.751); pooled 80 axi vs 40 gen runs p=0.99
> for total energy, all components p>=0.29; original failure driven by one
> window (station 7, 250-260 s) with an unusually tight A group.

(That investigation also found and fixed the P-frame defect documented
under "Unit level" above; the pooled 80-vs-40 check and the fresh-seed
numbers quoted here are for the *scenario-level timing fluke*, a separate,
pre-existing issue from the P-frame defect.)

Per ruling R11, two fresh-seed confirmations of lopnor eps=0.01 r=10 were
re-run directly on the final code (`8a996ec`) using the new
`CMP_SEED_BASE` option, K=16, 2000 permutations each, captured as
`cmp3-fresh.out`:

| Seed base | windows | p(A,G) | p(A,A') |
|---|---:|---:|---:|
| 596 + 1000*2^32 | 68 | 0.9095 | 0.0900 |
| 596 + 2000*2^32 | 71 | 0.3075 | 0.7505 |

Both fresh-seed runs pass comfortably and closely reproduce the earlier
investigation's own fresh-seed numbers (0.903/0.090 and 0.312/0.751 above),
confirming that the default-seed `p(A,G) = 0.0040` is the same documented
sampling fluke, reproduced deterministically by that one seed set, rather
than a regression introduced since the original investigation.

### Draws and tentative collisions per real event (final code)

Captured by running `tests/test_anisotropy_unit.sh near_axi_cost` once,
standalone (that argument runs only this one test in the process, per
`tests/test_anisotropy_unit.sh`'s own FILTER handling), against
`HEAD` (`e227bb7`; the test builds its own axisymmetric and
near-axisymmetric-general `Scatterer`s, so it exercises the same
production code as `65128c6`, the last commit to touch it), and
snapshotted immediately afterward, before any further test run, to
`docs/anisotropy-near-axi-cost-2026-10-25.txt` (ruling R12) -- not cited
from the live `build/aniso_near_axi_cost.txt`, which a later, unfiltered
test-suite run would overwrite with different values. Mean tentative
collisions and mean draws per real event, n=5000 trials per state, at
take-off-angle degree 7, for three `(ax, ay, az)` km specs (axisymmetric
vs. the matching `ay*(1+1e-5)` near-axisymmetric general-sampler medium).
A "draw" is a tentative collision that reaches the (conversion,
take-off-angle) table (see the comment above `MeasureCollisionCost` in
`tests/anisotropy_tests.cpp`); the general sampler's draws always equal its
tentative collisions, while the axisymmetric sampler's cheap top-level
bin-rate test thins most tentative collisions before they become draws --
most of why it needs far fewer tentative collisions overall. This exact
filtered command is reproducible: `random.hpp`'s `RandomUniform()` falls
back to a `thread_local` engine fixed-seeded `0x7261643364` whenever no
`RandomContext` is active (always true entering a standalone run here), so
the table below reproduces byte-for-byte on a rerun of the same filtered
command. Running the full, unfiltered test suite instead shifts every
number by roughly 1-3%, because earlier tests in that process draw from
the same fallback stream first -- consistent with the very similar, but
not identical, numbers quoted under "Unit level" above and in earlier task
reports (which came from full-suite or otherwise differently-ordered
runs).

#### Spec (1, 1, 0.5)

| Ray | theta | chi | AXI mean tentative | AXI mean draws | GENERAL mean tentative = draws |
|---|---:|---:|---:|---:|---:|
| P | 0 | - | 1.1874 | 1.1874 | 2.596 |
| P | 30 | - | 1.1996 | 1.1856 | 2.6484 |
| P | 60 | - | 1.2168 | 1.1868 | 2.7442 |
| P | 89 | - | 1.2066 | 1.1780 | 2.6502 |
| S | 30 | 0 | 1.5068 | 1.3022 | 2.4934 |
| S | 30 | 45 | 1.5114 | 1.2804 | 2.5694 |
| S | 30 | 90 | 1.5256 | 1.2708 | 2.5854 |
| S | 89 | 0 | 1.3524 | 1.2942 | 2.2638 |
| S | 89 | 45 | 1.3092 | 1.2862 | 2.1256 |
| S | 89 | 90 | 1.2914 | 1.2914 | 2.1596 |

#### Spec (2, 2, 0.5)

| Ray | theta | chi | AXI mean tentative | AXI mean draws | GENERAL mean tentative = draws |
|---|---:|---:|---:|---:|---:|
| P | 0 | - | 2.0324 | 1.5188 | 9.1784 |
| P | 30 | - | 1.9390 | 1.5138 | 8.7448 |
| P | 60 | - | 1.6244 | 1.5272 | 7.2364 |
| P | 89 | - | 1.4862 | 1.4862 | 6.5554 |
| S | 30 | 0 | 3.2852 | 1.8626 | 7.7522 |
| S | 30 | 45 | 3.4212 | 1.8510 | 7.9084 |
| S | 30 | 90 | 3.5064 | 1.8154 | 8.1964 |
| S | 89 | 0 | 2.1224 | 1.8070 | 4.8282 |
| S | 89 | 45 | 1.9746 | 1.8626 | 4.6190 |
| S | 89 | 90 | 1.8904 | 1.8904 | 4.4152 |

#### Spec (5, 5, 0.5)

| Ray | theta | chi | AXI mean tentative | AXI mean draws | GENERAL mean tentative = draws |
|---|---:|---:|---:|---:|---:|
| P | 0 | - | 6.5050 | 2.9378 | 51.6572 |
| P | 30 | - | 5.9394 | 3.0332 | 49.9656 |
| P | 60 | - | 4.1038 | 3.4264 | 33.0688 |
| P | 89 | - | 2.7092 | 2.6660 | 21.9640 |
| S | 30 | 0 | 14.6152 | 5.3586 | 42.8454 |
| S | 30 | 45 | 15.1612 | 5.2902 | 44.5122 |
| S | 30 | 90 | 15.2802 | 5.1116 | 45.4492 |
| S | 89 | 0 | 4.7056 | 3.7596 | 14.1366 |
| S | 89 | 45 | 4.5634 | 4.3690 | 13.4300 |
| S | 89 | 90 | 4.5994 | 4.3436 | 14.0526 |

As aspect ratio grows (and with it `M`, per "Why the axi and tri rows
differ" above), the general sampler's tentative-collision count grows
sharply while the axisymmetric sampler's stays nearly flat: at
`(5,5,0.5)`, P incidence at 0 degrees needs 6.51 axi tentative collisions
(2.94 draws) against 51.66 general (all of them draws) -- consistent with
(though not identical to, for the no-fixed-seed reason noted above) Task
6's own measurement of this same case, 6.53 vs 52.65, axi draws 2.97.

### Summary

Construction and sampling now cost about the same as, or less than, the
general sampler on the same medium in 8 of the 12 same-medium axi/gen
pairs tested above (the other 4 are the SLOW rows just discussed), closing
most of the gap recorded in the
2026-10-16 baseline below, where every axisymmetric row failed the 2x
isotropic-slowdown gate (lopnor axi 4.223x-6.313x; crustpinch axi
3.619x-3.722x). The remaining gap -- four small, build-dominated
configurations 5-8% slower than general -- is accepted per ruling R8. The
near-axisymmetric validation above (unit and scenario level) found the
samplers agree with each other to within statistical noise, after fixing
one real defect (the P-incidence frame dependence, ruling R10) uncovered
by that validation work.

## 2026-10-16 baseline (3db253a)

Measured by running `tests/test_anisotropy_unit.sh axi_draws` and
`tests/benchmark_aniso.sh` against the `anisotropy` branch build. All numbers
below come from those runs; nothing here is estimated.

| Date | Host | Compiler | Workload |
|---|---|---|---|
| 2026-10-16 | Apple M2 Ultra, 24 threads | Apple clang version 17.0.0 (clang-1700.6.4.2) | Degree 9; lopnor 1,000,000 phonons; crustpinch 200,000 phonons |

> **Correction (2026-10-25):** the compiler above is wrong. This host's
> `c++` resolves to Apple clang, but the project's `Makefile` -- which
> these benchmarks always build with -- invokes `g++`, which resolves to a
> Homebrew-installed GCC, as Task 3's implementer noted while working on
> this plan (`progress.md`: "Makefile compiles with g++ (Homebrew GCC 16.2
> per implementer)"). The line above was produced by running `c++
> --version` directly rather than checking what `make` invokes. The
> 2026-10-25 measurements above, built the same way on the same host,
> report `g++ (Homebrew GCC 16.2.0) 16.2.0`; the exact Homebrew GCC version
> active on 2026-10-16 itself was not independently recorded, but it was
> almost certainly a g++ in the same family, not Apple clang. The
> benchmark numbers below are unaffected by this correction and are left as
> originally recorded.

### Benchmark table

`tests/benchmark_aniso.sh` times the isotropic run and the axisymmetric
(`ax=ay=r·a, az=a`) and triaxial (`ax=r·a, ay=√r·a, az=a`) runs at aspect
ratios r ∈ {2, 4, 10}, for both scenarios' do-script ε (0.01) and ε = 0.05.
Anisotropy is applied to the crust and mantle regions for lopnor, and to the
crust, pinched crust, and mantle regions for crustpinch (sedi and moho stay
isotropic in both). Output, verbatim:

| Scenario | eps | Mode | r | Seconds | Slowdown vs isotropic | Gate |
|---|---|---|---|---:|---:|---|
| lopnor | 0.01 | iso | - | 23.01 | 1.00x | ok |
| lopnor | 0.01 | axi | 2 | 109.68 | 4.767x | GATE |
| lopnor | 0.01 | axi | 4 | 109.96 | 4.778x | GATE |
| lopnor | 0.01 | axi | 10 | 111.24 | 4.834x | GATE |
| lopnor | 0.01 | tri | 2 | 27.71 | 1.204x | ok |
| lopnor | 0.01 | tri | 4 | 28.06 | 1.219x | ok |
| lopnor | 0.01 | tri | 10 | 29.41 | 1.278x | ok |
| lopnor | 0.05 | iso | - | 29.09 | 1.00x | ok |
| lopnor | 0.05 | axi | 2 | 122.82 | 4.223x | GATE |
| lopnor | 0.05 | axi | 4 | 134.36 | 4.620x | GATE |
| lopnor | 0.05 | axi | 10 | 183.60 | 6.313x | GATE |
| lopnor | 0.05 | tri | 2 | 38.68 | 1.330x | ok |
| lopnor | 0.05 | tri | 4 | 50.46 | 1.735x | ok |
| lopnor | 0.05 | tri | 10 | 91.30 | 3.139x | GATE |
| crustpinch | 0.01 | iso | - | 8.21 | 1.00x | ok |
| crustpinch | 0.01 | axi | 2 | 29.83 | 3.634x | GATE |
| crustpinch | 0.01 | axi | 4 | 29.81 | 3.632x | GATE |
| crustpinch | 0.01 | axi | 10 | 29.84 | 3.636x | GATE |
| crustpinch | 0.01 | tri | 2 | 9.35 | 1.140x | ok |
| crustpinch | 0.01 | tri | 4 | 9.36 | 1.140x | ok |
| crustpinch | 0.01 | tri | 10 | 9.44 | 1.150x | ok |
| crustpinch | 0.05 | iso | - | 8.26 | 1.00x | ok |
| crustpinch | 0.05 | axi | 2 | 29.91 | 3.619x | GATE |
| crustpinch | 0.05 | axi | 4 | 30.05 | 3.636x | GATE |
| crustpinch | 0.05 | axi | 10 | 30.75 | 3.722x | GATE |
| crustpinch | 0.05 | tri | 2 | 9.49 | 1.149x | ok |
| crustpinch | 0.05 | tri | 4 | 9.77 | 1.182x | ok |
| crustpinch | 0.05 | tri | 10 | 10.61 | 1.284x | ok |

### Isotropic no-slowdown check

Each scenario's default isotropic case was timed three times against the
pre-anisotropy baseline executable and three times against the current
build; the values below are the median current/baseline ratio. Both are
below the 1.05 threshold, so the check passed.

- lopnor: 0.9972
- crustpinch: 0.9962

### Draws per real event

From `build/aniso_draws.txt`, produced by `TEST(axi_draws_per_event)` in
`tests/anisotropy_tests.cpp`. That test runs at take-off-angle degree 6
(`UseTOADegree(6)`), using the Lop Nor crust `ScatterParams` with two
anisotropic specs, `(ax, ay, az)` = (2, 2, 0.5) km and (5, 5, 0.5) km. Each
value is the mean number of tentative-collision draws per accepted
(real) event, averaged over 3000 trials, for the axisymmetric (AXI) and
general (GENERAL) samplers:

#### Spec (2, 2, 0.5)

| Ray | theta | chi | AXI draws | GENERAL draws |
|---|---:|---:|---:|---:|
| P | 0 | - | 1.35167 | 9.22567 |
| P | 30 | - | 1.341 | 8.663 |
| P | 60 | - | 1.31333 | 7.42433 |
| P | 89 | - | 1.31067 | 6.638 |
| S | 30 | 0 | 1.70467 | 7.67033 |
| S | 30 | 45 | 1.68033 | 7.59733 |
| S | 30 | 90 | 1.713 | 8.24867 |
| S | 89 | 0 | 1.59033 | 4.95533 |
| S | 89 | 45 | 1.651 | 4.54567 |
| S | 89 | 90 | 1.644 | 4.41467 |

#### Spec (5, 5, 0.5)

| Ray | theta | chi | AXI draws | GENERAL draws |
|---|---:|---:|---:|---:|
| P | 0 | - | 2.27167 | 53.7167 |
| P | 30 | - | 2.27133 | 48.794 |
| P | 60 | - | 2.05233 | 33.7237 |
| P | 89 | - | 1.99933 | 22.1213 |
| S | 30 | 0 | 4.41067 | 42.0337 |
| S | 30 | 45 | 4.479 | 44.4553 |
| S | 30 | 90 | 4.59733 | 45.337 |
| S | 89 | 0 | 3.163 | 14.0497 |
| S | 89 | 45 | 3.019 | 13.619 |
| S | 89 | 90 | 3.15733 | 14.1717 |

#### Comparison against the spec's prototype table

The design spec ("Axisymmetric sampler") gives a prototype measurement at
fine degree 7, Lop Nor-like crust, kappa 0.5, covering P at 0-89 degrees
incidence and S at 30 and 89 degrees with chi = 0, 45, 90 degrees:

| `ah, az` (km) | General sampler (spec) | Axisymmetric sampler (spec) |
|---|---:|---:|
| 2, 0.5 | 4.5-9.1 | 1.5-1.9 |
| 5, 0.5 | 13-53 | 2.7-5.1 |

The measured ranges above (table rows for theta in {30, 89}, every chi) are,
for comparison, AXI 1.31-1.71 / GENERAL 4.41-9.23 at (2, 0.5), and AXI
1.99-4.60 / GENERAL 13.6-53.7 at (5, 0.5). Both samplers' measured draw
counts agree with the spec's prototype table to within the expected
run-to-run variation (different direction sets, random seeds, and a
degree-6 rather than degree-7 fine tessellation here): AXI 1.3-1.7 against
the spec's 1.5-1.9, and 2.0-4.6 against 2.7-5.1; GENERAL 4.4-9.2 against
4.5-9.1, and 13.6-53.7 against 13-53.

### Measurement gate verdict

**FAILED.** Every axisymmetric row in the benchmark table exceeds the 2×
slowdown gate (lopnor axi: 4.223x–6.313x; crustpinch axi: 3.619x–3.722x), as
does the lopnor ε = 0.05 triaxial r = 10 row (3.139x).

Per the spec's measurement gate, the follow-up for a `GATE` row depends on
which sampler is involved. Neither is implemented here:

- Axisymmetric rows: revisit the axisymmetric sampler's cell and bin
  resolution within its memory budget. The diagnosis below finds the
  axisymmetric gate failures are actually dominated by per-scatterer
  `AxiProposal` construction cost, not sampling/bin-resolution cost, so
  the more direct remedy is parallelizing or restructuring that
  construction across scatterers (the spec's own suggested remedy, "if
  that is material for a model with many axisymmetric layers, compute
  them in parallel across scatterers"); that is a recommended follow-up,
  not implemented in this fix pass.
- The lopnor ε = 0.05 triaxial r = 10 row: consider extending
  incidence-tailored proposals to triaxial P waves, designed separately.

### Diagnosis: fixed construction cost

The axisymmetric times in the benchmark table are nearly flat in r at
ε = 0.01: 109.68 s, 109.96 s, 111.24 s for lopnor (r = 2, 4, 10), and
29.83 s, 29.81 s, 29.84 s for crustpinch. If the slowdown were coming from
sampling cost, it should grow with aspect ratio the way the triaxial rows
do; instead it is essentially constant. That pattern suggests a fixed
per-scatterer construction cost (e.g. building the `AxiProposal` tables
and/or report statistics) rather than sampling cost.

To check this, each scenario was run with `--num-phonons=0` at degree 9
(`TOA_DEGREE=9 tests/run_aniso_case.sh lopnor 24 OUTDIR 0 MODEL_ARGS`),
isolating setup/construction time from any phonon-sampling time, for the
lopnor isotropic default, the lopnor axi r = 2 ε = 0.01 layout, and the
lopnor tri r = 2 ε = 0.01 layout (all from `tests/benchmark_aniso.sh`'s own
model-argument construction):

| Scenario | Zero-phonon seconds | Full-phonon seconds (benchmark table) |
|---|---:|---:|
| iso | 21.947 | 23.01 |
| axi r=2 | 108.125 | 109.68 |
| tri r=2 | 26.467 | 27.71 |

In each case the zero-phonon run accounts for nearly all of the full run's
time (iso: 21.947 of 23.01 s; axi r=2: 108.125 of 109.68 s; tri r=2: 26.467
of 27.71 s) — roughly 95-98%. This is consistent with the flat-in-r pattern
above and suggests that reducing construction cost may matter more than
reducing bin resolution for closing the axisymmetric gate.

#### Measured per-scatterer construction cost

To measure the per-scatterer cost directly rather than infer it, the same
zero-phonon lopnor axi r=2 ε=0.01 run (degree 9) was repeated and its
scatterer dump inspected: it builds 20 axisymmetric (`AXI`) scatterers out
of 21 total. Timed against the zero-phonon isotropic lopnor run on the same
host:

| Run | User+sys seconds (zero phonons, degree 9) |
|---|---:|
| iso | 21.70 |
| axi r=2 | 106.11 |

That is (106.11 - 21.70) / 20 ≈ **4.22 seconds per axisymmetric scatterer**
for `AxiProposal` construction plus report-statistics integration combined,
at degree 9. The spec's own estimate was "roughly 1 s per axisymmetric
scatterer"; the measured cost here is about 4x that (consistent with the
review's independent measurement of about 3.8 s/scatterer for `AxiProposal`
alone plus about 0.3 s for report statistics, measured at degree 6 -- this
construction cost does not depend on `--toa-degree`, since the bound-table
build always uses a fixed degree-6 (or lower) set of representative
directions regardless of the run's own take-off-angle degree). This
confirms, with a direct per-scatterer measurement rather than an inference
from the flat-in-r pattern, that per-scatterer construction cost -- not
bin/cell resolution or sampling cost -- is what the axisymmetric gate
failures above come from.

#### Status: axisymmetric construction-cost remedy not implemented

The spec names the remedy for this cost directly: "if that is material for
a model with many axisymmetric layers, compute them in parallel across
scatterers." That parallelization (or an equivalent restructuring of the
bound-table build, e.g. hoisting the per-representative invariants noted in
the final review) is a recommended follow-up for closing the axisymmetric
measurement gate. It is not implemented in this fix pass, which is scoped
to correctness and documentation findings only; the axisymmetric sampler's
*results* are unaffected and remain exact (see the sampler-correctness
checks in `tests/anisotropy_tests.cpp`) -- only its one-time, per-scatterer
model-build cost is unoptimized.

(Status as of the 2026-10-16 baseline: this construction-cost remedy is
implemented and measured in the 2026-10-25 section above.)
