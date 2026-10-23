# Axisymmetric Sampler Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make axisymmetric anisotropic runs at least as fast as general-sampler runs on the same medium. Validate them against triaxial runs with `ax ≈ ay`, and document the take-off-angle quadrature bias.

**Architecture:**
- Restructure `AxiProposal`'s bound construction: hoist per-representative work, cache per-degree geometry, build the bins in parallel, and use a trig-free cap kernel. Cache the report-statistics `GSATO` rows and integrate the report states in parallel.
- Slim the axisymmetric sampling path: locate the bin once without trig, and accept in the ray frame.
- Validate with unit tests and a replicate-run permutation test. Measure with an extended benchmark.

**Tech Stack:** C++11 (`std::thread`, `std::atomic`, `std::mutex`); POSIX sh and awk test scripts; the existing `tests/anisotropy_tests.cpp` harness.

**Spec:** `docs/superpowers/specs/2026-10-17-axisymmetric-sampler-performance-design.md`. It extends `docs/superpowers/specs/2026-10-16-anisotropic-scale-lengths-design.md`; read both.

## Global Constraints

- Isotropic runs keep today's arithmetic and random-draw order. Fixed-seed non-image outputs are byte-identical to the baseline build of `aaec540`. This is enforced by `tests/test_aniso_isotropic_regression.sh` with `BASELINE_EXE=/tmp/r3d-baseline/main`.
- Random numbers come only from `RandomUniform()` and `RandomInteger()` in `random.hpp`. Model construction draws none.
- Scatterers are read-only once constructed. Propagation writes no shared state.
- Construction-time parallelism must never change a table, a report value, or any output, whatever the thread count.
- Anisotropic outputs are identical between `--threads 1` and `--threads N` at a fixed seed (`tests/test_aniso_scenarios.sh`).
- No tabulation or interpolation of anisotropic rates or patterns. Every bound is a rigorous upper bound. A violated bound throws `std::logic_error` and is never clamped.
- Memory: each axisymmetric scatterer keeps at most about 9.2 MB. Temporary construction memory stays at about 32 MB or less per scatterer.
- Normal runs print no new stdout lines.
- Do not edit the benchmark numbers in `README.md`, `doc/MANUAL.md`, or `docs/multithreaded-simulation-handoff.md`.
- Do not change `tests/run_lopnor_thread_case.sh` or `tests/test_lopnor_thread_equivalence.sh`.
- **Reference binaries:**
  - `BASELINE_EXE=/tmp/r3d-baseline/main` (built from `aaec540`);
  - `REF_EXE=/tmp/r3d-aniso-ref/main` (built from `3db253a`, the start of this plan).
  - Both are already built. Do not rebuild them.
- **Full check**, run at the end of every code task:
  `export BASELINE_EXE=/tmp/r3d-baseline/main; make && for t in tests/test_*.sh; do "$t" || exit 1; done`.

## Review Focus

1. **Thread counts at the edges:** `--threads 1`, automatic (0), and more threads than tasks (450 bins, 17 report states).
   - A `std::logic_error` thrown inside a worker must reach the caller after every thread joins, with no deadlock and no `std::terminate`.
   - Tests: `parallel_for_covers_and_rethrows`, `axi_tables_thread_invariant` (Task 3).
2. **Bin edges and poles in the trig-free locator:** `u_rel.z` exactly 1 or 0, χ exactly 0° or 90°, sinθ ≈ 0, and directions within 1e-9 rad of a bin edge.
   - Expected: the same bin as the `acos`/`atan2` reference, except within 1e-12 rad of an edge.
   - Test: `axi_locate_matches_reference` (Task 5).
3. **Degenerate kernel inputs:** `|m_rep| = 0` (k_in = k_out with t = ẑ), `m̂·u_c = ±1`, α_c ± δ outside [0, π], and r_i = 0 (degree ≤ 6).
   - Test: `axi_kernel_matches_reference` (Tasks 2 and 4).
4. **Curved mappings and vertical incidence in the ray-frame acceptance.** `AxiRatio` and the P shortcut must reproduce the fabric-frame values in `RAE_CURVED` and at `up ∥ dir`.
   - Tests: `axi_ratio_matches_fabric`, `axi_p_urel_matches_canonical` (Task 5); the existing `axi_curved_frame` and `axi_vertical_incidence`.
5. **Several take-off-angle degrees in one process.** Unit tests build degrees 2–7 in sequence, and the geometry cache must never serve another degree's geometry.
   - Test: `axi_geometry_cache_degrees` (Task 2); the existing non-icosahedral throw test.

## Pre-flight rulings for this plan

- **R-P1:** Tasks 2, 4 and 5 put frozen copies of `3db253a`'s `ComputeBound` and `Fold` formulas in `tests/anisotropy_tests.cpp`, named `ReferenceBound3db253a` and `ReferenceLocate3db253a`. They are regression oracles for the refactor and live only in test code, with a comment naming the commit they freeze. This is deliberate duplication. If wrong, the oracles drift harmlessly once the refactor is done.

---

### Task 1: Quadrature-bias future-fix note, manual corrections, seed rule

**Files:**
- Create: `docs/toa-quadrature-bias.md`
- Modify:
  - `doc/MANUAL.md` (the Validity-limits bias paragraph near line 283, and the `--seed` entry);
  - `scatterers.cpp` (the `ComputeMFPs` TODO comment, ~lines 290–301);
  - `tests/anisotropy_tests.cpp` (the S90SH reference comment, ~lines 846–859).

**Interfaces:** none (documentation and comments only).

- [ ] **Step 1: Write `docs/toa-quadrature-bias.md`.**
  - **Mechanism.** Each take-off-angle direction is the centre of one icosahedral subdivision triangle. The code averages over directions with equal weights, but the triangles' solid angles vary from 0.928× to 1.207× the mean (standard deviation 8.6%), identically at degrees 6–9. The equal-weight average therefore carries a fixed bias that does not shrink with degree.
  - **Isotropic effect.** The offset is the same for every incidence direction and consistent between rate and pattern. Examples: Lop Nor crust S +1.11%, P −0.04%; a = 2 km P +0.62%, S +0.78%.
  - **Anisotropic effect.** The bias depends on direction and polarization. For the (2, 2, 0.5) km crust (equal → area weight):

    | Quantity | Equal weight | Area weight |
    |---|---:|---:|
    | S90SH | 397.8 km | 374.7 km |
    | S90SV | 433.9 km | 442.8 km |
    | P90 | 1178.8 km | 1203.2 km |
    | P0 | 1610.4 km | 1602.7 km |
    | S0 | 759.9 km | 754.7 km |

    The area-weighted values match the original design's independent cross-check table, and Task 6 of the 2026-10-16 plan got 374.7 km from a 200k-point Fibonacci sphere.
  - **Provenance.** The tables above come from an independent re-implementation of the sampler maths (a numpy prototype, not in the repository). Say so in the note. The 397.8 km figure is reproduced in-repo by `test_report_directional_mfps`.
  - **Affected code**, with file:function for each item:
    - `PopulateProbDists`: weight each entry by `w_k = Ω_k·N/4π`;
    - `PopulateWholeProbs`, `ComputeMFPs`, `ComputeDipoles`;
    - `GetRandomScatteredRelativePhonon`, which needs no change once the tables are weighted;
    - `IntegrateRateAndDipole`, `ComputeReportStats`, `ReportDirs`;
    - `AxiProposal`: the cell sums inherit the weights, and the bounds on R do not change;
    - the general-sampler majorant;
    - the `ShearDislocation` radiation tables in `events.cpp`.
    - `S2::TesselSphere` must expose per-direction solid angles. Today the triangles are internal to `geom_s2.cpp`.
  - **Consequences.** Every output changes, isotropic included, which breaks byte-identity with `aaec540`. The change should therefore be opt-in (e.g. `--toa-quadrature=area`, default `equal`). The forced-isotropic and seven- vs five-value equivalences still hold within either mode.
  - **Validation for a future fix:**
    - weights sum to 4π within 1e-12;
    - isotropic MFPs match an exact 2-D integral within 1e-4;
    - the anisotropic values above are reproduced;
    - the sampler χ² tests pass against the weighted targets;
    - a scenario-level report shows how much outputs move.
- [ ] **Step 2: Correct `doc/MANUAL.md`.**
  - Replace the claim that the fixed-angle average "introduces no bias" for isotropic regions. It introduces a direction-independent offset of about 1%, e.g. +1.1% for Lop Nor crust S.
  - Link the paragraph to `docs/toa-quadrature-bias.md`.
  - Keep every other sentence and number in that paragraph.
- [ ] **Step 3: Add the seed rule to the manual's `--seed` entry.** Phonon i's random stream is seeded from `seed XOR i`. Two seeds whose XOR is smaller than the phonon count therefore reuse the same set of streams (596 and 597, for example). For independent replicate runs, use seeds that differ above bit 32, e.g. `596 + k·2^32`.
- [ ] **Step 4: Link the note from the `ComputeMFPs` TODO and the S90SH test comment.** Add one line each, without changing the code.
- [ ] **Step 5: Verify.** `make && tests/test_anisotropy_unit.sh` passes, and `git diff --stat` shows only the four files above.
- [ ] **Step 6: Commit.**

```bash
git add docs/toa-quadrature-bias.md doc/MANUAL.md scatterers.cpp tests/anisotropy_tests.cpp
git commit -m "Document take-off-angle quadrature bias and the seed-stream rule"
```

---

### Task 2: Bit-identical bound-kernel restructure and geometry cache

**Files:**
- Modify: `aniso_axi.hpp`, `aniso_axi.cpp`
- Create: `tests/compare_aniso_reference.sh`
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Produces, as public members of `AxiProposal`:

```cpp
struct BoundRep { R3::XYZ MHat; Real Mr, Rho, M2Lo, M2Hi, Numer; };
    // per (representative t, r_i, k_in, k_out): m_rep = k_out*t - k_in*zhat, MHat = m_rep.Unit() (zero vector if Mr == 0),
    // Mr = |m_rep|, Rho = k_out*2*sin(r_i/2), M2Lo = max(Mr-Rho,0)^2, M2Hi = (Mr+Rho)^2, Numer = 1 + aMin^2*M2Hi
struct BinCap { R3::XYZ Uc; Real Delta, CosDelta, SinDelta; };   // one per incidence bin
static BoundRep PrepareBoundRep(const R3::XYZ & t, Real r_i, Real k_in, Real k_out, Real aMin);
static BinCap   MakeBinCap(const R3::XYZ & u_c, Real delta);
static Real     BoundRatio(const BoundRep & rep, const BinCap & cap, Real ah, Real az);
    // (1 + aMin^2*m2_hi) / (1 + Q_lo), the quantity raised to (kappa+1.5) in today's ComputeBound
```

  A cell's bound is `min(M, M·pow(max over its representatives of BoundRatio, κ+1.5))`.
- Also produces:
  - a process-wide geometry cache, keyed by take-off-angle count and verified by a fingerprint of the set;
  - `tests/compare_aniso_reference.sh`, a manual script used again by Tasks 3–5.

- [ ] **Step 1: Write the failing tests** in `tests/anisotropy_tests.cpp`.
  - Add the frozen oracle `ReferenceBound3db253a(t, r_i, k_in, k_out, u_c, delta, ah, az, aMin, M, kappa)`: today's `ComputeBound` body made free-standing, with the specification values as arguments (ruling R-P1).

```cpp
test_axi_kernel_matches_reference: for t in {TesselSphere(TESS_ICO,3) centres} ∪ {zhat};
    r_i in {0, 0.003, 0.05}; (k_in,k_out) in {(1,1),(1,1.73),(1.73,1),(1.73,1.73)};
    caps: P-style bins b in {0,22,44} (delta = 1 deg) and S-style (b,c) in {(0,0),(22,4),(44,8)}
    (delta = 1 deg + sin(theta_hi)*5 deg), u_c as in the AxiProposal constructor;
    specs (ah,az,kappa) in {(2,0.5,0.3),(2,0.5,0.5),(0.5,2,0.5),(10,0.5,0.5)}, aMin=min(ah,az), M=ah*ah*az/aMin^3:
    min(M, M*pow(BoundRatio(PrepareBoundRep(t,r_i,k_in,k_out,aMin), MakeBinCap(u_c,delta), ah, az), kappa+1.5))
      == ReferenceBound3db253a(...)   (exact ==, every combination)
test_axi_geometry_cache_degrees: build AxiProposal for LopNorCrust(2,2,0.5) at degree 5, then degree 2, then degree 5
    again (each with BuildIsoPDists): the two degree-5 builds give identical MajorantRate(P,S) and identical
    AcceptanceBound(intype, u_c(bin), conv, cell) for every P bin x conv x cell (exact ==)
```

- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.** Expect a compile failure on `BoundRep`, `PrepareBoundRep` and `MakeBinCap`.
- [ ] **Step 3: Implement the restructure in `aniso_axi.cpp`.**
  - Precompute `BoundRep` arrays for the two (k_in, k_out) pairs each bin type uses: P bins (kP,kP) and (kP,kS); S bins (kS,kP) and (kS,kS).
  - Hoist the per-cell isotropic sums `S0`, `S1` out of the bin loop; they don't depend on the bin.
  - Per bin and conversion, track the largest `BoundRatio` per cell and apply `min(M, M·pow(·, κ+1.5))` once per cell.
  - In this task `BoundRatio` evaluates α_c with `acos` and the cap with `cos`, using exactly today's expressions. The tables are then bit-identical, because `pow` and the multiply are monotone.
  - Delete `ComputeBound`.
  - **Geometry cache.** Representative centres and radii, cell counts, and mirror permutations go into a function-local static map from take-off-angle count to geometry, guarded by a `std::mutex`. Before an entry is reused, a fingerprint of `toa` is checked: the sum of all direction components in index order. On a mismatch the geometry is rebuilt without caching.
  - Keep the non-icosahedral throw where it is.
- [ ] **Step 4: Write `tests/compare_aniso_reference.sh`.**
  - It is manual and not named `test_*`.
  - `REF_EXE` is required.
  - It runs the anisotropic `lopnor` (21-value) and `crustpinch` (42-value) layouts from `tests/test_aniso_scenarios.sh`. Runs use `THREADS ∈ {1, ${CMP_THREADS:-12}}` and `TOA_DEGREE ∈ {5, 7}` (degree 7 exercises representatives with radii). Each case runs once with `R3D_EXE=$ROOT/main` and once with `R3D_EXE=$REF_EXE`, both through `tests/run_aniso_case.sh` at default phonons.
  - Each pair must pass `compare_runs … strict`.
  - It prints one PASS/DIFF line per case and exits nonzero on any DIFF.
- [ ] **Step 5: Run the tests.** `tests/test_anisotropy_unit.sh` passes, and `REF_EXE=/tmp/r3d-aniso-ref/main tests/compare_aniso_reference.sh` shows PASS on every case.
- [ ] **Step 6: Measure.** Time one zero-phonon degree-9 Lop Nor axisymmetric r=2 build with both binaries:
  - layout `0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,1.0,1.0,0.5,0.3,1000,0.8,0.01,1.4,1.4,0.7,0.5,300`;
  - `TOA_DEGREE=9`, 0 phonons, default threads.

  Record both times in the report.
- [ ] **Step 7: Run the full check, then commit.**

```bash
git add aniso_axi.hpp aniso_axi.cpp tests/anisotropy_tests.cpp tests/compare_aniso_reference.sh
git commit -m "Hoist per-representative bound work and cache AxiProposal geometry"
```

---

### Task 3: Parallel construction of axisymmetric tables and report statistics

**Files:**
- Create: `parallel.hpp`
- Modify: `aniso_axi.hpp`, `aniso_axi.cpp`, `scatterers.hpp`, `scatterers.cpp`, `model.cpp`, `Makefile` (a `parallel_hpp` variable, plus prerequisites for `aniso_axi.o`, `scatterers.o` and the test target)
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Consumes: Task 2's `BoundRep`, `BinCap` and the cache.
- Produces:

```cpp
// parallel.hpp
template <class Body> void ParallelFor(Index begin, Index end, unsigned threads, Body body);
    // Calls body(i) exactly once for every i in [begin,end), on up to 'threads' std::threads
    // (threads <= 1 or a single index: inline on the caller, no thread created). Indices are claimed
    // from a std::atomic counter. If any body throws, remaining unclaimed indices are skipped, every
    // thread is joined, and the exception from the lowest throwing index is rethrown.
// Scatterer (public static)
static void     SetBuildThreads(unsigned n);   // construction-time worker count; 0 is treated as 1
static unsigned GetBuildThreads();             // default 1
// AxiProposal
AxiProposal(const ScatterParams & par, const S2::S2Set & toa, const std::vector<ProbDist> & pdists,
            unsigned threads = 1);
```

- [ ] **Step 1: Write the failing tests.**

```cpp
test_parallel_for_covers_and_rethrows: for threads in {0,1,3,64}: ParallelFor(0,1000,threads,...) increments
    an atomic counter per index -> every index exactly once; a body that throws std::logic_error("boom") at
    i==500 (threads 3 and 64) -> ParallelFor throws std::logic_error with what()=="boom", and no index > 999 is visited
test_axi_tables_thread_invariant: specs LopNorCrust(2,2,0.5) and LopNorCrust(0.5,0.5,2) at degree 5:
    AxiProposal(sp,toas,pdists,1) vs AxiProposal(sp,toas,pdists,7): MajorantRate(P), MajorantRate(S),
    BinRate at every P bin centre and S bin centre, AcceptanceBound for every (bin, conv, cell): exact ==
test_report_stats_thread_invariant: Scatterer::SetBuildThreads(1), build Scatterer(LopNorCrust(2,2,0.5)) and
    Scatterer(LopNorCrust(2,1,0.5)); SetBuildThreads(7), build both again: GetReportMFPs() equal element-wise
    (exact ==); restore SetBuildThreads(1)
```

- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.** Expect a compile failure on `ParallelFor`, `SetBuildThreads`, and the four-argument `AxiProposal`.
- [ ] **Step 3: Implement.**
  - **Bins.** The `AxiProposal` constructor builds the 45 P and 405 S bins with `ParallelFor` over bin indices. Each task writes only its own bin's table and rate. After the join, `mPMajorant` and `mSMajorant` are reduced serially in index order.
  - **Report rows.** `ComputeReportStats` builds one row per report direction holding the four `GSATO` values and cosψ, filled with `ParallelFor` over direction chunks.
  - **Report states.** The states are integrated with `ParallelFor` over states. There are 9 axis states, plus 8 angle states for axisymmetric specs.
  - **Shared integrator.** The public `IntegrateRateAndDipole` builds rows for its `dirs` and calls the same private row-based integrator, so the integration loop is not duplicated. Each state's sums keep today's order, so the values are bit-identical.
  - **Thread count.** `Scatterer`'s constructor passes `GetBuildThreads()` to `AxiProposal` and uses it in `ComputeReportStats`. `Model::Model` calls `Scatterer::SetBuildThreads(mNumThreads ? mNumThreads : GetAutomaticWorkerLimit())` after it records the thread count and before it builds any cells. That function is the existing static in `model.cpp`.
  - Print nothing new.
- [ ] **Step 4: Run the tests.** `tests/test_anisotropy_unit.sh` passes, and `REF_EXE=/tmp/r3d-aniso-ref/main tests/compare_aniso_reference.sh` passes every case.
- [ ] **Step 5: Measure** as in Task 2 Step 6. Add one more run with the triaxial layout `0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,1.0,0.707107,0.5,0.3,1000,0.8,0.01,1.4,0.989949,0.7,0.5,300`. Record the times.
- [ ] **Step 6: Run the full check, then commit.**

```bash
git add parallel.hpp aniso_axi.hpp aniso_axi.cpp scatterers.hpp scatterers.cpp model.cpp Makefile tests/anisotropy_tests.cpp
git commit -m "Build axisymmetric tables and report statistics in parallel"
```

---

### Task 4: Trig-free cap kernel

**Files:**
- Modify: `aniso_axi.cpp` (the `BoundRatio` body)
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Consumes: `BoundRatio`, `BinCap` and `MakeBinCap` (Task 2). The signatures are unchanged.

- [ ] **Step 1: Change the test.** In `test_axi_kernel_matches_reference`, replace exact `==` with relative agreement `|fast − ref| ≤ 1e-12·ref`, and add the edge cases:

```cpp
caps with u_c == rep.MHat (c = +1) and u_c == -rep.MHat (c = -1) for t = (0.6,0,0.8), k_in=1, k_out=1.73;
caps whose alpha_c is within delta of 0 and within delta of pi (u_c rotated 0.5*delta from +-MHat);
Mr == 0 (t = zhat, k_in = k_out = 1.73) with every cap above
```

- [ ] **Step 2: Run the test.** It still passes, because the `acos` kernel is exact.
- [ ] **Step 3: Rewrite `BoundRatio`.**
  - With `c = clamp(MHat·Uc, −1, 1)` and `s = √(1−c²)`:
    - `cos α_lo = (c ≥ CosDelta) ? 1 : c·CosDelta + s·SinDelta`;
    - `cos α_hi = (c ≤ −CosDelta) ? −1 : c·CosDelta − s·SinDelta`.
  - When `Mr == 0`, use `c = 1`.
  - Everything after s_lo and s_hi is unchanged.
  - No `acos`, `cos` or `sin` calls remain inside `BoundRatio`.
- [ ] **Step 4: Run the tests.**
  - `tests/test_anisotropy_unit.sh` passes, including `axi_bounds_valid`, `axi_low_degree` and `axi_extreme_ratio`.
  - Run `REF_EXE=/tmp/r3d-aniso-ref/main tests/compare_aniso_reference.sh` and record the result. Identity is expected but not required: the tables now differ from the reference at rounding level.
  - Measure as in Task 3 Step 5.
- [ ] **Step 5: Run the full check, then commit.**

```bash
git add aniso_axi.cpp tests/anisotropy_tests.cpp
git commit -m "Evaluate the AxiProposal cap bound without trigonometric calls"
```

---

### Task 4b: Parallelize the remaining serial AxiProposal setup (added 2026-10-23, ruling R3)

Profiling at `961ca45` put the degree-9 Lop Nor zero-phonon build at axi 25.8 s against tri 22.2 s. The 3.6 s gap breaks down as:
- ~2.2 s: the one-time geometry build (`BuildAxiGeometry`) runs serially. It computes 81,920 × 64 representative angles with `acos` and does a brute-force 1280² mirror match.
- ~0.4 s: the `PrepareBoundRep` arrays are built serially, per scatterer.
- ~0.7 s: the residual of the parallel bin build. The P table's 45 tasks and the 8- and 9-state report calls cannot fill 16 threads.
- ~0.3 s: thread spawn and join. That is four `ParallelFor` calls per scatterer, about 1,100 short-lived threads per run.

**Files:**
- Modify: `aniso_axi.hpp`, `aniso_axi.cpp`, `scatterers.cpp`
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Consumes: `ParallelFor`, `Scatterer::GetBuildThreads()` (Task 3), and the geometry cache and `BoundRep` arrays (Task 2).
- Produces:
  - `BuildAxiGeometry(toa, threads)`, which takes the build thread count. The cached geometry is identical for any thread count.
  - One `ParallelFor` over all 450 bins (P bins first, then S, in today's index order) instead of two calls.
  - One `ParallelFor` over all report states (9 axis states, plus 8 angle states for axisymmetric specs) instead of two.

- [ ] **Step 1: Write the failing test.**

```cpp
test_axi_geometry_thread_invariant: for degrees 5 and 7, build the AxiProposal geometry with threads 1 and 7
    (via the public AxiProposal constructor on fresh, distinct take-off-angle sets so each build misses the
    cache — e.g. the rotated set from test_axi_geometry_cache_mismatch_rotation plus the standard set, in an
    order that guarantees a miss — or via a test hook that bypasses the cache): MajorantRate and every
    AcceptanceBound for P bins x conv x cell are exact == between the two thread counts
```

  If no cache-bypass is possible without a hook, add a narrowly named test-only free function in `aniso_axi.cpp`'s public header section, `AxiGeometryChecksum(const S2::S2Set &, unsigned threads)`. It returns a checksum over every representative radius and mirror-permutation entry, and the test compares checksums across thread counts.

- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.** Expect a compile failure (new signature or hook).
- [ ] **Step 3: Implement.**
  - **Geometry.** Parallelize the representative-radius loop and the mirror-match loop with `ParallelFor`. Each index writes only its own element. Keep the "no match" throws, which `ParallelFor` rethrows from the lowest index.
  - Parallelize the `PrepareBoundRep` precompute loops.
  - Merge the P and S bin builds into one `ParallelFor` and the report states into one `ParallelFor`. Keep the serial majorant reduction and every per-state summation order.
  - Fix the stale `BoundRatio` doc comments in `aniso_axi.hpp` (~lines 175 and 224–226), which still mention `acos` and `cos` (deferred minor from Task 4).
- [ ] **Step 4: Run the tests.** `tests/test_anisotropy_unit.sh` passes, and `REF_EXE=/tmp/r3d-aniso-ref/main tests/compare_aniso_reference.sh` passes 8/8.
- [ ] **Step 5: Measure.** Time zero-phonon degree-9 Lop Nor builds, axi and tri, as in Task 3 Step 5, plus the same pair with `--threads=1`. Record all four.
- [ ] **Step 6: Run the full check, then commit.**

```bash
git add aniso_axi.hpp aniso_axi.cpp scatterers.cpp tests/anisotropy_tests.cpp
git commit -m "Parallelize AxiProposal geometry and representative setup"
```

---

### Task 5: Faster axisymmetric sampling path

**Files:**
- Modify:
  - `scatparams.hpp`, `scatparams.cpp` (`AxiRatio`);
  - `aniso_axi.hpp`, `aniso_axi.cpp` (`BinRef`, `Locate`, overloads);
  - `scatterers.cpp` (`AttemptAxisymmetricScatter`).
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Produces:

```cpp
// ScatterParams
Real AxiRatio(Real m2, Real mu2) const;
    // AnisoRatio for an axisymmetric spec from |m|^2 and (m.up)^2:
    // M * pow((1 + a^2*m2) / (1 + ax^2*max(m2-mu2,0) + az^2*mu2), kappa+1.5)
// AxiProposal (public)
struct BinRef { int ThetaBin; int ChiBin; bool Mx; bool My; };   // replaces the private FoldResult
BinRef Locate(raytype intype, const R3::XYZ & u_rel) const;     // no acos/atan2/sin
Real   BinRate(raytype intype, const BinRef & bin) const;
Draw   DrawProposal(raytype intype, const BinRef & bin) const;
// BinRate(intype,u_rel), DrawProposal(intype,u_rel), AcceptanceBound(...) remain, implemented via Locate()
```

- [ ] **Step 1: Write the failing tests.**
  - Add the frozen oracle `ReferenceLocate3db253a(intype, u)`: today's `Fold` made free-standing, returning `(thetaBin, chiBin, mx, my)` (ruling R-P1).

```cpp
test_axi_locate_matches_reference: P and S; u over 200000 Fibonacci-sphere directions, plus for every theta edge
    (2b deg) and chi edge (10c deg) the directions at the edge +-1e-9 rad (at chi 45 deg / theta 45 deg for the
    other angle), plus u = (0,0,+-1), (+-1,0,0), (0,+-1,0): Locate == ReferenceLocate3db253a, except that for
    directions within 1e-12 rad of an edge either adjacent bin is accepted
test_axi_ratio_matches_fabric: for mappings ENU_ORTHO and RAE_CURVED (loc 950 km along x on the surface),
    1000 deterministic (dir, toa direction t, conv) draws with spec LopNorCrust(2,2,0.5): build AA from
    (dir, pol), m = k_out*AA.Express(t) - k_in*dir; AxiRatio(m.MagSquared(), sq(m.Dot(up)))
    == AnisoRatio(ToFabric(loc,m)) within 1e-12 relative; restore ENU_ORTHO
test_axi_p_urel_matches_canonical: same mappings and locs, 1000 dirs plus dir == +-up: with cz = up.d,
    (-sqrt(max(0,1-cz^2)), 0, cz) equals (up.S1, up.S2, up.S3) of OrthoAxes(dir, CanonicalPPolarization(dir,up))
    within 1e-12 componentwise
```

- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.** Expect a compile failure on `AxiRatio`, `Locate` and `BinRef`.
- [ ] **Step 3: Implement.**
  - **`Locate`.** It folds as `Fold` did. The θ bin `b` satisfies `cos(2(b+1)°) < u_z ≤ cos(2b°)`, found by binary search over 46 precomputed edge cosines and clamped to [0, 44].
  - For S, after folding, `w = −u_x ≥ 0` and `v = u_y ≥ 0`. χ bin `c` is the largest c in [0, 8] with `v·cos(10c°) ≥ w·sin(10c°)`. When `w² + v² ≤ 1e-24`, c = 0.
  - **`AttemptAxisymmetricScatter`:**
    1. `up = ECS.GetUp(loc)`.
    2. **`u_rel`.** For P, `u_rel = (−√max(0, 1−cz²), 0, cz)` with `cz = clamp(up·d, −1, 1)`. For S, `u_rel` comes from `OrthoAxes(dir, pol)` as today.
    3. `bin = Locate(intype, u_rel)`, once.
    4. The thinning test uses `BinRate(intype, bin)`.
    5. **Each proposal:** `t = (*pTOA)[Toa].XYZ()`, `mu = k_out·t·u_rel − k_in·u_rel.z`, `m2 = k_out² + k_in² − 2·k_out·k_in·t.z`, and `p = AxiRatio(m2, mu²)·CellSum/Weight`.
    6. **On acceptance:** `InPol` is `CanonicalPPolarization(dir, up)` for P, or `pol` for S. Then `BuildRelativePhonon` as today.
  - **Kept as today:**
    - the `RandomUniform()` call sequence;
    - both `std::logic_error` checks (non-positive Weight, and `p > 1 + kBoundMargin`) with their messages;
    - the override loop;
    - `Draws` counting.
  - **Non-real attempts** carry the incoming `pol` in `InPol` (`Propagate` uses `InPol` only for real events). Check that no test reads `InPol` on a non-real attempt.
- [ ] **Step 4: Run the tests.**
  - `tests/test_anisotropy_unit.sh` passes, including `axi_matches_general`, `axi_bounds_valid`, `axi_draws_per_event` (thresholds unchanged), `axi_vertical_incidence`, `axi_curved_frame`, `axi_low_degree`, `axi_extreme_ratio` and `isolated_axi_overridemfp`.
  - `tests/test_aniso_scenarios.sh` passes.
  - **Measure** degree-9 Lop Nor total time at 1M phonons with ε=0.05 for the axisymmetric r=10 layout `0.8,0.05,0.5,0.5,0.5,0.2,50,0.8,0.05,5.0,5.0,0.5,0.3,1000,0.8,0.05,7.0,7.0,0.7,0.5,300`, with this build and with `REF_EXE`. Record both times.
- [ ] **Step 5: Run the full check, then commit.**

```bash
git add scatparams.hpp scatparams.cpp aniso_axi.hpp aniso_axi.cpp scatterers.cpp tests/anisotropy_tests.cpp
git commit -m "Locate incidence bins once and accept axisymmetric proposals in the ray frame"
```

---

### Task 6: Near-axisymmetric validation

**Files:**
- Modify:
  - `tests/aniso_lib.sh`: move `region5`, `region7`, `mode_mults`, `lopnor_model`, `crustpinch_model` and `scenario_model` out of `tests/benchmark_aniso.sh`, and add mode `gen`;
  - `tests/benchmark_aniso.sh`: source them instead.
- Create: `tests/compare_axi_general.sh`
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Consumes: `MeanDistanceToRealEvent`, the 4a in-states, the two-sample χ² helper used by `test_axi_matches_general`, `GetReportMFPs`, and `GetMajorantMeanFreePath`.
- Produces:
  - layout mode `gen` in `scenario_model`: the same as `axi`, except `ay = ax·(1+1e-5)`;
  - multiplied values printed with `%.9f` in every mode.

- [ ] **Step 1: Write the unit tests.** All use degree 5 unless noted. "Near" means `ay = ax·(1+δ)`.

```cpp
test_near_axi_rate_matches: specs (2,2,0.5),(5,5,0.5),(0.5,0.5,2) x kappa {0.3,0.5}; 4a in-states; for
    delta in {1e-6,1e-4,1e-2}: |IntegratedRate(near) - IntegratedRate(ax=ay)| <= 5*delta*IntegratedRate(ax=ay);
    and where |dR(1e-6)| > 1e-13*R: dR(1e-4)/dR(1e-6) in [50,200]
test_near_axi_sampler_matches: Scatterer(LopNorCrust(2,2,0.5)) (AUTO -> axisymmetric) vs
    Scatterer(LopNorCrust(2,2*(1+1e-5),0.5)) (AUTO -> general): for each 4a in-state, MeanDistanceToRealEvent
    n=20000 each: |difference| < 4*sqrt(sigma_a^2+sigma_g^2); for P Dir(60,30) and S along E pol U,
    200000 accepted (Conv, degree-2 cell) draws each: two-sample chi2/dof < 1.5 over bins with expected >= 20
test_near_axi_report_matches: GetReportMFPs, AXI (8) for ax=ay vs ANI (9) for near, within 1e-4 relative:
    P0~P_U; P90~P_E and P_N; S0~S_U/E and S_U/N; S90SV~S_E/U and S_N/U; S90SH~S_E/N and S_N/E
test_near_axi_cost (degree 7): specs (1,1,0.5),(2,2,0.5),(5,5,0.5); deterministic: axi
    GetMajorantMeanFreePath(rt) >= near general GetMajorantMeanFreePath(rt) * (1-1e-5) for P and S;
    measured over 5000 real events per state (P at theta 0,30,60,89; S at theta 30,89 with chi 0,45,90):
    tentative collisions (GetRandomPathLength calls) and draws per real event for both samplers, written to
    build/aniso_near_axi_cost.txt; axi mean tentative <= general mean tentative + 4*sigma_diff in every state
```

- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.** All pass. A failure here is a finding about the theory or the samplers: report it, do not tune it.
- [ ] **Step 3: Move the layout builders into `tests/aniso_lib.sh`** and source them from `tests/benchmark_aniso.sh`. For every (scenario, ε, mode, r) the benchmark used, the layouts must parse to the same numbers. Check by printing each layout before and after with awk-normalized numbers and diffing.
- [ ] **Step 4: Write `tests/compare_axi_general.sh`.** It is manual and not named `test_*`.
  - **Environment, with defaults:**
    - `CMP_DEGREE=7`, `CMP_LOPNOR_PHONONS=250000`, `CMP_CRUSTPINCH_PHONONS=50000`;
    - `CMP_REPLICATES=8`, `CMP_THREADS=auto`, `CMP_PERMUTATIONS=1000`;
    - `CMP_CONFIGS` (default all eight: {lopnor, crustpinch} × ε {do-script, 0.05} × r {4, 10}).
  - **Groups:** A is mode `axi`, G is mode `gen`, and A′ is mode `axi` with different seeds.
  - **Seeds:** replicate k of A uses `596 + k·2^32` for k = 1..K; G uses k = K+1..2K; A′ uses k = 2K+1..3K.
  - **Statistic.** From each run's `seis_traces_asc.dat`, take every station's energy components summed into 10-s windows (5 bins of 2 s). Keep windows where every replicate of both compared groups has positive energy. For each window, compute `z = (Ā−Ḡ)/√(s_A²/K + s_G²/K)` and let `T = mean z²`.
  - **p-value:** the fraction of `CMP_PERMUTATIONS` random relabellings of the 2K runs with T ≥ T_obs. Use a fixed-seed LCG in awk. Compute it for A vs G and for A vs A′.
  - **Timing:** record each run's wall-clock seconds.
  - **Output:** a Markdown table per config with columns scenario, ε, r, windows, p(A,G), p(A,A′), median s (A), median s (G), and the G/A time ratio.
  - **Exit nonzero** if any p < 0.01 or any run fails.
- [ ] **Step 5: Smoke-run the script.** Use `CMP_DEGREE=4 CMP_LOPNOR_PHONONS=20000 CMP_CRUSTPINCH_PHONONS=5000 CMP_REPLICATES=4 CMP_PERMUTATIONS=200 CMP_CONFIGS="lopnor:0.05:4"` (document the `CMP_CONFIGS` syntax in the header), and put the full output in the report. Do not run it at full size.
- [ ] **Step 6: Run the full check, then commit.**

```bash
git add tests/anisotropy_tests.cpp tests/aniso_lib.sh tests/benchmark_aniso.sh tests/compare_axi_general.sh
git commit -m "Validate the axisymmetric sampler against near-axisymmetric general runs"
```

---

### Task 7: Benchmark extension, measurements, and performance documentation

**Files:**
- Modify: `tests/benchmark_aniso.sh`, `docs/anisotropy-performance.md`, `doc/MANUAL.md` (the construction-cost paragraph in Validity limits)

**Interfaces:**
- Consumes: `scenario_model` with modes `iso`, `axi`, `tri` and `gen` (Task 6), `tests/compare_axi_general.sh` (Task 6), and `build/aniso_near_axi_cost.txt` (Task 6).

- [ ] **Step 1: Extend `tests/benchmark_aniso.sh`.**
  - **`gen` rows:** for every axi row, a matching `gen` row on the same medium (general sampler).
  - **Build column:** a zero-phonon build time per (scenario, ε, mode, r).
  - **Repeats:** `BENCH_REPEATS` (default 3). Report the median of the totals.
  - **Columns:** scenario, ε, mode, r, build s, total s, slowdown vs iso, gate vs iso (`ok` ≤ 2×, else `GATE`), and, on axi rows, the same-medium ratio axi/gen with its gate (`ok` ≤ 1.05, else `SLOW`).
  - Keep the isotropic current/baseline check.
  - **Exit nonzero** on a failed run, an isotropic ratio > 1.05, or any `SLOW` row. A `GATE` row is reported but does not affect the exit status.
  - Smoke-run it at small sizes, e.g. `BENCH_DEGREE=4 BENCH_LOPNOR_PHONONS=20000 BENCH_CRUSTPINCH_PHONONS=5000 BENCH_REPEATS=1`.
- [ ] **Step 2: Controller runs the measurements.** Run at default sizes:
  - `BASELINE_EXE=/tmp/r3d-baseline/main tests/benchmark_aniso.sh`;
  - `tests/compare_axi_general.sh`.

  Capture the output to the plan workspace. If any row is `SLOW`, the controller rules on a fallback, such as AUTO choosing GENERAL for axisymmetric specs below some M, before continuing.
- [ ] **Step 3: Update `docs/anisotropy-performance.md`** from the captured output only.
  - Add a new dated section with: host, compiler and workload; the benchmark table; build times; the same-medium verdict; the iso-gate verdict; the near-axisymmetric validation table (p-values and timing); and the tentative-collision and draw counts from `build/aniso_near_axi_cost.txt`.
  - Keep the earlier section, retitled as the `3db253a` baseline.
  - Explain that the axi and tri rows describe different media: M is ah²·az/a_min³ against ah·√r·a·az/a_min³.
- [ ] **Step 4: Update the manual's construction-cost paragraph** with the measured per-scatterer build cost and the remark that it is now parallel. Keep the quadrature paragraph from Task 1.
- [ ] **Step 5: Run the full check, then commit.**

```bash
git add tests/benchmark_aniso.sh docs/anisotropy-performance.md doc/MANUAL.md
git commit -m "Measure axisymmetric performance against the same-medium general sampler"
```
