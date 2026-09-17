# Parallel Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

## Status

This plan was completed on `parallel-performance`, fast-forwarded into local
`master` at `8b9713d`, and verified with the merged-result test suite. The
feature branch and temporary worktree were subsequently removed at the user's
request. This file is retained as an historical design and execution record;
current repository state and future follow-up are maintained in `PLANS.md` and
`docs/DEVELOPMENT.md`.

**Goal:** Remove shared reporting and seismometer accumulation from the normal multi-worker propagation path and reduce worker scheduling atomics without changing the CLI or seeded per-phonon streams.

**Architecture:** Workers receive explicit `SimulationReportContext` objects containing private counters, diagnostic flags, and flattened per-seismometer time bins. `DataReporter` keeps immutable geometry/configuration, writes into the worker context during propagation, and merges contexts after all workers join. Phonon indices are allocated in fixed-size chunks from one atomic counter so random stream indices remain stable while atomic traffic falls by roughly the chunk size.

**Tech Stack:** C++11, GCC/G++, `std::thread`, `std::atomic`, `std::vector`, GNU Make, shell regression tests, and assertion-based native tests.

**Spec:** `docs/superpowers/specs/2026-09-17-parallel-report-reduction-design.md`

## Global Constraints

- Keep the C++11 and `-pthread` build unless a deliberate migration is approved.
- Preserve serial defaults and backward-compatible command-line aliases.
- Use explicit seeded random streams for stochastic tests and comparisons.
- Keep model units and local vertical/radial coordinate conventions explicit.
- Run `make -j2` and `make test` after simulation or build changes.
- Integer counters and invalidity flags must agree between seeded serial and parallel runs.
- Seismometer floating-point bins may differ by normal parallel reduction order.
- Do not add MPI or worker-time mutation to shared probability tables.

## File Map

- Create `simulation_report.hpp` for the worker-local bin and counter data model.
- Modify `dataout.hpp` and `dataout.cpp` for worker-local collection, optional locked text output, and post-join reduction.
- Modify `phonons.hpp`, `phonons.cpp`, `events.hpp`, and `events.cpp` to pass an explicit report context through the worker propagation path while retaining compatibility overloads.
- Modify `model.hpp` and `model.cpp` for context ownership, chunked phonon scheduling, and post-join reduction.
- Create `tests/test_report_reduction.cpp` for pure context/bin reduction behavior.
- Create `tests/test_parallel_reproducibility.sh` for a small seeded process-level comparison.
- Modify `Makefile` to build/run both tests.
- Modify `docs/DEVELOPMENT.md`, `docs/MANUAL.md` only if public semantics need clarification, and `PLANS.md` with the completed architecture and benchmark findings.

---

### Task 1: Add the worker-local report data model

**Files:**
- Create: `simulation_report.hpp`
- Create: `tests/test_report_reduction.cpp`
- Modify: `Makefile:17-23, test prerequisites and recipes`

**Interfaces:**
- `struct SimulationReportBin` contains `Real mEnergyAxes[3]`, `Real mEnergyByType[RAY_NUMBASICTYPES]`, and `unsigned mCountByType[RAY_NUMBASICTYPES]`, with a zeroing constructor and `Add(const SimulationReportBin&)` method.
- `class SimulationReportContext` has constructors `SimulationReportContext()` and `SimulationReportContext(std::size_t seismometer_count, std::size_t bin_count)`, `Reset(std::size_t, std::size_t)`, `SimulationReportBin* BinsFor(std::size_t)`, `const SimulationReportBin* BinsFor(std::size_t) const`, `Merge(const SimulationReportContext&)`, `RecordLost()`, `RecordTimeout()`, `RecordInvalid(unsigned reason)`, and read-only count/size accessors.
- `Merge` adds every bin, sums counters, and ORs invalidity flags. It asserts matching seismometer and bin dimensions.

- [x] **Step 1: Write the failing reduction test.**

Add tests that construct two contexts with two seismometers and three bins, add distinct axis/type/count values to both contexts, record lost/timeout/invalid events, merge the second into the first, and assert every scalar and bin field has the expected sum.

- [x] **Step 2: Run the focused test to verify it fails.**

Run: `g++ tests/test_report_reduction.cpp -std=c++11 -pthread -O3 -g -Wall -I. -o /tmp/r3d-test-report-reduction && /tmp/r3d-test-report-reduction`

Expected: compilation fails because `simulation_report.hpp` and its context API do not yet exist.

- [x] **Step 3: Implement the minimal context and bin types.**

Use a flat `std::vector<SimulationReportBin>` of size `seismometer_count * bin_count`; `BinsFor(i)` returns the start of the `i`th seismometer range. Keep all fields private in the context and expose only the named methods.

- [x] **Step 4: Run the focused test to verify it passes.**

Run the command from Step 2. Expected: the executable prints its pass message and exits zero.

- [x] **Step 5: Add the test to the Makefile and commit.**

Build `tests/test_report_reduction` from its source and `simulation_report.hpp`, add it to the `test` prerequisite and recipe, then run `make test`. Commit with `feat: add worker-local report reduction state`.

### Task 2: Make seismometer accumulation worker-local

**Files:**
- Modify: `dataout.hpp:70-222, 228-447`
- Modify: `dataout.cpp:42-216, 523-625, 631-710`
- Modify: `tests/test_report_reduction.cpp` if a focused public helper assertion is needed

**Interfaces:**
- `Seismometer::CatchPhonon(const Phonon&, SimulationReportBin*) const` performs the existing geometry/time/bin calculation and writes to the supplied bins without touching `mTimeBins`.
- The existing `Seismometer::CatchPhonon(const Phonon&)` remains as a compatibility wrapper that calls the new overload with `mTimeBins`.
- `DataReporter::CreateWorkerContext() const` returns a context sized for the current seismometer count and static bin count.
- `DataReporter::MergeWorkerContext(const SimulationReportContext&)` adds worker bins to canonical `mTimeBins`, sums aggregate counters, and ORs invalidity diagnostics after workers have joined.
- Each report method gains a context overload, for example `ReportPhononCollected(SimulationReportContext&, const Phonon&)`; existing no-context overloads remain for compatibility.

- [x] **Step 1: Add the worker-context API declarations and compile-failing call sites.**

Include `simulation_report.hpp`, expose the bin record type needed by `Seismometer`, declare the new `CatchPhonon`, context factory/merge methods, and context report overloads. Do not alter the old overloads yet.

- [x] **Step 2: Run `make -j2` and record the expected compile failures.**

Run: `make -j2`

Expected: declarations compile only after matching definitions are supplied; any missing-definition or overload errors identify the exact implementation sites for the next step.

- [x] **Step 3: Move the existing collection arithmetic into the supplied-bin overload.**

Make the geometry/configuration reads const. Preserve the current time-window, radius, axis-energy, ray-type-energy, and count formulas exactly. The compatibility wrapper continues to update canonical bins for legacy callers.

- [x] **Step 4: Implement context creation and reduction.**

Construct a context using `mSeismometers.size()` and `Seismometer::NumberOfBins()`. During merge, iterate seismometers and bins, call `SimulationReportBin::Add`, then merge counters and diagnostics. Do not take `mReportMutex`; the caller invokes this only after all workers have joined.

- [x] **Step 5: Implement context report overloads without common-path locks.**

For context callbacks, update only context counters/diagnostics and private bins. Wrap textual `output_phonon_dataline` calls in `mReportMutex` only when the corresponding report flag is enabled. `ReportPhononCollected` may lock only around its optional text write; the seismometer scan must use `context.BinsFor(i)` outside the lock.

- [x] **Step 6: Run the existing native suite and commit.**

Run: `make -j2 && make test`

Expected: build succeeds, the existing parallel/anisotropy/script tests pass, and the reduction test remains green. Commit with `perf: accumulate reports per worker`.

### Task 3: Pass contexts through propagation and reduce scheduler atomics

**Files:**
- Modify: `events.hpp:80-90`, `events.cpp:111-128`
- Modify: `phonons.hpp:30-35, 315-325`, `phonons.cpp:514-700`
- Modify: `model.hpp` simulation declarations and worker state
- Modify: `model.cpp:630-677`

**Interfaces:**
- Add `struct SimulationReportContext;` forward declarations where needed.
- Add `ShearDislocation::GenerateEventPhonon(RandomEngine&, SimulationReportContext&)` and keep the existing overload as a compatibility wrapper.
- Add `Phonon::Propagate(RandomEngine&, SimulationReportContext&)` and keep `Propagate(RandomEngine&)` and `Propagate()` as compatibility overloads.
- Change the worker entry point to accept `std::atomic<long>& next_phonon` and `SimulationReportContext& worker_context`.

- [x] **Step 1: Add a seeded process-level regression command before wiring the worker.**

Create `tests/test_parallel_reproducibility.sh` that runs `./main` twice with `--reports=ALL_OFF`, `--num-phonons=512`, `--toa-degree=3`, `--source=EXPL`, `--source-loc=425.54,-169.53,-1.02`, `--frequency=2`, `--timetolive=10`, `--binsize=1`, `--grid-compiled=1`, `--range=1200`, `--flatten`, the three existing Lop Nor `--model-args` groups, the same fixed seed, and `--workers=1`/`--workers=2`. Capture each command's direct exit status and compare only the `Loss surfaces`, `Timeout`, and `Invalidity` summary lines.

- [x] **Step 2: Run the regression before implementation.**

Run: `make -j2` followed by `tests/test_parallel_reproducibility.sh`

Expected: the script either cannot build/run until its Makefile prerequisite is added or exposes the current implementation's stable-counter behavior; it must not treat a nonzero `main` exit as a comparison pass.

- [x] **Step 3: Add context-aware event generation and propagation.**

Have the worker path call `GenerateEventPhonon(rng, context)` and `P.Propagate(rng, context)`. Inside propagation, dispatch each event to the context overload. Keep the old path calling the old report overloads so existing callers remain source-compatible.

- [x] **Step 4: Replace one-phonon scheduling with chunked index allocation.**

Use `const long kPhononWorkChunk = 256`. Each worker calls `fetch_add(kPhononWorkChunk)` once per chunk, processes `[begin, min(begin + kPhononWorkChunk, mNumPhonons))`, and derives every RNG stream from the unchanged absolute phonon index. Allocate one context per worker before launching threads, join all workers, then call `dataout.MergeWorkerContext` for each context before post-simulation output.

- [x] **Step 5: Add the regression to `make test` and run it.**

Make `main` and `tests/test_parallel_reproducibility.sh` prerequisites, mark the script executable, and run `make test`. Expected: both seeded runs exit zero and the three stable summary lines match.

- [x] **Step 6: Commit the worker wiring.**

Run: `make -j2 && make test`

Commit with `perf: remove shared reporting from worker hot path`.

### Task 4: Verify numerical behavior and document the concurrency boundary

**Files:**
- Modify: `docs/DEVELOPMENT.md` execution flow, parallel execution, testing, and limitations sections
- Modify: `docs/MANUAL.md` reporting/parallel notes only if the user-visible report-order or floating-point caveat is absent
- Modify: `PLANS.md` current baseline and prioritized follow-up sections

- [x] **Step 1: Add documentation for worker-local reduction.**

Document that seeded per-phonon random streams remain stable, integer summaries are expected to match, report ordering is not guaranteed, and multi-worker floating-point seismometer bins may differ in the last bits because worker-local bins are reduced after propagation.

- [x] **Step 2: Run the full native verification.**

Run: `make -j2` and `make test`. Record compiler warnings separately from test failures; both commands must exit zero.

- [x] **Step 3: Run the plotting regression when available.**

Run: `make test-plotting`

Expected: the plotting regression remains independent of the threading changes; if optional Octave/gnuplot tools are unavailable, record that limitation rather than changing the simulation implementation.

- [x] **Step 4: Run the matched performance benchmark.**

From the branch root, run the same 10M-phonon command used by `do-lopnor-big.sh` with `--workers=1`, then with `--workers=2`, then with `--workers=4`, all using seed `0x5eedc0de12345678`. Measure only the `main` command interval from each logfile, exclude figure generation, and record machine/compiler/build revision, worker count, seed, phonons, summary counters, elapsed seconds, and speedup relative to one worker.

- [x] **Step 5: Review the diff and commit the documentation.**

Run: `git diff --check; git status --short --branch; git diff --stat`

Commit with `docs: document parallel report reduction and benchmark results`.

### Task 5: Final verification and handoff

**Files:**
- No new files; verify the branch and committed history.

- [x] **Step 1: Run the complete verification set on the final tree.**

Run: `make -j2`, `make test`, and `make test-plotting` when optional plotting dependencies are installed.

- [x] **Step 2: Confirm the worktree and branch state.**

Run: `git status --short --branch` and `git log --oneline --decorate -8`. The
worktree must be clean on merged `master`; generated benchmark data must remain
ignored.

- [x] **Step 3: Report the implementation and measured result.**

Include the branch name, commit identifiers, tests run, benchmark configuration,
serial/parallel elapsed times, speedups, any numerical caveat, and the
historical isolated-worktree path when relevant.
