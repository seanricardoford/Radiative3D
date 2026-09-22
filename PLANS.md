# Radiative3D Continuation and Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Preserve the merged multithreaded and anisotropic-scattering baseline while making future numerical, performance, and integration work reproducible and reviewable.

**Architecture:** `main.cpp` parses options into `ModelParams`; `Model` constructs the compiled-in Earth model and owns the simulation worker boundary; `Phonon`, `Scatterer`, `ScatterParams`, and `RandomEngine` implement stochastic propagation; `DataReporter` owns canonical output plus worker-local report reduction. New work should extend these boundaries without introducing hidden global state or worker-time mutation of shared probability tables.

**Tech Stack:** C++11, GCC/G++, `std::thread`, `std::atomic`, `std::mutex`, `std::mt19937_64`, GNU Make, shell run scripts, and native assertion-based C++ regression tests.

**Spec:** `docs/DEVELOPMENT.md`

## Global Constraints

- Keep the C++11 and `-pthread` build unless a deliberate migration is approved.
- Preserve serial defaults and backward-compatible command-line aliases.
- Use explicit seeded random streams for stochastic tests and comparisons.
- Keep model units and local vertical/radial coordinate conventions explicit.
- Run `make -j2` and `make test` after simulation or build changes.
- Update `docs/MANUAL.md`, `docs/DEVELOPMENT.md`, and this file when public behavior or project state changes.

---

## Current baseline

The capability and performance baseline is on local `master` at commit
`8b9713d`. The `parallel-performance` branch was fast-forwarded into `master`
and deleted after merged-result verification. Local `master` is currently ahead
of `origin/master`; it has not been pushed unless a later session changes that
state. Always re-check Git state because branch pointers and remote tracking
information change over time.

## Completed capabilities

### Multi-worker simulation

- `--workers=N`, `-W N`, and `--threads=N` select shared-memory worker count.
- The default is `std::thread::hardware_concurrency()`, falling back to one.
- `--seed=VALUE` sets the base 64-bit seed; decimal and `0x` hexadecimal forms
  are accepted by the existing command-line parser.
- Each phonon gets a stream derived from the base seed and stable phonon index
  through `RandomEngine::SeedForStream()`.
- `Model::RunSimulation()` caps workers at the requested phonon count, uses an
  atomic next-index counter to assign 256-phonon chunks, and joins all workers
  before post-simulation output.
- Stochastic APIs have explicit `RandomEngine&` overloads; no-argument overloads
  remain as thread-local compatibility wrappers.
- `SimulationReportContext` stores worker-local counters, diagnostics, and
  flattened seismometer bins. `DataReporter` merges contexts after workers
  join; only enabled text reports use the mutex in the worker hot path. Lazy
  probability integration is completed before workers start.

Primary files: `model.cpp`, `model.hpp`, `probability.cpp`, `probability.hpp`,
`phonons.cpp`, `phonons.hpp`, `sources.cpp`, `sources.hpp`, `dataout.cpp`,
`dataout.hpp`, and `Makefile`.

### Anisotropic scattering

- `--scatter-horizontal=L_H` and `--scatter-vertical=L_V` configure a global
  ellipsoidal correlation spectrum.
- Both options are required together for a nonzero override; values must be
  finite and positive. The all-zero pair or omission retains each cell's
  isotropic `hs.a()` length.
- `ScatterParams` stores horizontal and vertical correlation lengths and uses
  them in its directional power spectral density.
- `AxisymmetricScatteringKernel` caches 65 incoming polar bins and fixed
  48-by-64 outgoing quadrature/CDF tables before workers start, so unequal
  anisotropic events no longer scan the full legacy sphere.
- `Scatterer` computes direction-dependent MFP and angular scattering weights
  using the incoming direction and `ECS.GetUp(location)` as the local vertical,
  then converts the sampled 3D direction and polarization into the existing
  incoming-phonon transform frame.
- Equal horizontal and vertical lengths retain the existing isotropic path;
  focused tests check PSD symmetry and transformed-trajectory equivalence.

Primary files: `scatparams.cpp`, `scatparams.hpp`, `axisymmetric_scattering.cpp`,
`axisymmetric_scattering.hpp`, `scatterers.cpp`, `scatterers.hpp`,
`sources.cpp`, `sources.hpp`, `phonons.cpp`, `model.cpp`, `model.hpp`,
`cmdline.cpp`, `cmdline.hpp`, and `docs/MANUAL.md`.

### Current Octave/gnuplot compatibility

- The do-scripts invoke Octave through a checked `--no-gui` wrapper and no
  longer force `GNUTERM=dumb`.
- Shared plotting helpers use current Octave cell expansion, current
  `colorbar` argument ordering, and figure annotations for paper-space text.
- `make test-plotting` verifies the headless plotting path and generated PDF
  artifact without making Octave/gnuplot a dependency of `make test`.

## Verification baseline

Run from the repository root:

```bash
make -j2
make test
make test-plotting
```

`make test` currently builds and runs:

- `tests/test_parallel_features.cpp`: default parameter checks, deterministic
  engine sequences, stream derivation, and seeded `ProbDist` selection.
- `tests/test_anisotropic_scattering.cpp`: isotropic/an-isotropic parameter
  state, directional PSD difference, and invalid-length rejection.
- `tests/test_axisymmetric_scattering.cpp`: reduced quadrature agreement,
  cache/CDF validity, local-frame conversion, MFP overrides, and no-deflection.
- `tests/test_report_reduction.cpp`: worker-context arithmetic and diagnostic
  reduction.
- `tests/test_seismometer_worker_bins.cpp`: worker-local seismometer capture.
- `tests/test_parallel_context_api.cpp`: explicit context overloads used by
  event generation and propagation.
- `tests/test_parallel_reproducibility.sh`: a process-level one-versus-two
  worker summary comparison with a fixed seed and reports disabled.
- `tests/test_do_capability_scripts.sh`: syntax and option checks for the
  capability-focused Lop Nor recipes.

`make test-plotting` separately runs
`tests/test_octave_plotting.sh`, which checks current Octave cell expansion,
gnuplot colorbar invocation, figure annotations, and PDF generation.

## Completed and merged into master

### Worker-local report reduction and chunked scheduling

- [x] Add `SimulationReportContext` for worker-local counters, diagnostics,
  and flattened seismometer bins.
- [x] Add context-aware event generation and propagation APIs while retaining
  compatibility overloads.
- [x] Replace per-phonon atomic scheduling with 256-phonon dynamic chunks.
- [x] Preserve random streams by deriving them from the absolute phonon index.
- [x] Merge contexts after worker join and keep enabled text reports serialized.
- [x] Add focused reduction, seismometer, API, and process-level regression
  tests.

### Performance comparison

- [x] Use the existing `do-lopnor-big.sh` and
  `do-lopnor-big-parallel.sh` recipes as matched serial/parallel workload
  definitions.
- [x] Measure one, two, and four workers on a 10M-phonon Lop Nor run.
- [x] Validate the merged result with one, four, and eight workers on the same
  10M-phonon Lop Nor workload.
- [x] Record machine/compiler/workload details and stable output summaries in
  `docs/DEVELOPMENT.md`.
- [x] Document that report ordering and floating-point reduction order are not
  equivalence criteria.

### Axisymmetric anisotropic scattering (completed in this branch)

- [x] Correct directional `GSATO` geometry around the local vertical symmetry
  axis with a stable parallel-incidence fallback.
- [x] Replace unequal-anisotropy full-sphere event scans with the immutable
  cached reduced kernel while preserving 3D trajectories and polarization.
- [x] Verify independent quadrature accuracy, seeded one/two-worker
  reproducibility, explicit MFP overrides, and no-deflection behavior.
- [x] Document the fixed numerical resolution and the two-axis/global-model
  limitations; retain analytic azimuth reduction and fully tensorial
  anisotropy as future work.

## Prioritized follow-up

### Task 1: Exercise command-line validation at the process boundary

**Files:**
- Create or extend: a focused command-line regression test.
- Modify: `cmdline.cpp`, `main.cpp`, or `Makefile` only if the test reveals an actual parsing/validation gap.
- Modify: `docs/MANUAL.md` for any corrected behavior.

**Cases:** `--workers=0`, a negative worker count, only one anisotropic length,
zero/negative/NaN/Infinity anisotropic lengths, equal anisotropic lengths, and
valid decimal/hexadecimal seeds. Keep error text stable enough for users but
avoid brittle tests that require unrelated banner text.

- [ ] Add process-boundary cases for invalid worker counts and anisotropy pairs.
- [ ] Add decimal and hexadecimal seed cases.
- [ ] Verify the all-zero anisotropy pair retains isotropic behavior.
- [ ] Run `make test`.

### Task 2: Establish continuous verification

**Files:**
- Create: CI configuration appropriate for the hosting forge.
- Modify: `Makefile` only when needed to expose the same local build/test entry points in CI.
- Modify: `docs/DEVELOPMENT.md` with supported compiler/platform assumptions.

**Acceptance:** CI runs `make -j2` and `make test` on every change. Add an
available sanitizer job when the runtime links correctly; document the
arm64-GCC ThreadSanitizer limitation rather than marking an unavailable job
as a false pass.

- [ ] Add the CI build and focused test jobs.
- [ ] Confirm CI uses the same Makefile commands as local development.
- [ ] Add sanitizer coverage only on a runner with a working runtime.

### Task 3: Decide the next parallelism layer

Before adding MPI or distributed-memory execution, write a design note covering
phonon ownership, seeded stream partitioning, report reduction, failure
propagation, and output ordering. The current chunked atomic worker boundary is
a good shared-memory seam but is not an MPI protocol.

## Domain limitations to keep visible

- Earth grids are currently compiled into `user.cpp`; file grid import is not implemented.
- The code uses class-static configuration for frequency, global scattering overrides, and several reporting/model settings; one active `Model` is the supported usage pattern.
- Multi-worker report ordering is nondeterministic even when per-phonon random streams are deterministic.
- Anisotropy is global per model invocation, not a per-cell CLI feature.
- Unequal anisotropy currently uses an axisymmetric two-correlation-length
  (2.5-D) kernel; three unequal axes and per-cell CLI symmetry axes are not
  supported.
- `--overridemfp` intentionally bypasses computed directional MFP behavior.
- ThreadSanitizer was not linkable in the arm64 GCC environment used for the current feature verification; rerun with a working sanitizer runtime before making stronger concurrency claims.
