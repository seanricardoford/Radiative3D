# Radiative3D Continuation and Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Preserve the merged multithreaded and anisotropic-scattering baseline while making future numerical, performance, and integration work reproducible and reviewable.

**Architecture:** `main.cpp` parses options into `ModelParams`; `Model` constructs the compiled-in Earth model and owns the simulation worker boundary; `Phonon`, `Scatterer`, `ScatterParams`, and `RandomEngine` implement stochastic propagation; `DataReporter` serializes shared reports and seismometer accumulation. New work should extend these boundaries without introducing hidden global state or worker-time mutation of shared probability tables.

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

At the time this handoff was written, `master` contains the feature-code tip
`d4b8039` and both requested capabilities. Documentation commits may be newer
than that feature-code tip:

```text
ae81634  origin/master
  |
aeb3cb8  initial thread coordination
  |
4386340  deterministic multi-worker simulation support
  |
d4b8039  ellipsoidal anisotropic scattering (feature-code tip)
```

The local `multithread` branch points to `4386340`; the local
`anisotropic-scattering` branch points to `d4b8039`. They are retained as
historical feature pointers. Local `master` is ahead of `origin/master`; the
merge was local and has not been pushed unless a later session changes that
state. Always re-check Git state because this snapshot will age.

## Completed capabilities

### Multi-worker simulation

- `--workers=N`, `-W N`, and `--threads=N` select shared-memory worker count.
- The default is `std::thread::hardware_concurrency()`, falling back to one.
- `--seed=VALUE` sets the base 64-bit seed; decimal and `0x` hexadecimal forms
  are accepted by the existing command-line parser.
- Each phonon gets a stream derived from the base seed and stable phonon index
  through `RandomEngine::SeedForStream()`.
- `Model::RunSimulation()` caps workers at the requested phonon count, uses an
  atomic remaining counter, and joins all workers before post-simulation output.
- Stochastic APIs have explicit `RandomEngine&` overloads; no-argument overloads
  remain as thread-local compatibility wrappers.
- `DataReporter` protects reports, counters, and seismometer accumulation with a
  mutex. Lazy probability integration is completed before workers start.

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
- `Scatterer` computes direction-dependent MFP and angular scattering weights
  using the incoming direction and `ECS.GetUp(location)` as the local vertical.
- Equal horizontal and vertical lengths are intended to recover isotropic
  behavior; the focused test checks the corresponding PSD symmetry.

Primary files: `scatparams.cpp`, `scatparams.hpp`, `scatterers.cpp`,
`scatterers.hpp`, `model.cpp`, `model.hpp`, `cmdline.cpp`, `cmdline.hpp`, and
`docs/MANUAL.md`.

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
- `tests/test_octave_plotting.sh`: current Octave cell expansion, gnuplot
  colorbar invocation, figure annotations, and PDF generation.

The current unit suite does not replace an end-to-end simulation comparison.
The next parallel regression should run the same small compiled-in model with
`--workers=1` and `--workers=2` (or more), the same `--seed`, reports disabled,
and compare stable summary/seismometer results. Do not compare raw report line
order from multi-worker runs.

## Prioritized follow-up

### Task 1: Add a small end-to-end serial/parallel regression

**Files:**
- Create: `tests/test_simulation_reproducibility.cpp` or a repository-supported equivalent integration test.
- Modify: `Makefile` to build and run it through `make test`.
- Modify: `docs/DEVELOPMENT.md` with the exact invocation and comparison rule.

**Acceptance:** A small deterministic compiled-in model runs with one and
multiple workers using the same seed; reports are disabled; the test compares
stable simulation outputs or counters and does not depend on scheduling order.

- [ ] Add the smallest supported compiled-in model invocation.
- [ ] Run it with `--workers=1 --seed=0x123456789abcdef0`.
- [ ] Run it again with the same seed and multiple workers.
- [ ] Compare stable counters or seismometer values, not report ordering.
- [ ] Run `make test`.

### Task 2: Exercise command-line validation at the process boundary

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

### Task 3: Add reproducible performance measurements

**Files:**
- Create: a short documented benchmark recipe under `docs/` or `scripts/`.
- Modify: `docs/DEVELOPMENT.md` with machine, compiler, model, phonon count,
  seed, worker counts, and timing interpretation.

**Acceptance:** The recipe uses fixed repository-relative paths, does not alter
scientific defaults, and reports throughput for `--workers=1` and at least
two larger worker counts. It must distinguish wall-clock speedup from
stochastic-output equivalence.

- [ ] Record compiler, machine, model, phonon count, seed, and worker counts.
- [ ] Measure serial and at least two multi-worker configurations.
- [ ] Record output-equivalence checks separately from timing results.

### Task 4: Establish continuous verification

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

### Task 5: Decide the next parallelism layer

Before adding MPI or distributed-memory execution, write a design note covering
phonon ownership, seeded stream partitioning, report reduction, failure
propagation, and output ordering. The current atomic worker boundary is a good
shared-memory seam but is not an MPI protocol.

## Domain limitations to keep visible

- Earth grids are currently compiled into `user.cpp`; file grid import is not implemented.
- The code uses class-static configuration for frequency, global scattering overrides, and several reporting/model settings; one active `Model` is the supported usage pattern.
- Multi-worker report ordering is nondeterministic even when per-phonon random streams are deterministic.
- Anisotropy is global per model invocation, not a per-cell CLI feature.
- `--overridemfp` intentionally bypasses computed directional MFP behavior.
- ThreadSanitizer was not linkable in the arm64 GCC environment used for the current feature verification; rerun with a working sanitizer runtime before making stronger concurrency claims.
