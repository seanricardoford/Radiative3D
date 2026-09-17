# Radiative3D development guide

This document is the durable technical handoff for contributors and coding
agents. It describes the current merged baseline, how the program is assembled,
where the recent parallel and anisotropic features live, how the optimized
worker path is verified, and what remains to be investigated.

## Project shape

Radiative3D is a C++11 seismic radiative-transport simulator. It traces
elastic P/S phonons through a compiled-in Earth model and uses stochastic
scattering to represent statistically described small-scale structure. The
repository has no external package manager or generated model importer: models
are selected from functions compiled into `user.cpp` and its `user_*_inc.cpp`
fragments.

The normal build is GCC/G++ plus GNU Make:

```bash
make -j2
make test
# Optional plotting regression; requires GNU Octave and gnuplot.
make test-plotting
```

The executable is `main`. Build artifacts and focused test executables are
ignored by Git. `docs/MANUAL.md` is the authoritative user-facing option
reference; the scripts beginning with `do-` are larger experiment recipes and
may require Bash, Octave, GMT, SAC, or other local tools.

## Source map

| Area | Files | Responsibility |
| --- | --- | --- |
| Entry point and CLI | `main.cpp`, `cmdline.cpp`, `cmdline.hpp`, `params.hpp` | Parse options, populate `ModelParams`/`MissionParams`, select actions. |
| Model orchestration | `model.cpp`, `model.hpp`, `grid.cpp`, `grid.hpp`, `media.cpp`, `media.hpp`, `media_cellface.*` | Build the compiled-in Earth model, create cells, initialize global simulation state, run phonons. |
| Geometry and materials | `geom_*`, `geom.hpp`, `elastic.*`, `ecs.*`, `tensors.hpp`, `raytype.hpp` | Coordinates, tessellations, elastic properties, coordinate transforms, and ray types. |
| Sources and propagation | `sources.*`, `events.*`, `phonons.*`, `raypath.*`, `rtcoef.*` | Generate event/scattered phonons, propagate through cells, handle reflection/refraction and ray state. |
| Randomness and discrete choices | `probability.*` | `RandomEngine` and lazy relative/cumulative probability distributions. |
| Scattering physics | `scatparams.*`, `scatterers.*` | Sato/Fehler scattering parameters, PSD/G functions, MFPs, outgoing direction/type sampling. |
| Output | `dataout.*` | Micro-reports, counters, seismometer bins, post-simulation traces and metadata. |
| User models | `user.cpp`, `user_*_inc.cpp` | Compiled-in model constructors and model-selection logic. |
| Tests | `tests/test_parallel_features.cpp`, `tests/test_anisotropic_scattering.cpp`, `tests/test_report_reduction.cpp`, `tests/test_seismometer_worker_bins.cpp`, `tests/test_parallel_context_api.cpp`, `tests/test_parallel_reproducibility.sh`, `tests/test_octave_plotting.sh` | Focused native, process-level, and plotting regressions invoked by `make test` or `make test-plotting`. |
| User docs and recipes | `README.md`, `docs/MANUAL.md`, `do-*.sh`, `scripts/`, `vis/` | Usage, experiment setup, post-processing, and visualization. |

## Octave and gnuplot workflows

The repository's `do-*.sh` recipes and visualization helpers invoke Octave in
non-GUI mode with `octave --no-gui`. They do not force the `GNUTERM` environment
variable; the plotting code selects gnuplot explicitly where headless output
needs it. This avoids the legacy `--no-window-system` and `GNUTERM=dumb`
combination, which is noisy and unreliable with current gnuplot releases.

The plotting helpers use current Octave cell-expansion syntax (`{:}`), pass
the colorbar location in the current argument order, and implement paper-space
labels as figure annotations rather than a second axes. The latter avoids
gnuplot multiplot stream warnings when writing PNG and PDF output. Run the
focused compatibility check with:

```bash
make test-plotting
```

The check requires GNU Octave and gnuplot but is intentionally separate from
`make test`, so native simulation tests remain usable on systems without the
optional visualization toolchain.

## Execution flow

1. `main()` creates `MissionParams` and `ModelParams`, suppresses reports by
   default, and parses `argv` through `CmdOpt`.
2. `process_option()` writes user values into `ModelParams`. Worker counts,
   seeds, and global horizontal/vertical scattering lengths are validated at
   this boundary or during `Model` construction.
3. `Model::Model()` initializes frequency-dependent static state, configures
   the global scattering override, constructs the selected grid and media
   cells, creates the event source/seismometers, and prepares probability
   distributions before simulation.
4. `Model::RunSimulation()` starts up to the requested number of C++11 worker
   threads. Each worker claims a 256-phonon chunk through an atomic next-index
   counter, derives a `RandomEngine` from each absolute phonon index, creates
   one event phonon, and calls `Phonon::Propagate(rng, context)`.
5. `Phonon::Propagate()` advances to cell boundaries or scattering events. It
   passes the explicit RNG into path-length, direction, polarization,
   reflection, and transmission choices. Event counters and seismometer bins
   are written to the worker's `SimulationReportContext`.
6. Workers join before `DataReporter` merges the contexts and before progress
   completion and post-simulation summaries are emitted. Enabled text reports
   are serialized by the reporter mutex; default report-off runs avoid that
   synchronization entirely.

The code historically expects one active `Model` at a time because several
subsystems use class-static configuration. Do not make model construction
concurrent without first designing ownership for those statics.

## Multi-worker implementation

### User interface

```text
-W, --workers=N       Number of shared-memory workers
--threads=N            Compatibility alias
--seed=VALUE           Base 64-bit random seed
```

`WorkerCount` defaults to detected hardware concurrency, with a one-worker
fallback when the platform reports zero. `RunSimulation()` caps the number of
workers at `NumPhonons`; `--workers=1` is the serial reference. Workers claim
256-phonon chunks, reducing atomic coordination from once per phonon to once
per chunk while retaining dynamic load balancing for variable propagation
costs.

### Deterministic stream design

`RandomEngine` wraps `std::mt19937_64`. The simulation uses:

```text
stream_seed = RandomEngine::SeedForStream(base_seed, phonon_index)
```

The phonon index is assigned atomically but the stream is derived from the
index, not from worker identity or scheduling order. Therefore worker-count
changes should not change the stochastic sequence assigned to an individual
phonon. This does not make output text order deterministic: multiple workers
can report events in different interleavings.

The explicit-RNG overloads are the preferred internal APIs. The no-argument
overloads use a thread-local compatibility engine for older callers and should
not be used to create shared global randomness in new code.

`ProbDist` lazily converts relative weights to cumulative weights. The
simulation calls `PrepareForSimulation()` on sources/scatterers before workers
start so worker threads only read the prepared tables.

`SimulationReportContext` owns the worker-local counters, invalid-diagnostic
bitfield, and flattened `[seismometer][time-bin]` report bins. The worker
context overloads of `DataReporter`, event generation, and propagation update
these objects without taking the reporter mutex. After joining, one
`MergeWorkerContext()` call per worker reduces them into the canonical
seismometer traces and summary counters. The legacy no-context APIs remain for
single-threaded callers and retain their synchronized behavior.

When micro-reports are enabled, only the text stream operation is protected by
`mReportMutex`; line order is intentionally nondeterministic with multiple
workers. `--reports=ALL_OFF` is the appropriate setting for performance
measurements. Seismometer and summary reductions use floating-point addition,
so future changes to reduction order may produce tiny round-off differences;
bit-for-bit equality is not a scientific requirement for this stochastic code.

## Anisotropic scattering implementation

The command-line options are:

```text
--scatter-horizontal=L_H
--scatter-vertical=L_V
```

Both must be supplied together for a nonzero override.
`ScatterParams::SetGlobalCorrelationLengths()` rejects non-finite or
non-positive nonzero values. The all-zero pair disables the global override and
leaves each cell's `Elastic::HetSpec::a()` as its isotropic correlation length.

For an override, `ScatterParams` stores `ah=L_H` and `av=L_V`. The directional
`GSATO(incoming, vertical, toa, ...)` path expresses the scattering geometry in
the local frame whose vertical axis is supplied by `ECS.GetUp(mLoc)`. The
`Scatterer` then uses directional G-values both for the exponential path
length (direction-dependent MFP) and for the sampled outgoing conversion/take-
off direction. `--overridemfp` intentionally takes precedence over computed
MFPs.

Equal horizontal and vertical lengths are the isotropic limit. The focused
anisotropic test checks that the isotropic PSD has equal horizontal/vertical
samples while unequal lengths produce a directional difference. It does not
yet validate a full propagated model statistically.

## Command-line and experiment workflow

For the complete option list, read `docs/MANUAL.md` or run:

```bash
make -j2
./main --help
```

For a real compiled-in model experiment, start from a `do-*.sh` recipe such as
`do-halfspace.sh`. Those scripts set model selection, source, seismometer,
frequency, phonon count, output directory, and post-processing parameters. For
parallel reproducibility, add explicit values such as:

```text
--workers=1 --seed=0x123456789abcdef0
--workers=4 --seed=0x123456789abcdef0
```

Use a small phonon count for development. Large example defaults (often `10M`
or more) are research runs, not unit tests.

The repository includes three focused Lop Nor recipes for exercising the merged
capabilities directly:

```bash
./do-lopnor-big.sh big-demo
./do-lopnor-parallel.sh parallel-demo
./do-lopnor-anistropic.sh anisotropic-demo
```

`do-lopnor-big.sh` runs the complete, visualized isotropic Lop Nor workflow as
a serial `10M`-phonon baseline. `do-lopnor-parallel.sh` uses the same model,
workload, and fixed seed with four explicit workers; compare the recorded run
times in their logs to estimate shared-memory speedup. The two speed recipes
are scientifically comparable because worker count is their intended runtime
variable. `do-lopnor-anistropic.sh` remains a smaller serial reference with a
fixed seed and a global ellipsoidal scattering override of 0.25 horizontally
and 1.25 in the local vertical direction. Its filename follows the requested
example name.

## Performance benchmark

The optimized worker/reduction path was measured with the Lop Nor compiled-in
model using the same 10,000,000-phonon workload, seed
`0x5eedc0de12345678`, 320 seismometers, and `--reports=INV`; plotting was
excluded. The benchmark ran on arm64 Darwin with Apple Clang 17.0.0. These are
single-run wall-clock measurements for this machine and are a comparison
record, not a portable performance guarantee:

| Workers | Wall time | Speedup vs. one worker |
| ---: | ---: | ---: |
| 1 | 428 s | 1.000x |
| 2 | 229 s | 1.869x |
| 4 | 130 s | 3.292x |

The previous implementation measured 635 s with one worker and 898 s with
four workers on the same workload. The new serial path was 32.6% shorter in
that comparison, and four workers provided a real 3.29x speedup. The output
summary was stable for all three worker counts: 9,789,917 loss-surface exits,
210,083 timeouts, and zero invalid phonons. The sampled seismometer files also
matched. Compare counters, traces, and other scientific outputs separately
from timing; report text ordering is not a valid equivalence check.

A fresh validation after merging the optimization into `master` used the same
machine, model, seed, 10M-phonon workload, and 320 seismometers. The serial
time is the second-resolution interval recorded by `do-lopnor-big.sh`; the
four-worker time is the full `do-lopnor-parallel.sh` wall time; the eight-worker
run used the exact parallel recipe command with only `--workers=8` substituted.
The serial figure-generation stage was excluded:

| Workers | Simulation wall time | Speedup vs. one worker |
| ---: | ---: | ---: |
| 1 | ~634 s | 1.00x |
| 4 | 191.46 s | 3.31x |
| 8 | 108.59 s | 5.84x |

The eight-worker run was a further 1.76x faster than four workers, reducing
simulation time by 43.3%. All three runs reported 9,789,917 loss-surface
exits, 210,083 timeouts, and zero invalid phonons. All 320 seismometer files
matched byte-for-byte across the comparisons. Generated benchmark directories
are ignored and may not be present in a fresh clone; the reproducible recipe
inputs are the checked-in do-scripts and the command recorded in their output
logs.

## Test coverage and gaps

`make test` currently runs five focused executables and two shell checks:

1. `test_parallel_features` checks default worker/seed parameters,
   repeatability of `RandomEngine`, stable stream derivation, and seeded
   `ProbDist` selection.
2. `test_anisotropic_scattering` checks anisotropic state, isotropic-limit PSD
   symmetry, unequal-length PSD difference, and rejection of an invalid length.
3. `test_report_reduction`, `test_seismometer_worker_bins`, and
   `test_parallel_context_api` cover worker-local reduction arithmetic,
   seismometer-bin accumulation, and the explicit context APIs.
4. `test_parallel_reproducibility.sh` runs the compiled-in model with 512
   phonons at one and two workers and compares stable summary counters for the
   same seed. It deliberately ignores report ordering.
5. `test_do_capability_scripts.sh` checks that the three capability-focused Lop
   Nor recipes are executable, syntactically valid, and pass their intended
   workload, worker, seed, and scattering options.

The current integration regression compares summary counters rather than full
floating-point traces. Process-boundary tests for malformed worker and
anisotropy options remain useful future coverage.

## Verification and toolchain notes

Fresh baseline verification is:

```bash
make -j2
make test
git status --short --branch
```

AddressSanitizer/UndefinedBehaviorSanitizer verification passed during the
feature work. ThreadSanitizer could not link in the arm64 GCC environment
because required runtime symbols were missing. This is an environment/toolchain
limitation; it is not a race-freedom proof. Repeat sanitizer checks with a
working runtime before making stronger concurrency claims.

## Current Git integration state

The optimized implementation is merged into local `master` at `8b9713d`. The
temporary `parallel-performance` branch and worktree were removed after the
fast-forward merge and merged-result test run. No remote feature branch was
created; local `master` is ahead of `origin/master` until a later session
pushes it. Check `git status --short --branch` and `git branch -vv` before
starting new work.

## Known limitations

- Grid-file import is not implemented; model geometry is compiled into the
  source tree.
- MPI/distributed-memory execution is not implemented.
- Class-static setup means one active `Model` is the supported pattern.
- Anisotropic correlation lengths are global per invocation, not per cell.
- Explicit MFP overrides bypass directional MFP calculation by design.
- Multi-worker report text is not ordered deterministically.
- Build/test infrastructure is GNU Make plus native executables; there is no
  CI configuration in the repository at this snapshot.
