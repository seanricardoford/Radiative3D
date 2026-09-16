# Radiative3D development guide

This document is the durable technical handoff for contributors and coding
agents. It describes the current merged baseline, how the program is assembled,
where the recent parallel and anisotropic features live, and what remains to be
verified.

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
| Tests | `tests/test_parallel_features.cpp`, `tests/test_anisotropic_scattering.cpp` | Focused assertion-based regression executables invoked by `make test`. |
| User docs and recipes | `README.md`, `docs/MANUAL.md`, `do-*.sh`, `scripts/`, `vis/` | Usage, experiment setup, post-processing, and visualization. |

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
   threads. Each worker claims a phonon index through the atomic remaining
   counter, creates a `RandomEngine` for that index, generates one event
   phonon, and calls `Phonon::Propagate(rng)`.
5. `Phonon::Propagate()` advances to cell boundaries or scattering events. It
   passes the explicit RNG into path-length, direction, polarization,
   reflection, and transmission choices. `DataReporter` receives event calls
   during propagation.
6. Workers join before progress completion and post-simulation summaries are
   emitted. Report and seismometer writes are protected by the reporter mutex.

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
workers at `NumPhonons`; `--workers=1` is the serial reference.

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

`DataReporter` uses `mReportMutex` around report streams, counters, and
seismometer accumulation. Any new shared accumulator must follow the same
ownership rule or use a worker-local reduction followed by a join-time merge.

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

## Test coverage and gaps

`make test` currently runs two focused executables:

1. `test_parallel_features` checks default worker/seed parameters,
   repeatability of `RandomEngine`, stable stream derivation, and seeded
   `ProbDist` selection.
2. `test_anisotropic_scattering` checks anisotropic state, isotropic-limit PSD
   symmetry, unequal-length PSD difference, and rejection of an invalid length.

The repository still needs an automated end-to-end regression that runs a small
compiled-in model with one and multiple workers and compares stable summary or
seismometer results. Report line order must not be used as the comparison
because it is schedule-dependent. CLI process-boundary tests for malformed
worker/anisotropy options would also strengthen coverage.

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

The requested feature branches were fast-forwarded into local `master`; the
feature-code tip is `d4b8039` and documentation commits may follow it:

```text
origin/master -> aeb3cb8 -> 4386340 -> d4b8039 -> documentation commits (master)
                 multithread     anisotropic-scattering
```

The exact remote pointers can change. At the handoff snapshot, local `master`
was three commits ahead of `origin/master`, and neither feature branch was
deleted. Check `git branch -vv` before pushing or integrating more work.

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
