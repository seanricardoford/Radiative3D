# Radiative3D agent instructions

## Start here

Before changing code, read these files in order:

1. `AGENTS.md` — repository rules.
2. `PLANS.md` — current baseline, completed work, and prioritized follow-up.
3. `docs/DEVELOPMENT.md` — architecture, runtime behavior, test coverage, and known limitations.
4. `docs/MANUAL.md` — user-facing command-line semantics.

Then inspect the worktree and branch before making assumptions:

```bash
git status --short --branch
git branch -vv
git log --oneline --decorate -8
```

Do not rewrite history, delete feature branches, or push to a remote unless the
user explicitly requests it. Create a new feature branch from the current
`master` tip for new implementation work.

## Build and test

- Use the repository Makefile for normal builds: `make` or `make -j2`.
- Run the focused regression suite with `make test`.
- Run both commands after changes to simulation, scattering, threading, output,
  command-line parsing, or build rules.
- The test programs are native C++ executables built by the Makefile; they use
  assertions and return nonzero on failure.
- `make clean` removes object files and `main` but keeps the generated floating
  point configuration. `make cleanall` also removes that generated file.
- Do not claim a build or test passes without running the command on the current
  tree and checking its exit status and output.

## Reproducibility and stochastic code

- Keep random state explicit in simulation APIs. Use `RandomEngine` from
  `probability.hpp`; do not introduce global `rand()` state.
- Preserve the stable phonon-index mapping used by `Model::SimulationThread()`.
  A phonon's stream is derived from `ModelParams::RandomSeed` and its index, so
  changing worker count must not change that phonon's random sequence.
- When changing numerical behavior, add or update a regression test before the
  implementation change where practical.
- Compare serial and multi-worker runs with the same `--seed` when changing
  propagation code. Per-phonon stochastic results are intended to be stable;
  report-line ordering is not guaranteed when multiple workers are active.
- Preserve documented defaults unless a behavior change is intentional, tested,
  and documented in `docs/MANUAL.md`.

## Parallel execution rules

- The simulation uses portable C++11 `std::thread` workers and an atomic
  remaining-phonon counter; it is not MPI and must not acquire an MPI dependency
  without a separate design decision.
- Shared reporting and seismometer accumulation are protected by the
  `DataReporter` mutex. Any new shared mutable state needs equivalent
  synchronization or a documented thread-local/reduction design.
- Probability distributions use lazy integration. Call or preserve the existing
  `PrepareForSimulation()` freeze point before worker threads start; do not add
  worker-time mutation to shared `ProbDist` objects.
- Validate worker counts at input boundaries. `--workers=1` is the serial
  reference; `--threads` is its compatibility alias.

## Anisotropic scattering rules

- `--scatter-horizontal` and `--scatter-vertical` are a pair. A nonzero override
  requires both values to be finite and positive; the all-zero pair disables the
  global override and is equivalent to omitting both options.
- Lengths use the model's length unit. Horizontal means perpendicular to the
  local model vertical/radial direction; vertical means parallel to it.
- With both values omitted, each cell uses its existing isotropic heterogeneity
  correlation length. Equal horizontal and vertical overrides must retain the
  isotropic result.
- The current implementation is global for a model. Do not imply that the CLI
  provides per-cell anisotropy unless that capability is added explicitly.
- Keep `--overridemfp` semantics intact: an explicit MFP override is diagnostic
  and takes precedence over the computed directional MFP.

## Scientific software practices

- Validate physical parameters at input boundaries and document units and
  coordinate frames in both code comments and the user manual.
- Keep model construction and simulation initialization order explicit. The
  `Model` constructor initializes class-static frequency, scattering, and
  seismometer state; the code historically expects one active `Model` at a time.
- Prefer small, focused regression tests over broad golden-output fixtures when
  testing stochastic or floating-point behavior.
- Record the random seed and worker count in diagnostics whenever adding a new
  execution path or output format.
- Do not silently change floating-point precision, compiler standard, or model
  units. The project currently builds as C++11 with GCC/G++ and `-pthread`.

## Documentation and scripts

- Update `docs/MANUAL.md` for every user-visible option or behavior change.
- Update `docs/DEVELOPMENT.md` when architecture, testing, limitations, or
  integration state changes.
- Keep `PLANS.md` current after completing or superseding roadmap work.
- Prefer direct project recipes and fixed repository-relative or documented
  absolute paths in project scripts.
- Do not add environment-variable overrides unless explicitly requested.
- Do not add explicit missing-file guard blocks unless the script truly needs
  recovery behavior.
- Do not parse logs to infer success when the underlying command can be checked
  directly.
- Keep scripts linear and readable over making them highly configurable.

## Verification caveat

AddressSanitizer/UndefinedBehaviorSanitizer verification has passed for the
current feature work. ThreadSanitizer could not link in the arm64 GCC
environment because required runtime symbols were unavailable. Treat that as a
toolchain limitation, not as evidence that concurrency is race-free; rerun it
when a working Clang/GCC sanitizer runtime is available.
