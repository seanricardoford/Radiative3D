# Radiative3D development guidance

## Build and test

- Use the repository Makefile for normal builds: `make`.
- Run the focused regression suite with `make test`.
- Keep stochastic tests explicitly seeded.
- Verify serial and multi-worker runs with the same seed when changing simulation code.

## Script style

- Prefer direct project recipes and fixed repository-relative or documented absolute paths.
- Do not add environment-variable overrides unless explicitly requested.
- Do not add explicit missing-file guard blocks unless recovery is genuinely required.
- Do not parse logs to infer success when the underlying command can be checked directly.
- Keep scripts linear and readable rather than highly configurable.

## Scientific-computing practices

- Keep random state explicit in simulation APIs; do not use global `rand()` state.
- Preserve backward-compatible defaults unless a documented option changes behavior.
- Validate physical parameters at input boundaries and report units and coordinate frames.
- Add a regression test before changing numerical behavior.
- Record the random seed and parallel worker count in run diagnostics.
