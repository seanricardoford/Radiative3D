# Radiative3D capability plan

## `multithread`

The branch adds portable shared-memory parallelism. Worker count is selected with
`--workers` (or `--threads`) and defaults to detected hardware concurrency. The
`--seed` option controls deterministic per-phonon random streams. Stochastic APIs
receive an explicit `RandomEngine`; worker scheduling must not change a phonon's
random sequence. Shared reporting state is synchronized, and simulation progress
is emitted after workers join.

Future MPI work should build on the same stable phonon-index and reduction boundary
without adding an MPI dependency to this branch.

## `anisotropic-scattering`

Branch from the verified `multithread` tip. Add globally configured horizontal and
vertical correlation lengths, expressed relative to the model's local vertical or
radial axis. Extend the scattering power spectral density and angular sampling to
use the ellipsoidal correlation model while preserving isotropic behavior when the
two lengths are equal.

Required validation includes positive finite lengths, isotropic regression, and
directional statistical checks for horizontal and vertical propagation.
