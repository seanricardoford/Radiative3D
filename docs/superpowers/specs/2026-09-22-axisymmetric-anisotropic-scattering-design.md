# Axisymmetric Anisotropic Scattering Design

**Date:** 2026-09-22  
**Status:** Approved for implementation  
**Branch:** `axisymmetric-anisotropic-scattering`

## Problem

The current anisotropic implementation evaluates the direction-dependent Sato
and Fehler scattering functions over the full degree-9 take-off-angle sphere
for every scattering event. That sphere contains 5,242,880 directions. The
implementation therefore preserves the desired anisotropic formula but is not
usable for realistic unequal horizontal and vertical correlation lengths.

The approved solution must retain fully three-dimensional phonon trajectories
and polarization transforms while exploiting the rotational symmetry of a
medium whose horizontal correlation lengths are equal. It must preserve the
existing isotropic result when the horizontal and vertical lengths are equal,
and it must continue to honor explicit MFP overrides.

## Scientific model and nomenclature

The local model vertical/radial unit vector is the symmetry axis, `up`. The
correlation lengths are:

* `ah`: horizontal correlation length, perpendicular to `up`;
* `av`: vertical correlation length, parallel to `up`.

The directional von Karman power spectral density is already represented by

```text
Phi(q) = C / (1 + ah^2 |q_h|^2 + av^2 q_v^2)^(kappa + 3/2)
```

where `q_v = q dot up` and `|q_h|^2 = |q|^2 - q_v^2`. This is the
axisymmetric (`a=b`) specialization of the textured-medium spectrum described
in Margerin, *Tectonophysics* 416 (2006), including the paper's exponential
case when the von Karman parameter is one half.

For an incoming direction `p`, the kernel depends on its polar angle relative
to `up`, represented by `mu_in = p dot up`, but not on its absolute azimuth.
The outgoing relative direction still has two coordinates:

* `psi`: scattering angle relative to the incoming trajectory;
* `zeta`: relative azimuth around the incoming trajectory, measured from the
  incoming/up meridian.

The outgoing direction is sampled in this local relative frame and passed to
the existing `Phonon::Transform()` method. Thus the propagation path remains
fully 3D even though the scattering kernel is reduced to axisymmetric
coordinates.

The kernel's relative azimuth is measured from the incoming/up meridian. The
existing `Phonon::Transform()` API instead interprets a relative phonon's
azimuth in the incoming polarization frame. The anisotropic scattering API
must therefore receive the incoming polarization angle and convert each
sampled direction and S polarization vector from the symmetry-axis meridian
into the transform frame before constructing the relative phonon. This keeps
the anisotropic kernel tied to `up` and prevents the current S polarization
orientation from silently rotating the statistical medium.

## Goals and non-goals

### Goals

1. Remove full `nTOA` angular scans from the anisotropic simulation hot path.
2. Compute direction-dependent mean free paths from the reduced solid-angle
   integral.
3. Sample P/P, P/S, S/P, and S/S conversion and outgoing directions from the
   same reduced kernel used for the MFP calculation.
4. Preserve S/S polarization computed by `GSATO`.
5. Freeze all shared tables before worker threads start.
6. Validate the reduced calculation against an independent high-resolution
   reference.

### Non-goals

* A fully anisotropic three-axis correlation tensor with no symmetry axis.
* An anisotropic elastic stiffness tensor or wave-speed birefringence.
* Per-cell symmetry-axis or correlation-length command-line options.
* An analytic special-function implementation of every azimuthal integral in
  the first version.
* Changes to the isotropic scattering algorithm or its stochastic defaults.

## Approaches considered

### Cached axisymmetric quadrature — selected

Use a one-dimensional incoming-angle table and a tensor-product quadrature over
`cos(psi)` and `zeta`. For each incoming-angle bin, cache the directional
inverse MFP, conversion totals, and conditional CDFs for outgoing nodes. At
runtime, interpolate between neighboring incoming bins and sample from the
cached distributions.

This directly represents the axisymmetric radiative-transfer kernel, gives
deterministic accuracy, and makes event-time work logarithmic in the reduced
table size. It also keeps the existing 3D direction and polarization transform
untouched.

### Analytic azimuthal reduction

Integrate the azimuth analytically using the quadratic `q_v` dependence and
the low-order trigonometric factors in the Sato radiation patterns, then use a
one-dimensional polar quadrature. This could be faster, but conditional
azimuth sampling and S/S polarization would still need a numerical inversion or
rejection step. Special-function edge cases would add scientific and testing
risk without being necessary to remove the current bottleneck.

### Runtime rejection sampling

Draw candidates from a simple sphere distribution and accept them according to
the anisotropic kernel. This is easy to prototype, but acceptance becomes poor
for elongated or flattened scatterers, and the directional MFP still requires a
separate integration. It does not provide a reliable performance floor.

## Architecture

### `ScatterParams`

Keep `PowerSpectralDensity(q, up)` and the directional `GSATO` interface as the
single source of scattering physics. Add comments documenting `ah`, `av`,
`q_h`, `q_v`, and the local symmetry-axis convention. No caller will construct
the reduced kernel by duplicating the PSD or radiation-pattern equations.

The directional `GSATO(incoming, up, toa, ...)` implementation will construct a
stable orthonormal basis from `incoming` and `up`: the first transverse axis is
the incoming/up meridian and the second is its right-handed perpendicular. It
will use that basis to form the outgoing wavevector in the PSD calculation; it
will not assume that the simulation-frame global Z axis is the local symmetry
axis. Near parallel incoming/up directions, a deterministic perpendicular
fallback will avoid division by a vanishing transverse projection.

### `Scatterer`

Add a private axisymmetric-kernel cache containing:

* incoming-angle nodes or bins over `mu_in` in `[-1, 1]`;
* per-bin inverse MFPs for P and S;
* per-bin unnormalized totals for the four conversion types;
* per-bin cumulative direction tables for each conversion type;
* shared outgoing relative nodes `(psi, zeta)` and their S/S polarization
  angles.

The cache is constructed in the `Scatterer` constructor for unequal lengths and
is immutable during simulation. The legacy take-off-angle distributions remain
available for isotropic sources and scatterers; the anisotropic `Scatterer`
uses the compact cache instead and does not scan or sample full-sphere arrays.
`PhononSource` will therefore gain an internal compact-construction mode that
leaves the existing default/event-source allocation and behavior unchanged.
In compact mode the anisotropic scatterer does not populate `mPDists`; its
directional methods and diagnostics use the axisymmetric cache directly.

The first numerical configuration will use fixed private resolutions selected
by convergence tests: a moderate incoming-angle table, Gauss-Legendre nodes in
`cos(psi)`, and periodic trapezoidal nodes in `zeta`. The quadrature weights
will be normalized so their sum represents a unit average over `4*pi`, matching
the existing MFP convention. The implementation will not expose resolution
as a user option.

### Runtime path length

For unequal anisotropy and no explicit MFP override:

1. Compute `mu_in = incoming dot vertical`.
2. Find the neighboring cached incoming-angle bins.
3. Interpolate the inverse MFP between those bins.
4. Return `-log(1-U) / inverse_mfp`.

When `--overridemfp` is active, retain the existing explicit P/S MFP values
and bypass the computed directional MFP. The no-deflection override continues
to return an unmodified relative phonon.

### Runtime scattering direction

For unequal anisotropy and normal deflection:

1. Use the same neighboring incoming bins and blend their unnormalized
   conversion totals.
2. Select the allowed conversion for the incoming P or S type.
3. Select one of the neighboring cached direction CDFs in proportion to that
   bin's contribution to the selected conversion.
4. Sample a reduced quadrature node from that CDF.
5. Construct `Phonon(S2::ThetaPhi(psi, zeta), output_type)` and assign the
   cached S/S polarization angle when applicable, after converting both
   direction and polarization from the symmetry-axis meridian into the current
   incoming polarization frame.
6. Let `Phonon::Transform()` express the relative direction in the current
   global 3D frame.

This mixture interpolation is equivalent to linearly interpolating the
unnormalized conditional kernel while avoiding per-event CDF construction.
Zero-total or numerically degenerate bins fall back to a compact isotropic
reference sampler rather than producing an invalid phonon. The anisotropic
no-geometry diagnostic overload will likewise use a defined cached reference
incoming direction; it will never index an unpopulated `mPDists` array.

When `ah == av`, retain the current isotropic path exactly. This is both a
behavioral compatibility requirement and a useful regression oracle.

## Initialization and thread safety

All quadrature nodes, radiation values, totals, CDFs, and diagnostic summaries
must be built before `Model::RunSimulation()` starts worker threads. No worker
may lazily populate or mutate a shared axisymmetric table. The explicit
`RandomEngine&` APIs remain the only simulation path for random choices.

The cache is owned by the `Scatterer` and shared read-only by cells that reuse
matching scatter parameters. No new global random state or worker lock is
allowed.

## Diagnostics and compatibility

The existing scattering statistics remain available. Since anisotropic MFP and
dipole values are functions of incoming angle, the scalar values printed by
the legacy report will be defined as solid-angle averages over incoming
directions, while the cache retains the directional values used in simulation.
The report should identify the correlation lengths as before and document this
summary convention in the development guide.

The existing scatter-pattern diagnostic will use a defined reference incoming
direction and local vertical when it exercises the anisotropic sampler. It must
not accidentally call the isotropic no-geometry overload with an uninitialized
full-sphere distribution.

The geometry-aware scattering overload will carry the incoming polarization
angle required by the meridian-to-transform-frame conversion. Existing
no-geometry overloads remain compatibility wrappers using zero polarization.

## Verification strategy

### Numerical kernel tests

Extend the anisotropic test coverage to verify:

* equal lengths retain isotropic PSD and sampler behavior;
* unequal lengths produce distinct vertical and horizontal PSD values;
* rotating the incoming direction around the symmetry axis leaves directional
  MFP and conversion totals unchanged;
* rotating the incoming S polarization while holding the incoming direction and
  local vertical fixed changes only the transform-frame representation, not the
  sampled physical outgoing direction or S polarization vector;
* reduced quadrature MFPs and conversion totals agree with an independent,
  higher-resolution two-angle reference within documented tolerances;
* increasing quadrature resolution produces a converged result for the Lop Nor
  parameters and a more strongly anisotropic stress case;
* all cached CDFs are monotone, finite, and normalized.

### Propagation and reproducibility tests

Run a small seeded anisotropic simulation through the existing process-level
path and verify completion, zero invalid-phonon diagnostics, and a finite
summary. Compare one-worker and multi-worker runs with the same seed using the
existing reproducibility conventions. Confirm that explicit MFP overrides and
no-deflection behavior remain intact.

### Performance test

Repeat the 100-phonon unequal-anisotropy probe that previously stalled during
event-time full-sphere scans. The test should complete without scanning the
5,242,880-direction table per event. Record construction and simulation timing
as diagnostic evidence rather than making a fragile machine-specific wall-time
assertion.

### Required repository verification

After implementation, run:

```text
make -j2
make test
make test-plotting
```

The plotting check remains optional when its external tools are unavailable,
but any unavailable toolchain must be reported explicitly.

## Documentation and integration

Update `docs/DEVELOPMENT.md` with the axisymmetric cache, quadrature variables,
thread-freeze point, diagnostic averaging convention, and known numerical
limits. Update `docs/MANUAL.md` only for user-visible semantics or explicitly
documented performance/accuracy behavior. Keep `PLANS.md` current as the
implementation milestones are completed.

The feature branch will be merged back into the pre-existing
`anisotropic-lopnor-script` branch after review and verification. The temporary
feature branch will then be deleted; no history will be rewritten.
