# Take-off-Angle Quadrature Bias

Radiative3D's anisotropic directional mean free paths and dipole moments --
both the scatterer-dump report values and the rates the simulation actually
realizes -- are averaged over the run's icosahedral take-off-angle set with
equal weight per direction. Per design ruling R9 this bias is accepted, not
fixed. `doc/MANUAL.md`'s Validity-limits section gives the user-facing
summary; this note records the mechanism and a plan for a future fix.

## Mechanism

Each take-off-angle direction is the center of one icosahedral subdivision
triangle (`S2::TesselSphere`, `TESS_ICO`). The code averages over directions
with equal weights, but the triangles' solid angles vary from 0.928× to
1.207× the mean (standard deviation 8.6%), identically at degrees 6–9:
further subdivision refines the triangles without equalizing them. The
equal-weight average therefore carries a fixed bias that does not shrink
with degree.

## Isotropic effect

The offset is the same for every incidence direction and consistent between
rate and pattern. Examples: Lop Nor crust S +1.11%, P −0.04%; `a` = 2 km P
+0.62%, S +0.78%.

## Anisotropic effect

For an anisotropic region the bias also depends on direction and
polarization, which can distort the anisotropy signal itself. For the
(2, 2, 0.5) km crust (equal weight against area weight):

| Quantity | Equal weight | Area weight |
|---|---:|---:|
| S90SH | 397.8 km | 374.7 km |
| S90SV | 433.9 km | 442.8 km |
| P90 | 1178.8 km | 1203.2 km |
| P0 | 1610.4 km | 1602.7 km |
| S0 | 759.9 km | 754.7 km |

The area-weighted values match the original design's independent
cross-check table, and Task 6 of the 2026-10-16 plan got 374.7 km from a
200k-point Fibonacci sphere.

## Frame orientation for P incidence

Incoming P has no physical polarization-rotation angle, but a sampler
must still pick some frame to build its take-off-angle grid in, and the
equal-weight quadrature's bias makes the realized P rate depend on that
otherwise-meaningless choice. Scanning the frame angle about the ray, for
a Lop Nor r=10-like model (unchanged at `--toa-degree` 5, 7, and 8):

- Crust at 60 degrees from vertical: MFP 1442–1518 km across frame angles.
- Mantle at 89 degrees from vertical: MFP 742–826 km across frame angles.

Both the axisymmetric and general samplers use the canonical orientation
(`Scatterer::CanonicalPPolarization()`), the highest-rate (shortest-MFP)
end of each range above, so their P results agree, including for an
already-scattered P phonon that carries an arbitrary polarization-
rotation angle from its previous scattering event. Before task 6b of the
2026-10-17-axisymmetric-sampler-performance plan (ruling R10), the
general sampler used that carried angle to build its take-off-angle grid
instead of the canonical one, so an already-scattered P phonon's realized
rate could land anywhere in the range above; `test_p_frame_consistency`
(`tests/anisotropy_tests.cpp`) checks that the two samplers now agree at
a representative non-canonical carried angle. An area-weighted quadrature
(see below) would remove the orientation dependence itself, rather than
just pin both samplers to its highest point.

## Provenance

The tables above come from an independent re-implementation of the sampler
maths (a numpy prototype, not in this repository). The 397.8 km equal-weight
S90SH figure is reproduced in-repo by `test_report_directional_mfps`
(`tests/anisotropy_tests.cpp`). The frame-orientation ranges above come
from the same investigation's diagnostic scan over 24 frame angles
(2026-10-17-axisymmetric-sampler-performance, `debug-scratch.patch`'s
`scratch_p_rate_vs_frame`, not in this repository); `test_p_frame_consistency`
covers the resulting agreement between samplers, not a re-derivation of the
ranges themselves.

## Affected code

A future area-weighted fix would touch:

- `PopulateProbDists` (`scatterers.cpp`): weight each entry by
  `w_k = Ω_k·N/4π`.
- `PopulateWholeProbs`, `ComputeMFPs`, `ComputeDipoles` (`scatterers.cpp`).
- `GetRandomScatteredRelativePhonon` (`scatterers.cpp`), which needs no
  change once the tables are weighted.
- `IntegrateRateAndDipole`, `ComputeReportStats`, `ReportDirs`
  (`scatterers.cpp`).
- `AxiProposal` (`aniso_axi.cpp`): the cell sums inherit the weights, and
  the bounds on R do not change.
- the general-sampler majorant, computed in the `Scatterer` constructor
  (`scatterers.cpp`).
- the `ShearDislocation` radiation tables (`events.cpp`).
- `S2::TesselSphere` (`geom_s2.cpp`/`geom_s2.hpp`), which must expose
  per-direction solid angles. Today the triangles are internal to
  `geom_s2.cpp`.

## Consequences

Every output changes, isotropic included, which breaks byte-identity with
`aaec540`. The change should therefore be opt-in (e.g.
`--toa-quadrature=area`, default `equal`). The forced-isotropic and seven-
vs five-value equivalences documented in `doc/MANUAL.md` still hold within
either mode.

## Validation for a future fix

- weights sum to 4π within 1e-12;
- isotropic MFPs match an exact 2-D integral within 1e-4;
- the anisotropic values above are reproduced;
- the sampler χ² tests pass against the weighted targets;
- a scenario-level report shows how much outputs move.
