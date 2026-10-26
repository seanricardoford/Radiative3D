# Take-off-Angle Quadrature Bias

Radiative3D's directional mean free paths and dipole moments -- both the
scatterer-dump report values and the rates the simulation actually
realizes -- are averaged over the run's icosahedral take-off-angle set.
`doc/MANUAL.md`'s Validity-limits section gives the user-facing summary;
this note records the mechanism, the two available quadratures' effects,
and how each is validated.

## Status

Opt-in, on branch `quadrature` (cut from `anisotropy` at commit
`319c8dd`): `--toa-quadrature=equal|area`, default `equal`. `equal` mode
is unchanged -- byte-identical to the behavior this note has always
described, for isotropic and anisotropic runs alike
(`.superpowers/sdd/2026-10-25-toa-area-quadrature/global-constraints.md`).
Design ruling R9 (the `2026-10-16`/`2026-10-17` specs) accepted the
mechanism below as a known, documented bias rather than fixing it; this
branch implements the fix found feasible by the `2026-10-25` design
spec's spike as the new `area` mode
(`docs/superpowers/specs/2026-10-25-toa-area-quadrature-design.md`,
decisions 1-9). `area` is not the default -- that remains a separate
decision for project review (the same spec's non-goals).

## Mechanism

Each take-off-angle direction is the center of one icosahedral
subdivision triangle (`S2::TesselSphere`, `TESS_ICO`). In `equal` mode
the code averages over directions with equal weights, but the triangles'
solid angles vary from 0.928× to 1.207× the mean (standard deviation
8.6%), identically at degrees 6–9: further subdivision refines the
triangles without equalizing them. The equal-weight average therefore
carries a fixed bias that does not shrink with degree.

Equal weighting also distorts *which* directions are actually sampled or
emitted, separately from that rate bias: every node gets the same vote
regardless of the solid angle it represents, so a node's effective share
of scattered or radiated phonons is modulated by mean(Ω)/Ω_k relative to
a uniform density on the sphere -- a factor that ranges from 0.829 to
1.078 across the icosahedral node set at degrees 6-8
(`docs/toa-quadrature-runs-2026-10/toa_quadrature_solid_angles.txt`,
computed as 1/max_over_mean to 1/min_over_mean). `area` mode's
w_k = Ω_k·N/4π weighting is exactly proportional to solid angle, so it
removes both the rate bias above and this direction-density distortion;
see Isotropic effect, Anisotropic effect, and Event source, below.

## Isotropic effect

The offset is the same for every incidence direction and consistent
between rate and pattern. Examples: Lop Nor crust S +1.11%, P −0.04%;
`a` = 2 km P +0.62%, S +0.78%.

`area_iso_mfp_exact` (`tests/anisotropy_tests.cpp`) measures this same
kind of offset directly against an exact 2-D numerical integral, for
both quadratures, at take-off-angle degrees 5-7
(`docs/toa-quadrature-runs-2026-10/toa_quadrature_iso.txt`):

| a | kappa | degree | P equal | P area | S equal | S area |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.3 | 5 | −0.0442% | −0.000049% | +1.0893% | +0.000543% |
| 0.5 | 0.3 | 6 | −0.0441% | −0.000012% | +1.0878% | +0.000135% |
| 0.5 | 0.3 | 7 | −0.0441% | −0.000003% | +1.0874% | +0.000034% |
| 2 | 0.5 | 5 | +0.6221% | +0.005743% | +0.7928% | +0.006867% |
| 2 | 0.5 | 6 | +0.6180% | +0.001435% | +0.7864% | +0.001728% |
| 2 | 0.5 | 7 | +0.6169% | +0.000359% | +0.7848% | +0.000433% |

(relative error against the exact integral, per ray type.) The test
checks the area column at (a=0.5, kappa=0.3, degree 5) and
(a=2, kappa=0.5, degree 6) -- requiring ≤1e-4 absolute relative error at
both (actual: ≤5.4e-6 and ≤1.7e-5 respectively) -- and additionally
requires the *equal* column to exceed 3e-3 at the second point (actual
0.618% P, 0.786% S), confirming the two quadratures are actually
discriminated rather than both happening to be small. Area weighting's
error keeps falling with degree because, unlike the fixed per-direction
bias equal weighting leaves in place, it is a pure discretization
residual of the Ω_k computation itself.

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
200k-point Fibonacci sphere. `area_report_directional_mfps`
(`tests/anisotropy_tests.cpp`) now reproduces this table in-repo, at
take-off-angle degree 6, checking P0/P90/S0/S90SV/S90SH against these
values within 2e-4 relative; the full 8-state comparison (also P45,
S45SV, S45SH, which have no independent reference here, same as the
equal-mode table above) is logged to
`docs/toa-quadrature-runs-2026-10/toa_quadrature_report.txt`.

## Frame orientation for P incidence

Incoming P has no physical polarization-rotation angle, but a sampler
must still pick some frame to build its take-off-angle grid in, and the
equal-weight quadrature's bias makes the realized P rate depend on that
otherwise-meaningless choice. Scanning the frame angle about the ray,
for a Lop Nor r=10-like model (unchanged at `--toa-degree` 5, 7, and 8):

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
rate could land anywhere in the range above; `p_frame_consistency`
(`tests/anisotropy_tests.cpp`) checks that the two samplers now agree at
a representative non-canonical carried angle.

Area weighting does not remove this dependence outright, but reduces it
to a discretization residual that shrinks roughly 4x with each
additional take-off-angle degree, measured at the same two cases (crust
60 degrees, mantle 89 degrees; design spec's feasibility-spike table,
`docs/superpowers/specs/2026-10-25-toa-area-quadrature-design.md`):

| Degree | Crust (60°) area spread | Mantle (89°) area spread |
|---:|---:|---:|
| 5 | 0.054% | 0.117% |
| 6 | 0.013% | 0.028% |
| 7 | 0.003% | 0.007% |

(equal-mode spreads stay flat at 5.28% and 11.2% respectively, at every
degree -- the same spike table.) `area_p_frame_residual`
(`tests/anisotropy_tests.cpp`) reproduces the degree-6 row in-repo,
logging every sampled frame angle's equal- and area-weighted MFP, plus
the frame-independent exact value, to
`docs/toa-quadrature-runs-2026-10/toa_quadrature_pframe.txt`: the
measured spread (maximum/minimum − 1, across the 24 sampled angles) is
0.013% (crust) and 0.028% (mantle) under area weighting, against 5.28%
and 11.21% under equal weighting -- matching the spike figures above and
the test's `< 1e-3` / `> 4e-2` thresholds; every area-mode value also
lands within 2e-4 of the frame-independent exact rate.

## Event source

Event-source radiation tables (`ShearDislocation`, `events.cpp`) are
weighted the same way as the scatterer tables: in `area` mode, each take-
off direction's stored P/SH/SV squared-amplitude is multiplied by
`TOAWeight(k)` before being stored, in a separate code path from the
unweighted `equal`-mode arithmetic (so GCC's default
`-ffp-contract=fast` cannot fold the weight into an equal-mode sum).
`area_event_source_tables` (`tests/anisotropy_tests.cpp`) confirms this
entry-by-entry, for two moment tensors across all three ray types, to
within 1e-11 of each ray type's table total.

Icosahedral symmetry integrates these radiation patterns' low-order
polynomials exactly under either quadrature -- both modes share the same
direction set, so only the per-direction weight differs -- so the event
source's overall partition between ray types does not move: the same
test confirms the P-fraction of the total radiated energy is 0.4 within
1e-12 in both modes, for both tested moment tensors. Only the
*distribution* of radiated phonons across direction changes, by the same
mean(Ω)/Ω_k modulation described in Mechanism, above.

## Implementation

`area` mode's implementation touches:

- `geom_s2.hpp`/`geom_s2.cpp` -- `S2::TesselSphere` gained an optional
  `withSolidAngles` constructor argument and a parallel `SolidAngles()`
  vector; `Triangle::populate_lists()` computes each leaf's Ω via the Van
  Oosterom-Strackee formula, `2*atan2(|a.(bxc)|, 1 + a.b + b.c + c.a)`, in
  the same recursive subdivision that produces the directions, so the
  directions themselves never change. Off by default (`HasSolidAngles()`
  false, `SolidAngles()` empty), so `equal`-mode runs pay nothing --
  confirmed by `tessel_solid_angles`.
- `sources.hpp`/`sources.cpp` -- `PhononSource` gained the
  `toa_quadrature_e {TOA_EQUAL, TOA_AREA}` enum and now owns the mode:
  `Set_TOA_Array(TesselSphere*, toa_quadrature_e = TOA_EQUAL)` sets the
  array and mode together, rejecting `TOA_AREA` (leaving the previous
  array, mode, and weights untouched) if the array lacks solid angles;
  `TOAWeight(k)` returns 1.0 in `TOA_EQUAL` mode and
  `SolidAngles()[k] * (N/4π)` in `TOA_AREA` mode -- `quadrature_mode_switch`
  confirms both the rejection behavior and that the weights sum to N
  either way.
- `scatterers.cpp`/`scatterers.hpp` -- `PopulateProbDists` (the scatterer
  tables) and `IntegrateRows`/`IntegrateRateAndDipole`/`IntegratedRate`
  (the reporting integrals) branch on the mode and multiply by
  `TOAWeight(k)` (tables) or an explicit `solidAngles` argument
  (integrals) in a code path separate from the unweighted `TOA_EQUAL`
  arithmetic. `ReportDirs()` now always builds its fixed degree-6 set
  with solid angles, and `ComputeReportStats()` passes them through in
  `area` mode. `PrintAllScatteringStats()`, in the same file, prints the
  third provenance line, `#  TOA quadrature: area`, in the scatterer-dump
  header.
- `events.cpp` -- `ShearDislocation`'s constructor reads
  `Get_TOA_Quadrature()` once and, in `area` mode, multiplies each ray
  type's per-direction squared amplitude by `TOAWeight(ctr)` before
  `SetRelativeProb()`; see Event source, above.
- `aniso_axi.cpp`/`aniso_axi.hpp` -- unchanged: `AxiProposal`'s per-cell
  sums and in-cell draws read `mPDists`, which already carries the
  weighting applied by `PopulateProbDists` above, so they inherit it with
  no code change; its bounds (`BoundRatio()`/`AcceptanceBound()`, a ratio
  of a cell's proposal weight to that cell's `mPDists` sum) stay
  numerically unchanged in either mode. `area_axi_bounds_valid` and
  `area_axi_matches_general` confirm this.
- `GetRandomScatteredRelativePhonon` (`scatterers.cpp`) and the general
  sampler's majorant -- also unchanged code, inheriting the weighted
  tables the same way.
- `cmdline.hpp`, `main.cpp`, `model.hpp`/`model.cpp` -- add the
  `--toa-quadrature=equal|area` option (`OPT_TOAQUAD`) and two of the
  three provenance lines (`TOA Quadrature: area` in the parameter echo,
  and `(Area-weighted take-off-angle quadrature.)` under the take-off-
  angle count; the third is the `scatterers.cpp` bullet's
  `PrintAllScatteringStats()` line above -- see Consequences, below, for
  all three together). `--event-test` (`main.cpp`) always builds its own
  degree-3 `TesselSphere` in `TOA_EQUAL` mode regardless of the
  command-line option, since it displays per-direction pattern values
  that a weighted table would distort.

## Consequences

Area mode is opt-in (`--toa-quadrature=area`); the default `equal` mode
is unchanged, keeping `aaec540`/`319c8dd` byte-identity (global
constraints, both specs). Every output changes in `area` mode, isotropic
regions included -- this is why the change is opt-in rather than a fix
to the default.

The forced-isotropic and seven- vs five-value equivalences documented in
`doc/MANUAL.md` hold within `area` mode as well as `equal`:
`area_iso_forced` repeats `general_iso_forced`'s check under
`AreaQuadratureScope`, and `tests/test_toa_quadrature.sh`'s seven- vs
five-value check reproduces `test_aniso_inputs.sh`'s isotropic 15-value/
21-value Lop Nor layouts under `--toa-quadrature=area` and confirms
matching output.

Provenance appears only in `area` mode, as exactly three lines: `TOA
Quadrature: area` in the parameter echo, `|          (Area-weighted
take-off-angle quadrature.)` under the take-off-angle count, and `#  TOA
quadrature: area` in the scatterer-dump header; `equal` mode prints
nothing new, and `out_mparams.octv` is unchanged in both modes.
`tests/test_toa_quadrature.sh` checks all three lines verbatim, and that
the default (no-flag) run contains zero mentions of "quadrature".

## Validation

Each item below was a target of the original future-fix validation list;
all are now implemented and checked:

- **Weights sum to 4π within 1e-12.** `tessel_solid_angles`
  (`tests/anisotropy_tests.cpp`) checks this at every degree 0-8
  (relative error ≤1e-12 in every case; ≤2.02e-13 at degree 8, the worst
  case recorded in `toa_quadrature_solid_angles.txt`). PASS.
- **Isotropic MFPs match the exact integral within 1e-4.**
  `area_iso_mfp_exact` checks this at (a=0.5, kappa=0.3, degree 5) and
  (a=2, kappa=0.5, degree 6): both points' area-mode P and S relative
  errors are ≤1e-4 (actual: 4.9e-7/5.4e-6 at the first point, 1.4e-5/
  1.7e-5 at the second), while equal mode's error at the second point
  (0.618% P, 0.786% S) exceeds the test's required 3e-3 threshold,
  confirming the two quadratures are actually discriminated there. PASS.
- **The anisotropic table is reproduced.** `area_report_directional_mfps`
  checks the (2,2,0.5) km crust's P0/P90/S0/S90SV/S90SH against the
  design table within 2e-4 relative, at degree 6. PASS.
- **The sampler's rate and distribution checks pass against weighted
  targets.** `area_general_rate_matches_integration` confirms the mean
  distance to a real scattering event still matches
  `1/IntegratedRate(...)` under area weighting; `area_general_
  distribution_matches_target` confirms the sampled direction
  distribution matches the area-weighted target (χ²/dof < 1.5, the same
  threshold the equal-mode tests use); `area_axi_matches_general` and
  `area_axi_bounds_valid` repeat the axisymmetric sampler's equivalence
  and rigorous-bound checks under area weighting. PASS (all four).
- **The P frame residual is checked.** `area_p_frame_residual`: see
  Frame orientation for P incidence, above. PASS.
- **Deterministic table checks confirm each entry is w_k times its
  equal-mode entry.** `area_tables_weighted` (scatterer tables, within
  1e-11 of each table's total) and `area_event_source_tables`
  (event-source tables, same tolerance) confirm this entry-by-entry;
  `quadrature_mode_switch` confirms `TOAWeight(k)` itself equals
  `SolidAngles()[k]*(N/4π)` to 1e-14 and sums to N to 1e-12. PASS.
- **Outputs are thread-invariant, and the seven- vs five-value
  equivalence holds in area mode.** `area_tables_thread_invariant` and
  `area_report_thread_invariant` check this in-process;
  `tests/test_toa_quadrature.sh`'s determinism checks (THREADS=1
  reproducibility, and THREADS=1 vs THREADS=12 agreement after
  address-masking, for an isotropic layout and the lopnor_21/
  crustpinch_42 anisotropic layouts) and its seven- vs five-value check
  do so out-of-process, over full command-line runs. PASS.
- **A scenario-level report shows how much outputs move.**
  `tests/compare_toa_quadrature.sh`; see Scenario-level effect, below.

Checks added once the option itself existed, beyond the original list:
`--overridemfp` under area mode (`isolated_area_overridemfp`);
`--toa-quadrature=equal` strictly identical to no flag, the three
provenance lines verbatim, and `--event-test` ignoring the option
(`tests/test_toa_quadrature.sh`); default-mode byte-identity against the
`319c8dd` reference build (`tests/compare_aniso_reference.sh` against
`REF_EXE`, global constraints). PASS.

## Scenario-level effect

A manual script, `tests/compare_toa_quadrature.sh`, compares replicate
runs (`equal` against `area`) of Lop Nor and crust pinch, isotropic and
axisymmetric r=4, at take-off-angle degree 7 (K=4 replicates per mode,
seed base 596, 24 threads). Full output is
`docs/toa-quadrature-runs-2026-10/compare-toa-quadrature.out`; see that
directory's `README.md` for the complete provenance and run parameters.

**Energy table.** For each configuration and station, r_s = mean(station
energy, area)/mean(station energy, equal) across the K=4 replicates,
with SE(r_s) from first-order error propagation on the two means'
standard errors (`tests/compare_toa_quadrature.sh`'s header comment), and
z = (r_s − 1)/SE(r_s). A matching "all stations" row, R, applies the same
ratio to each replicate's whole-scenario total instead of one station.

At the whole-scenario level, none of the four configurations' "all
stations" ratios are resolved at 2 SE: lopnor:iso z=+1.67, lopnor:axi4
z=+1.66, crustpinch:iso z=−1.89, crustpinch:axi4 z=−1.79 (all |z|<2).
With only K=4 replicates, this run cannot distinguish a genuine few-
percent scenario-wide energy shift from replicate noise.

At the individual-station level, most of the 80 station-rows (16+16+24+
24, one set per configuration) are likewise consistent with noise, but
two are not: lopnor:iso station 1 moved from a mean energy of 0.2718 to
0.3266 (r=1.202, z=+5.61), and crustpinch:iso station 23 moved from
1.207e-4 to 5.991e-5 (r=0.496, z=−3.98) -- both low-amplitude stations,
where a few phonons changing direction bin is a large fraction of a
small total. No other station-row, in any of the four configurations,
reaches |z|=3 (the next-largest magnitudes are 2.93 in lopnor:axi4 and
2.75 in crustpinch:axi4). These two individual moves are real at this
replicate count, but they do not establish a scenario-wide pattern: most
stations' z stay under 2, and the "all stations" rows above show no
resolved net shift.

**MFP table.** From the first run of each mode (not replicated, so no
SE): across the isotropic layers, P mean free paths move by −0.33% to
+0.06% and S by −1.33% to −0.17%; across the r=4 axisymmetric layers, P
moves by +0.39% to +1.45% and S by −1.90% to −1.28%
(`compare-toa-quadrature.out`, MFP table). The direction of these
per-region shifts (S always shrinking, P usually growing slightly) is
consistent with the isotropic and anisotropic per-degree findings above.

**Degree-9 build check.** Lop Nor at degree 9, 24 threads, median of 3
runs: `area` mode's build time is 0.29 s slower for an isotropic model
(3.810 s to 4.100 s) and 0.32 s slower for the full 20-region
axisymmetric model (4.470 s to 4.790 s) -- both "a few tenths of a
second," as expected (design spec feasibility spike: "about 0.2 s" for
the solid-angle computation alone). The exact extra storage at degree 9
is 41.9 MB (20·4^9 = 5,242,880 directions × 8 bytes;
`docs/toa-quadrature-runs-2026-10/README.md`); the build check's measured
maximum-RSS deltas (+78.0 MB for the isotropic variant, +28.6 MB for the
20-region axisymmetric one) are noisier than this exact figure, swamped
by several-GB-process allocator noise at only 3 replicates per cell --
see that README's "Fix round 1" section for the full accounting. The
41.9 MB figure, not the measured RSS deltas, is the authoritative storage
cost.

## Provenance

The equal-mode tables above come from an independent re-implementation
of the sampler maths (a numpy prototype, not in this repository). The
397.8 km equal-weight S90SH figure is reproduced in-repo by
`report_directional_mfps` (`tests/anisotropy_tests.cpp`). The
frame-orientation ranges above come from the same investigation's
diagnostic scan over 24 frame angles
(2026-10-17-axisymmetric-sampler-performance, `debug-scratch.patch`'s
`scratch_p_rate_vs_frame`, not in this repository); `p_frame_consistency`
covers the resulting agreement between samplers, not a re-derivation of
the ranges themselves.

The area-mode numbers above come from one of three places, named at each
use: an in-repo test in `tests/anisotropy_tests.cpp`, which during a
normal `tests/test_anisotropy_unit.sh` run also writes
`toa_quadrature_solid_angles.txt` (`tessel_solid_angles`),
`toa_quadrature_report.txt` (`area_report_directional_mfps`),
`toa_quadrature_pframe.txt` (`area_p_frame_residual`), and
`toa_quadrature_iso.txt` (`area_iso_mfp_exact`) to `build/`, from where
they were copied into `docs/toa-quadrature-runs-2026-10/`; the manual
scenario script `tests/compare_toa_quadrature.sh`, captured to
`compare-toa-quadrature.out` in the same directory; or the 2026-10-25
design spec's feasibility spike
(`docs/superpowers/specs/2026-10-25-toa-area-quadrature-design.md`), a
throwaway program outside this repository, cited only for the per-degree
P-frame-residual trend (degrees 5 and 7, alongside the in-repo degree-6
figures above). See
`docs/toa-quadrature-runs-2026-10/README.md` for the full capture
procedure, host details, and an md5 confirmation that this capture's own
fix round (storage reservation and reporting detail only) never moved a
computed solid-angle, report-MFP, P-frame-MFP, or isotropic-exactness
value.
