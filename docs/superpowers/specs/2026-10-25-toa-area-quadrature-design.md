# Area-weighted take-off-angle quadrature (opt-in)

Implements the future-fix plan in `docs/toa-quadrature-bias.md`. Branch
`quadrature`, cut from `anisotropy` at `319c8dd`. Every constraint of the
`2026-10-16` and `2026-10-17` specs still applies to the default mode.

## Request

Implement the findings of `docs/toa-quadrature-bias.md`: an area-weighted
take-off-angle quadrature, applied consistently to the isotropic and
anisotropic paths, the event source and the scatterer reports. It must be
opt-in and validated against the note's list.

## Feasibility check (throwaway spike, 2026-10-25)

A scratch program outside the repository replicated `TesselSphere`'s
subdivision. It reproduced the directions bit for bit at degrees 2–5, then
computed each triangle's solid angle Ω_k.

| Check | Result |
|---|---|
| ΣΩ_k − 4π | 1.1e-12 at degree 7, 7e-12 at degree 9 (≤ 6e-13 relative) |
| Ω_k / mean, degrees 6–9 | 0.9277 to 1.2066, rms 8.64% (the note's figures) |
| Isotropic MFP, area vs exact 2-D integral | a = 0.5: 5e-7 (P), 5e-6 (S) at degree 5; a = 2: 5.7e-5 (P), 6.9e-5 (S) at degree 5, 1.4e-5 / 1.7e-5 at degree 6, ~4e-6 at degree 7 |
| Isotropic MFP, equal vs exact | a = 0.5, κ = 0.3: −0.044% (P), +1.09% (S); a = 2: +0.62%, +0.78%, at every degree |
| (2, 2, 0.5) report values, area weights at degree 6 | P0 1602.68, P90 1203.18, S0 754.67, S90SV 442.83, S90SH 374.71 km: the note's table, and within 8.3e-5 of the exact integral |
| P frame-angle spread, r = 10 crust at 60° | equal 5.28% at every degree; area 0.054% at degree 5, 0.013% at 6, 0.003% at 7 |
| P frame-angle spread, r = 10 mantle at 89° | equal 11.2%; area 0.117%, 0.028%, 0.007% at degrees 5, 6, 7 |
| Double-couple and explosion radiation totals | exact (4/15 P, 2/5 S per unit M0²) under both weightings |
| Solid-angle cost, degree 9 | about 0.2 s (TesselSphere itself takes 0.6 s) |

Two consequences shape the design:

- **Area weights also fix the direction density.** With equal weights the
  sampler picks node k with probability proportional to the pattern value
  alone. Nodes are about 30% denser near the icosahedron's 12 vertices, so
  scattered and emitted directions follow the pattern times a smooth
  icosahedral modulation. That modulation is mean(Ω)/Ω_k, which spans 0.829
  to 1.078. The note currently describes only the integrated-rate bias.
- **The radiation totals are exact under both weightings.** Icosahedral
  symmetry integrates the low-order radiation polynomials exactly, so the
  event source's P/S partition does not move. Only its direction
  distribution changes.

## Decisions

1. **Opt-in option `--toa-quadrature=equal|area`, default `equal`.**
   - Equal mode stays byte-identical to `319c8dd`: traces, `.octv` files
     and normalized stdout, for isotropic and anisotropic runs.
   - Isotropic equal-mode runs stay byte-identical to `aaec540`.
2. **Weights.**
   - Ω_k is the solid angle of the icosahedral subdivision triangle whose
     normalized vertex centroid is direction k. It is computed with the Van
     Oosterom–Strackee formula, 2·atan2(|a·(b×c)|, 1 + a·b + b·c + c·a).
   - The quadrature weight is w_k = Ω_k·N/4π.
   - Σw_k = N up to rounding, so every existing "÷N" stays unchanged.
3. **`S2::TesselSphere` can compute solid angles.**
   - It computes Ω_k on request, in the same recursion that produces the
     directions. The directions do not change.
   - This is off by default, so equal-mode runs pay nothing.
4. **`PhononSource` owns the mode.**
   - `Set_TOA_Array(TesselSphere *, toa_quadrature_e = TOA_EQUAL)` sets the
     array and the mode together.
   - Area mode requires solid angles. A call that fails changes nothing.
5. **Weights are applied once, where the tables are built:**
   `Scatterer::PopulateProbDists` and the `ShearDislocation` tables.
   - Everything downstream inherits them with no code change:
     `PopulateWholeProbs`, `ComputeMFPs`, `ComputeDipoles`,
     `GetRandomScatteredRelativePhonon`, the general sampler (proposals and
     majorant), and `AxiProposal` (cell sums, majorants and in-cell draws).
   - `AxiProposal` stays exact whatever the weights. Its bounds stay valid
     as long as weighted cell sums are mirror-symmetric within the 1e-6
     margin. The tessellation is x- and y-mirror symmetric, so they are;
     a test confirms it.
6. **Reports.**
   - `IntegrateRows` takes optional solid angles.
     `IntegrateRateAndDipole` and `IntegratedRate` take them as an
     optional trailing argument.
   - `ComputeReportStats` passes `ReportDirs()`'s solid angles in area mode.
   - `ReportDirs()` is always built with solid angles. At degree 6 that
     costs about 5 ms once per process.
7. **Provenance appears only in area mode,** as three lines:
   - `TOA Quadrature: area` in the parameter echo;
   - `|          (Area-weighted take-off-angle quadrature.)` under the
     take-off-angle count;
   - `#  TOA quadrature: area` in the scatterer-dump header.

   Equal mode prints nothing new, and `out_mparams.octv` is unchanged in
   both modes.
8. **The canonical P frame stays in both modes,** because it keeps the two
   samplers in agreement. Area weighting turns the P frame-orientation
   dependence into a residual that shrinks about 4× per degree. The
   documentation states it that way, not as "removed".
9. **`--event-test` ignores `--toa-quadrature`.** It displays per-direction
   pattern values, and a weighted value would distort that display.
10. **Validation,** with the thresholds in "Design calls" below:
    - weights sum to 4π;
    - isotropic MFPs match the exact integral;
    - the anisotropic table is reproduced;
    - the sampler χ² and mean-distance tests pass against weighted targets;
    - the P frame residual is checked;
    - deterministic table checks confirm each entry is w_k times its
      equal-mode value;
    - outputs are thread-invariant and the seven- vs five-value equivalence
      holds in area mode;
    - default-mode byte-identity is checked against `319c8dd`.
11. **Scenario report:** a manual script compares replicate runs (equal
    against area) of Lop Nor and crust pinch, isotropic and axisymmetric
    r = 4, at degree 7. It reports per-station energy ratios and per-region
    MFP changes, plus a degree-9 build time and memory check. The raw
    output is kept in `docs/toa-quadrature-runs-2026-10/`.
12. **Documentation:**
    - the manual gets the option entry, an updated Validity-limits paragraph
      and a scatterer-dump note;
    - `docs/toa-quadrature-bias.md` gets the new status, results, the
      direction-density finding and the scenario effect;
    - code comments that describe the old TODO are updated.

## Constraints

- In equal mode the `2026-10-17` Global Constraints hold unchanged.
- In both modes:
  - construction draws no random numbers;
  - construction parallelism never changes a table, report value or output;
  - `--threads 1` and `--threads N` outputs are identical at a fixed seed.
- Area mode changes table contents only. The random-draw structure is the
  same in both modes.
- Memory: area mode adds one `Real` per take-off direction to the global set
  (about 42 MB at degree 9). No scatterer's memory changes.
- Equal mode prints no new stdout lines. Area mode prints exactly the three
  lines in decision 7.

## Non-goals

- Making `area` the default. That stays the user's call after review.
- Other node sets (Fibonacci, HEALPix), higher-order rules, or changing the
  tessellation.
- Re-measuring the existing performance documents.
- Weighting `--event-test`.

## Design calls

- The option name and values, and `equal` as the default.
- Thresholds:
  - |ΣΩ − 4π| ≤ 1e-12·4π at degrees 0–8;
  - Ω_k/mean within ±1e-4 of 0.9277 and 1.2066 at degrees 6–8;
  - isotropic MFP within 1e-4 of the exact integral (a = 0.5 at degree 5,
    a = 2 at degree 6), with the equal-weight error required to exceed 3e-3
    for a = 2 so the test discriminates;
  - the (2, 2, 0.5) table within 2e-4 relative;
  - χ²/dof < 1.5;
  - area P frame spread < 1e-3 at degree 6, against > 4e-2 for equal
    weights;
  - table entries equal to w_k times the equal-mode entry within 1e-9
    relative;
  - the event-source P fraction 0.4 within 1e-12 in both modes.
- Weights are computed on the fly as Ω_k·(N/4π). There is no second
  per-direction array.
- The reference integral is Gauss–Legendre in cos ψ times the midpoint rule
  in φ, in the ray frame, with the same integrand as `IntegrateRows`.
