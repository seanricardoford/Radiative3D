# Anisotropic Heterogeneity Scale Lengths Design

## Purpose

Radiative3D describes small-scale heterogeneity with one von Kármán scale
length `a` per region, so scatterers are statistically isotropic. Sanborn et
al. (2017, Appendix A2.1) list anisotropy of heterogeneity scale lengths, after
Margerin (2006) and Jing et al. (2014), as a known extension that was never
implemented. Stretched or flattened heterogeneity changes coda shapes and
amplitude coherence depending on the angle between the wavefront and the
stretch axes. Horizontally stretched heterogeneity, for example, has been
needed to model Pn and Sn codas.

This work lets each region carry three scale lengths `ax, ay, az` and simulates
scattering in the resulting anisotropic random medium. Three modes are named:

| Mode | Scale lengths |
|---|---|
| Isotropic | `ax = ay = az` (today's behavior) |
| Axisymmetric | `ax = ay ≠ az` |
| Anisotropic | any other combination |

## References

- Sanborn et al. 2017, `docs/references/Sanborn17.pdf`, Appendix A2.1.
- Margerin 2006, `docs/references/Margerin06.pdf`. Scalar radiative transfer
  in anisotropic random media. Eqs. 3–6 generalize the spectrum, Eq. 14 gives
  the direction-dependent mean free path, Eq. 30 a cross-section that depends
  on both incoming and outgoing directions, and Sec. 6 a Monte Carlo method
  with direction-dependent free paths and rejection sampling.
- Sato 2019, `docs/references/Sato19.pdf`, Eq. 3: the anisotropic von Kármán
  power spectrum used here. Sato19 cites Jing et al. (2014) for "anisotropic
  scattering", which there means a non-isotropic scattering pattern, not
  anisotropic scale lengths.
- Sato 1995 (`Sato95.pdf`), Jing et al. 2014 (`Jing14.pdf`), and Zeng 2017
  (`Zeng17.pdf`), all in `docs/references/`. Energy-transport theory for
  multiple non-isotropic scattering in statistically isotropic media. See
  "Consistency with the literature" below.

## Decided with the user

- Three modes as tabled above, inferred from the values given. No mode flag.
- The scale-length axes are fixed: z is vertical, and the fabric can't be
  tilted or rotated. Dipping deterministic structure carries the same
  vertically aligned fabric as everything else; the fabric is imprinted on
  the structure, not aligned with it.
- "Vertical" means the local vertical. In curved Earth mappings the fabric
  follows Earth curvature rather than the internal `z` axis.
- Keep the Sato & Fehler elastic Born radiation-pattern terms unchanged and
  swap only the spectrum, evaluated at the full 3-D exchange wavevector.
  Direction-dependent mean free paths, direction-dependent P↔S conversion, and
  polarization-dependent S-wave scattering rates follow from that and are
  intended.
- Sample the way Margerin does: draw from the existing isotropic tables and
  accept or reject with an exact correction. No per-direction pattern tables.
- For axisymmetric media, tailor the proposal to incidence angle. This follows
  the user's observation that the isotropic table, being defined relative to
  the incoming ray, already serves every incidence angle.
- Isotropic regions keep today's code path.
- Acceptance is based on the `do-lopnor.sh` and `do-crustpinch.sh` scenarios,
  non-image outputs only. Octave figures don't need to match.
- Input is new seven-value-per-region `--model-args` layouts for Lop Nor
  (baseline and Moho variant) and crust pinch. Other compiled models stay
  isotropic-only.

## Goals

1. Each region's heterogeneity spectrum holds `ax, ay, az`; the mode is
   inferred from the values.
2. Anisotropic regions scatter according to the physics model below, sampled
   exactly up to the existing take-off-angle discretization. No coarse
   tabulation or interpolation of rates or patterns.
3. Isotropic regions behave exactly as today: same arithmetic and same
   random-number draw order. Fixed-seed runs produce byte-identical
   non-image outputs.
4. Users enter anisotropic values through new `--model-args` layouts; every
   existing layout keeps working.
5. Anisotropic runs repeat exactly under a fixed seed and give identical
   traces across thread counts, matching the guarantees of the multithread
   work.

## Non-goals

- User-specified tilt or rotation of the fabric.
- Background velocity anisotropy. Margerin shows the phase-velocity effect of
  scale-length anisotropy is below 0.5% at 5% fluctuations, and the code
  already neglects velocity shifts.
- Earth-flattening of scale lengths. `EarthCoords::FlattenHetSpec()` stays a
  no-op, as it is for isotropic models today. The flattening transform
  stretches vertical and horizontal lengths by the same factor, so it
  introduces no anisotropy.
- Anisotropic input layouts for compiled models other than Lop Nor (baseline,
  Moho) and crust pinch. Upthrust, Halfspace, SphereEarth, ToySphere, Lop Nor
  Moho2, and ScatParamsStudy stay isotropic-only, though hand-coded C++ models
  can still build anisotropic specs.
- Validation in spherical-shell models.
- Grid-file input, which the code doesn't implement for any parameter.
- Matching Octave or image outputs.
- Performance work on the anisotropic sampler beyond the method below, unless
  measurement shows a problem.

## Anti-goals

These would make the work a failure even if anisotropic runs work:

- Any change to isotropic results, including last-bit floating-point
  differences. Past studies must stay reproducible.
- Slower isotropic runs.
- An approximation that silently biases anisotropic rates or patterns, such
  as nearest-neighbor lookup of a tabulated mean free path.
- Breaking existing do-scripts or the scripts that parse Radiative3D's stdout
  by column position.

## Constraints

- C++11, the existing Makefile, and no new dependencies.
- Scatterer objects are shared read-only across simulation workers. The
  anisotropic path must not add shared mutable state, and it must draw random
  numbers only from the per-phonon stream (`RandomUniform()`).
- Scatterer tables are already large, so per-incoming-direction pattern
  tables are not feasible. The sampling method below needs no tables beyond
  the ones an isotropic scatterer already has.
  - The binary's default take-off-angle degree is 7 (327,680 directions).
    Waveform runs through the do-scripts use degree 9 (5,242,880 directions),
    set by `PopDefaults` in `scripts/do-fundamentals.sh`.
  - Lop Nor builds about 21 scatterers, one per layer. At degree 9, model
    build already takes about 22 s and 4.7 GB today. Isotropic runs must not
    add to either.
  - Axisymmetric scatterers may add about 4.4 MB each, about 2% for Lop Nor,
    plus their bound-table build. The bound tables take roughly 1–2 s per
    scatterer single-threaded. If that is material for a model with many
    axisymmetric layers, compute them in parallel across scatterers.

## Physics model

### Fabric frame

`ax, ay, az` are correlation lengths along local East, North, and Up at the
scattering point, as returned by the global `ECS` object (`GetEast`,
`GetNorth`, `GetUp`).

- In flat and flattened mappings (`ENU_ORTHO`, `RAE_ORTHO`), these are the
  internal `+x, +y, +z` axes everywhere. Lop Nor uses this case.
- In curved mappings (`RAE_CURVED`, `RAE_SPHERICAL`), the frame follows the
  local vertical, so horizontally stretched heterogeneity stays horizontal
  across the model. Crust pinch uses `RAE_CURVED`. At its 950 km range the
  local vertical is tilted about 8.5° from internal `z`.

### Spectrum

Sato19 Eq. 3, with the exchange wavevector `m` expressed in the fabric frame
as `(mE, mN, mU)`:

```
P(m) = 8 π^(3/2) ε² ax ay az Γ(κ+3/2) / [ Γ(κ) (1 + ax² mE² + ay² mN² + az² mU²)^(κ+3/2) ]
```

With `ax = ay = az = a` this is exactly today's `ScatterParams::PSATO()`.

### Scattering coefficients

Today's Sato & Fehler (4.50–4.52) terms are kept, with `l = ω/Vs` (`el` in the
code), `γ0 = Vp/Vs`, wavenumbers `kP = l/γ0` and `kS = l`, and `(ψ, ζ)` the
scattered direction relative to the incoming ray and its polarization:

```
gPP = (l⁴/4π) XPP²(ψ,ζ) P(m)
gPS = (1/γ0) (l⁴/4π) XPS²(ψ,ζ) P(m)
gSP = γ0 (l⁴/4π) XSP²(ψ,ζ) P(m)
gSS = (l⁴/4π) [XSSψ²(ψ,ζ) + XSSζ²(ψ,ζ)] P(m)
```

The spectrum is evaluated at the exchange wavevector

```
m = k_out ê_out − k_in ê_in
```

where `ê_in` and `ê_out` are the global incoming and scattered directions and
the `k` values follow the in/out ray types. In the isotropic case `|m|`
reduces to the arguments `PSATO()` receives today.

What follows from this:

- The total scattering rate `g0 = (1/4π) ∮ Σ g dΩ` depends on the incoming
  direction. For S waves it also depends on polarization, because `XSP` and
  `XSS` depend on `ζ`, which is measured from the polarization direction.
- P→S and S→P fractions depend on direction.
- The SV/SH split of S→S scattered energy (`spol`) is still given by the
  relative-frame `X` terms, because the spectrum scales both components by the
  same factor.

### Consistency with the literature

Sato (1995) formulates multiple non-isotropic scattering in an energy-transport
equation. Jing et al. (2014) and Zeng (2017) apply it to S-coda envelopes from
the 2008 Wells, Nevada, sequence. All three assume a statistically isotropic
medium, so the scattering coefficient depends only on the scattering angle ψ,
which is axially symmetric about the incident ray. None of them treats
anisotropic scale lengths.

These papers support the design in three ways:

- **Same formulas as the code.** Jing14's S→S Born coefficient (Eqs. 19–21,
  with μ in place of ν) is term for term the code's `XSATO()`/`GSATO()` S→S
  term. Its normalization matches as well: the total coefficient is
  `g0 = (1/4π) ∮ g dΩ` and the mean free path is `1/g0`, which is what
  `Scatterer::ComputeMFPs()` computes.
- **The Monte Carlo already solves the elastic version.** Radiative3D is a
  Monte Carlo solution of the elastic generalization of Sato95's transport
  equation, with P↔S conversions; Zeng17 and Jing14 consider S→S only. Each
  phonon scatters at rate `g0` into directions in proportion to `g`, and
  energy is conserved because scattering never removes phonons. The
  anisotropic sampler keeps all of this. Null collisions change nothing about
  the phonon.
- **Anisotropy breaks one of their assumptions, as expected.** With
  anisotropic scale lengths the coefficient depends on both the incoming and
  outgoing directions, not just ψ, and `g0` depends on the incoming direction.
  That breaks the axial-symmetry assumption Sato95's spherical-harmonic
  solution relies on. Margerin (2006, Eqs. 29–30) gives the generalized
  transport equation, `σ(p̂, p̂′)` with a direction-dependent mean free path,
  and the sampler draws from it exactly. Evaluating the anisotropic spectrum
  at the exchange wavevector inside the elastic coefficients is valid because
  the `X` terms come from the form of the elastic perturbation, not from the
  medium's statistics.

The papers' validity limits carry over and go in the user documentation:

- Born and Bourret validity requires `k² a_max² ε² ≪ 1`, using the largest
  scale length. For Lop Nor values (`k ≈ 3.6 /km`, `a_max = 5 km`,
  `ε = 0.01`) it is 0.03.
- Margerin cautions that for extremely flat scatterers (`az → 0`,
  `ah → ∞`) the medium acts as fine layering, where transport theory breaks
  down.

## Sampling method

The method is delta tracking (null collisions), with the scatterer's own
isotropic tables as the proposal. It generalizes Margerin's rejection approach
to cover the free path as well as the direction. Triaxial scatterers use this
general sampler as described here. Axisymmetric scatterers use the same
framework with a proposal tailored to incidence angle; see "Axisymmetric
sampler" below.

For the general sampler:

- `a_min = min(ax, ay, az)` and `M = ax ay az / a_min³`.
- Reference tables are exactly the tables an isotropic scatterer with
  `a = a_min` builds, with the same ν, ε, κ, velocities, and frequency. They
  give the proposal: conversion type from the whole-probability arrays, then a
  take-off angle from the chosen distribution. This is today's
  `GetRandomScatteredRelativePhonon()` logic.
- The majorant mean free path per incoming ray type is the isotropic mean free
  path at `a_min` divided by `M`.
- Tentative collision distances along the path are exponential with the
  majorant mean free path.
- At a tentative collision at distance `s`:
  1. Take the phonon's position and direction at `s`. Keep the polarization
     angle unchanged, as `Phonon::Move()` already does along curved rays.
  2. Draw a proposal (conversion type and take-off angle).
  3. Compute `ê_out` with the same `OrthoAxes` construction
     `Phonon::Transform()` uses.
  4. Form `m` in the fabric frame and accept with probability

     ```
     p = [ (1 + a_min² |m|²) / (1 + ax² mE² + ay² mN² + az² mU²) ]^(κ+3/2)
     ```

     This is never above 1, because `a_min` is the smallest scale length.
- If accepted, the event is a real scattering event: move the phonon to `s`,
  apply the proposal's relative phonon through `Transform()`, and report the
  scatter event, all as today.
- If rejected, it's a null collision. The phonon isn't moved, nothing is
  reported, and the move count doesn't change. Draw the next tentative
  distance from `s`. If the accumulated distance reaches the cell boundary,
  handle the boundary exactly as today.

Why this is exact: the rate of accepted events into conversion type `t` and
solid angle `dΩ` is the majorant rate times the proposal probability times
`p`, which equals `g_t(ê_out | incoming ray)/4π` for the anisotropic medium.
Direction- and polarization-dependent rates, conversion fractions, and
patterns therefore all come out right without new tables. Because the
majorant bounds every direction, this also holds along curved rays, where the
direction changes within a segment.

An anisotropic scatterer under the general sampler uses the same memory as an
isotropic one.

Isotropic scatterers never enter this path. The propagation loop keeps
today's call and random-draw sequence for them.

### Mode detection

A spec is isotropic if and only if `ax == ay == az` exactly. It is
axisymmetric if and only if `ax == ay ≠ az` exactly, and those specs get the
axisymmetric sampler. All other specs are anisotropic and get the general
sampler.

### Axisymmetric case

With `ax = ay = ah`, the spectrum depends only on the horizontal and vertical
parts of `m`:

```
P(m) = 8 π^(3/2) ε² ah² az Γ(κ+3/2) / [ Γ(κ) (1 + ah² mH² + az² mU²)^(κ+3/2) ]
```

Here `mH² = mE² + mN²`. Equivalently, it is an isotropic spectrum with a
direction-dependent scale length `a(θm)² = ah² sin² θm + az² cos² θm`, where
`θm` is the angle of `m` from vertical (Margerin 2006, Eq. 4).

It follows that, in the frame of the incoming ray, the scattering pattern
depends only on where the local vertical sits in that frame:

- For P, only the incidence angle `θ` from vertical matters.
- For S, `θ` and the polarization angle `χ` from the vertical plane matter;
  `χ` runs from SV (0°) to SH (90°).

The mean free path follows the same dependence. Axisymmetry is used three
ways:

- **Reporting.** Axisymmetric rows in the scatterer dump report mean free
  path against incidence angle (see "Output and reporting") instead of
  duplicated East and North values.
- **Validation.** Mean-free-path curves against incidence angle are cheap to
  compute and check. See check 4e.
- **Sampling.** See the next section.

### Axisymmetric sampler

**Why the general sampler is slow here.** The isotropic table is built in the
frame of the incoming ray. An isotropic medium looks the same from every
direction, so that one table is exact for every incidence angle: in effect,
it is already tailored to incidence. An axisymmetric medium is not the same
from every direction. Its pattern in the incoming-ray frame changes with `θ`
(and with `χ` for S). The isotropic proposal can't follow that change, so
proposals get rejected in proportion to the mismatch.

Measured at degree 9, the cost per draw is set by the table lookup, not the
physics:

| Step | Cost |
|---|---:|
| Drawing a proposal (binary search over a 5.24M-entry table) | ~156 ns |
| Evaluating the spectrum and the acceptance test | ~14 ns |

So the number of draws per real event is what matters. The isotropic
proposal's mismatch is fundamental: the bound `M` is reached along the whole
short-axis direction of `m`, not just at `m = 0`.

**The axisymmetric sampler.** It uses incidence-tailored tables, as the
isotropic case does in effect, factored to fit in memory. Full-resolution
tables per incidence bin are out of the question: each would cost as much as
the whole isotropic table set, about 200 MB per scatterer at degree 9. Instead
the proposal is the product of two pieces:

- the existing full-resolution isotropic table at `a_min`, shared as today;
- a small, coarse, per-incidence-bin weight table that tilts proposals toward
  the directions the axisymmetric medium favors at that incidence.

The exact accept/reject step then removes all remaining mismatch, so sampling
stays exact at full angular resolution.

Details:

- **Coarse cells.** The take-off-angle set is generated by recursive
  subdivision of the icosahedron and stored in hierarchical order (see
  `TesselSphere` in `geom_s2.cpp`). Each degree-3 triangle (1,280 of them)
  therefore owns a contiguous block of `4^(D−3)` fine directions, where `D` is
  the run's degree. Below degree 3, each fine direction is its own cell.
  - The sum of the isotropic table over a cell is the difference of two
    entries of the existing cumulative array.
  - A draw within a cell is a binary search over that cell's range.
- **Incidence bins.** `θ` every 2° and, for S, `χ` every 10°.
  - Symmetries fold the bins to `θ, χ ∈ [0°, 90°]`.
  - The ratio `R = P_aniso(m)/P_iso(|m|; a_min)` depends on `(m·Û)²`, so it is
    unchanged when the local vertical `Û` is replaced by `−Û`.
  - The tessellation and the isotropic P and S patterns are symmetric under
    reflection through the two coordinate planes of the incoming-ray frame
    that contain the ray. These reflections map coarse cells to coarse cells
    and are applied as precomputed cell permutations.
- **Bin tables.** For each bin, conversion type `c`, and cell `j`:
  - Store an upper bound `B(c, j)` on `R`, valid over every fine direction in
    the cell and every incidence state in the bin, capped at `M`.
  - Store the cumulative proposal weights `B(c, j) · S(c, j)`, accumulated
    over both conversion types of the incoming ray type, where `S(c, j)` is
    the cell's isotropic-table sum.
  - The cell-choice probabilities and the bound used in the acceptance test
    must come from the same stored numbers, so storage rounding can't bias
    results.
  - Memory is about 4.4 MB per axisymmetric scatterer for one `float` value
    per entry: 0.4 MB for P and 4.0 MB for S. That's about 2% on top of
    today's Lop Nor tables at degree 9, or up to twice that depending on
    precision.
- **Rigorous bounds.**
  - Compute the bounds at model build as an exact maximum over mid-level
    representative directions, using the run's own directions if its degree
    is 6 or less and degree-6 representatives otherwise.
  - Add an analytic slack that covers each representative's descendants: a
    descendant's `m` lies within `k_out` times the chord to the
    representative. The slack also covers the bin's spread of vertical
    directions, a cap with the bin's half-diagonal as radius.
  - A relative margin of about 1e-9 covers rounding in the symmetry folds.
  - Build cost is about 75M bound evaluations, roughly 1 s per axisymmetric
    scatterer, whatever the run's degree.
- **P-incidence frame.** P has no physical polarization, so the sampler sets
  the P phonon's polarization angle so that the frame's first axis lies along
  the local `θ̂`, the downward direction in the vertical plane. The local
  vertical then sits at a fixed azimuth in the frame, and P needs only `θ`
  bins. Scattered P→S polarization lies in the scattering plane whatever the
  frame, so this changes no physics.
- **Rates.** Each bin's majorant is `g*_bin = (1/nTOA) Σ B·S` over its
  entries, which bounds the true rate for every incidence state in the bin.
  Tentative collisions use `g*_max`, the largest bin majorant for the
  incoming ray type. At a tentative collision:
  1. Locate the bin.
  2. Thin with probability `g*_bin/g*_max`. This costs no table draw.
  3. Draw a cell from the bin's coarse cumulative weights, then a fine
     direction within the cell.
  4. Accept with probability `R/B`.

  This is exact on straight and curved rays.
- **Fail loudly.** `R/B` is computed for every proposal. A value above
  `1 + margin` means a bound is wrong; the run aborts with an error rather
  than clamping and silently biasing results.

Expected draws per real event, from the prototype (fine degree 7, Lop Nor-like
crust, κ 0.5; P at 0°–89° incidence, S at 30° and 89° with χ = 0°, 45°, 90°).
In every case tested, the bounds were confirmed to be at or above the exact
ratio at every fine direction.

| `ah, az` (km) | General sampler | Axisymmetric sampler |
|---|---:|---:|
| 2, 0.5 | 4.5–9.1 | 1.5–1.9 |
| 5, 0.5 | 13–53 | 2.7–5.1 |
| 0.5, 2 (vertical cigars) | 2.1–3.1 | 1.2–1.6 |

Each draw costs about what a draw costs today. One draw from the bin's coarse
table, which is small and cache-friendly, replaces the top of today's binary
search.

The same idea doesn't extend to triaxial media within memory. Their P
incidence family is two-dimensional and their S family three-dimensional,
which would need hundreds of megabytes per scatterer for S. Triaxial media
keep the general sampler.

### Interaction with existing overrides

- `--overridemfp`: anisotropic scatterers draw path lengths from the override
  mean free path directly, with no null collisions. At each event they redraw
  proposals until one is accepted, so the deflection pattern still depends on
  the incoming direction and polarization.
- `--nodeflect`: delta tracking runs as normal, so event density reflects the
  anisotropic rates, but accepted events don't deflect.
- Both together: the existing collapse to a single dummy isotropic scatterer
  stays.

### Scatterer matching

`Scatterer::GetScattererMatchingParams()` and
`ScatterParams::CompareRoughly()` must distinguish all three scale lengths. An
anisotropic scatterer must never merge with an isotropic one, or with another
anisotropic one unless all parameters are identical. Today's exact-equality
matching for isotropic specs must stay unchanged.

## Input

- `Elastic::HetSpec` holds `ax, ay, az`; the existing isotropic constructors
  and schema set all three to `a`. Add a schema with argument order
  `(nu, eps, ax, ay, az, kappa)`.
- Scale lengths must be positive. Otherwise model construction fails with a
  clear message.
- New `--model-args` layouts use seven values per region:
  `nu, eps, ax, ay, az, kappa, Q`. As in today's five-value layouts, `Q` is
  `Q_s` with `Q_kappa` infinite.

| Model | New count | Regions, in order | Mirrors |
|---|---:|---|---:|
| Lop Nor baseline (`LopNorCylinder`) | 21 | sediments, crust, mantle | 15 |
| Lop Nor Moho (`LopNorCylinderMoho`) | 28 | sediments, crust, mantle, Moho transition | 20 |
| Crust pinch (`CrustPinchWCG`) | 35 | sediments, crust, pinched crust, Moho, mantle | 25 |
| Crust pinch | 38 | as above, then `SediThick, CrustThick, MohoThick` | 28 |
| Crust pinch | 42 | as above, then `PinchFrac, DepthFrac, SediFrac, MohoFrac` | 32 |

- Model-index 1 dispatch in `user.cpp` selects the Moho model when the count is
  20 or 28 and the baseline otherwise. Today it sends any count of 20 or more
  to the Moho model, which would misroute 21.
- Existing counts keep working unchanged, and unknown counts still fail with
  the existing error.
- Example axisymmetric crust, stretched 4:1 horizontally:
  `SCAT2=0.8,0.01,2.0,2.0,0.5,0.3,1000`.
- `do-lopnor.sh` and `do-crustpinch.sh` gain comments documenting the
  seven-value form, with no functional change. `doc/MANUAL.md` documents the
  layouts, the modes, and the fabric frame.

## Output and reporting

Runs whose regions are all isotropic must produce today's stdout exactly,
apart from addresses, the revision string, and timings. Downstream scripts
parse stdout by column position:

- The do-scripts feed the `--dump-grid` table to Octave `load`, which reads
  columns 1–15.
- `do-halfspace-nearsrc50.sh` reads column 11 of the grid dump and computes
  column 5 × column 9 (`el × MFP(S)`) from the scatterer dump.

So anisotropic output only appends columns and never moves existing ones:

- **Grid dump.** If any node carries an anisotropic spec, every node row,
  including the `***` rows, gets two numeric columns, `ay` and `az`, after
  `kappa`, and the existing `a` column holds `ax`. Otherwise the format is
  unchanged.
- **Scatterer dump.** Isotropic rows are unchanged. Anisotropic rows keep
  every existing column in place:
  - The `a` column holds `ax`.
  - `MFP (P)` and `MFP (S)` hold the mean free path for the scattering rate
    averaged over incidence along the three fabric axes. For S, the average
    also covers both transverse polarizations.
  - The `DM` columns are averaged the same way.
  - Appended at the end of each row: `ay`, `az`, a mode label (`AXI` or
    `ANI`), and directional mean free paths:
    - `ANI` rows give P along E, N, and U, plus S along each axis for each of
      the other two polarization axes. For example, S traveling east
      polarized up is SV, and polarized north is SH.
    - `AXI` rows give P at incidence angles of 0°, 45°, and 90° from
      vertical, then S at 0° (horizontally polarized, so SV and SH coincide),
      and SV and SH at 45° and at 90°.
  - These values come from direct numerical integration of the anisotropic
    coefficients at model build. The integration uses a fixed degree-6
    direction set (81,920 directions) whatever `--toa-degree` is, so it costs
    about 0.1 s per anisotropic scatterer. Integrating over the run's own
    degree-9 set would add several seconds per scatterer, about a minute per
    Lop Nor build with an anisotropic mantle. The dump prints on every run,
    so that cost would recur every time.
- Seismogram files, `seis_traces_asc.dat`, `out_mparams.octv`, and report
  formats are unchanged.

## Validation

Test scripts follow `tests/run_lopnor_thread_case.sh`: simulation only, fixed
`--seed`, explicit `--threads`, no Octave.

- Scenario parameters come from the do-scripts, including the 600 s recording
  time, with two reductions: fewer phonons, and take-off-angle degree 5 (as in
  the existing harness) instead of the do-scripts' 9.
- The scripts accept overrides for degree and phonon count, so the same
  comparisons can be run at full do-script fidelity on demand.
- Degree doesn't change which code paths run, so byte-identity at degree 5
  exercises the same isotropic code as at degree 9.
- Every scenario run must show that real scattering occurred in its
  scattering regions, for example from scatter-event report counts.
- The existing thread-equivalence case uses a 0.02 s time-to-live. Its traces
  are all zero, so it neither scatters nor records anything.

The design review confirmed these checks are achievable:

- A 600 s isotropic Lop Nor run with 400,000 phonons gave byte-identical traces
  and `seis_*.octv` files at 1 and 12 threads.
- A 20,000-phonon crust pinch run with the do-script parameters scattered in
  every depth band.

1. **Isotropic regression.** Lop Nor and crust pinch with today's parameters.
   The test builds the baseline revision (the commit this branch starts from,
   `aaec540`) in a temporary git worktree, runs both builds with identical
   arguments, and compares byte for byte:
   - `seis_traces_asc.dat`
   - every `seis_*.octv`
   - `out_mparams.octv`
   - stdout, after masking addresses, the revision string, and timing lines

   A live baseline avoids depending on platform floating-point details, which
   committed golden files would.
2. **Seven-value isotropic equivalence.** Each scenario with seven-value
   layouts where `ax = ay = az = a` must be byte-identical to the five-value
   layout. This exercises parsing and mode detection.
3. **Anisotropic scenarios.** Each scenario with at least one axisymmetric
   region and one triaxial region must:
   - run to completion;
   - repeat byte-identically with the same seed;
   - give identical traces with `--threads 1` and `--threads 12`, as in
     `tests/test_lopnor_thread_equivalence.sh`;
   - list the anisotropic scatterers in the scatterer dump with the given
     scale lengths and mode labels.
4. **Sampler correctness.** For representative anisotropic scatterers, run
   these checks for incoming states covering P along each fabric axis plus one
   oblique direction, and S along each axis with both transverse
   polarizations:
   - a. Mean distance to an accepted event agrees with `1/g0` from direct
     numerical integration over the scatterer's own take-off-angle set (not
     the degree-6 reporting set), within statistical tolerance.
   - b. The distribution of accepted conversion types and take-off angles
     agrees with the directly computed target distribution.
   - c. With `ax = ay = az` forced through the anisotropic path, every
     proposal is accepted (`p = 1` within rounding) and the mean free paths
     equal the isotropic ones. Scenario outputs on the forced path match the
     isotropic path statistically; an exact match is acceptable but not
     required.
   - d. Directional mean free paths have the expected symmetries: invariance
     under reflection of the incidence direction, and for axisymmetric specs
     invariance under rotation about Up.
   - e. For axisymmetric specs, P and both S polarizations: mean free path
     against incidence angle every 5° from 0° to 90°, from direct
     integration. These curves are recorded for review against the
     qualitative trends of Margerin's Figs. 2–4. Margerin's curves are
     scalar, so only the trends are comparable, not the values.
   - f. For axisymmetric scatterers, the axisymmetric sampler and the general
     sampler agree statistically on mean distance to an accepted event and on
     the distribution of accepted conversion types and take-off angles. Both
     target the same distribution, so the general sampler is the reference.
   - g. Bound validity: `B ≥ R` at every fine direction, for a dense set of
     incidence states covering bin edges and every folded quadrant. This
     covers several values of `ah`, `az`, κ, and frequency, including
     vertically elongated media (`az > ah`).
   - h. Draws per real event are recorded for both samplers and compared with
     the prototype values under "Axisymmetric sampler".
5. **Performance.**
   - `tests/benchmark_lopnor_threads.sh` (isotropic) shows no slowdown beyond
     run-to-run noise.
   - Anisotropic overhead for both scenarios is measured and recorded. There
     is no fixed budget.

## Expected cost

These numbers are for the general sampler. The axisymmetric sampler's are
under "Axisymmetric sampler". They come from a throwaway numerical estimate
with Lop Nor-like crust: 2 Hz, Vp 6.2 km/s, Vs 3.6 km/s, ν 0.8, ε 0.01, κ 0.5,
and 200,000 Fibonacci-sphere directions. The range covers vertical, 45°, and
horizontal incidence, and both S polarizations.

| `ax, ay, az` (km) | `M` | Tentative collisions per real event, general sampler |
|---|---:|---:|
| 1, 1, 0.5 | 4 | 2.1–2.7 |
| 2, 2, 0.5 | 16 | 4.2–9.0 |
| 5, 5, 0.5 | 100 | 13–53 |
| 2, 1, 0.5 | 8 | 2.6–4.9 |

Approximate directional mean free paths for `2, 2, 0.5`, a rough cross-check
for the integration code (the code's take-off-angle set will differ
slightly):

| Ray | Mean free path (km) |
|---|---:|
| P, vertical | 1603 |
| P, horizontal | 1203 |
| S, vertical | 755 |
| S, horizontal, SV | 443 |
| S, horizontal, SH | 375 |

At ε = 0.01 the mean free paths run from hundreds to thousands of kilometers,
so a phonon sees few events in 600 s and null-collision overhead is
negligible for these scenarios. Overhead scales with event count in strongly
scattering models.

Each tentative collision costs:

- one exponential draw;
- one proposal draw, a binary search over the take-off-angle table;
- one direction transform;
- one `pow`;
- for curved rays, one `AdvanceLength`.

Throwaway micro-benchmarks put the proposal draw at about 42 ns, 68 ns, and
156 ns at degrees 5, 7, and 9, and the spectrum plus acceptance test at about
14 ns. So a tentative collision at degree 9 costs roughly 0.2 µs. Under the
general sampler that's about 1–2 µs per real event at aspect ratio 4, and
about 3–11 µs at aspect ratio 10. Under the axisymmetric sampler it's about
0.3–0.4 µs and 0.5–1 µs.

**Measurement gate.** Once the anisotropic scenario tests pass, time the
anisotropic versions of both scenarios at production settings:

- degree 9, with phonon counts chosen for a reasonable run time;
- axisymmetric and triaxial regions at aspect ratios 2, 4, and 10;
- the do-scripts' ε and a strongly scattering ε of 0.05.

Record the results next to the matching isotropic runs, with draws per real
event for each sampler. If an anisotropic run is more than 2× slower than its
isotropic match:

- with axisymmetric regions, revisit the axisymmetric sampler's cell and bin
  resolution within its memory budget;
- with triaxial regions, consider extending incidence-tailored proposals to
  triaxial P waves, designed separately.

Otherwise nothing more is needed.

## Design calls (not discussed with the user)

- Delta tracking is the concrete form of the agreed accept/reject sampling.
  There are no mean-free-path tables and no fixed performance budget for
  anisotropic runs.
- The axisymmetric sampler's design and parameters are mine, worked out in a
  prototype: degree-3 cells, 2° × 10° incidence bins, degree-6 representatives
  for the bounds, `float` storage, and aborting on any bound violation.
  Triaxial media keep the general sampler.
- The 2× threshold in the measurement gate is my choice.
- The isotropic regression compares against a live baseline build, not
  committed golden files.
- Anisotropic rows in the output dumps append columns. Their legacy
  mean-free-path columns show axis-averaged values.
- The forced-anisotropic check with equal scale lengths is deterministic
  (every proposal accepted) on top of the agreed statistical comparison.
- Scale lengths must be positive.
- Automated tests run at take-off-angle degree 5 rather than the do-scripts'
  9, with overrides available for full-fidelity runs.
- Reported anisotropic mean free paths are integrated over a fixed degree-6
  direction set, not the run's own set.

## Left to the builder

- Names, data layout, and code placement. Expected touch points:
  `elastic.hpp`, `scatparams.*`, `scatterers.*`, `phonons.cpp`, `user.cpp`,
  `user_LopNorCyl_inc.cpp`, `user_LopNorCylMoho_inc.cpp`, `user_NSCP_inc.cpp`,
  `grid.cpp`, `doc/MANUAL.md`, the two do-scripts' comments, and `tests/`.
- The mechanism for the sampler checks: a hidden diagnostic mode in `main` or
  a separate test executable linked against the same objects.
- Exact order and formatting of appended report columns, within the
  constraints above.
- Stdout masking details for the regression comparison.
- Test region values, phonon counts, and statistical tolerances.
- Micro-optimizations, such as precomputed take-off-angle unit vectors or
  skipping `AdvanceLength` for straight-ray cells, provided exactness and the
  isotropic path are untouched.
- Axisymmetric sampler internals: bin indexing, construction of the
  reflection permutations, computing the P-frame angle in curved mappings,
  whether to parallelize the bound-table build, and retuning cell, bin, or
  representative resolution. Exact bounds and the memory budget must hold.
- Order of work. Building and validating the general sampler first gives the
  axisymmetric sampler its reference (check 4f).
