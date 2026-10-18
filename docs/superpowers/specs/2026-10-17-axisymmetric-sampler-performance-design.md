# Axisymmetric sampler performance and near-axisymmetric validation

Follow-up to `2026-10-16-anisotropic-scale-lengths-design.md`, whose
constraints all still apply. Branch `anisotropy`, starting at `3db253a`.

## Requests

1. Make the axisymmetric capability at least as efficient as the triaxial
   (general-sampler) capability.
2. Test the theory by running triaxial mode with `ax ≈ ay` and comparing
   results and performance with the axisymmetric case.
3. Document the take-off-angle quadrature bias (ruling R9). The bias is
   accepted, but it must be documented for users and as a note for a
   possible future fix.

## Measured starting point

Measurements are from `docs/anisotropy-performance.md` and the follow-up
investigation. The workload is Lop Nor at degree 9, 24 threads, Apple M2
Ultra.

| Configuration | Total | Zero-phonon build |
|---|---:|---:|
| Isotropic | 23.0 s | 21.9 s |
| Triaxial r=2 (general sampler) | 27.7 s | 26.5 s |
| Axisymmetric r=2 | 109.7 s | 108.1 s |

Construction dominates the axisymmetric runs:

- **AxiProposal bound loop: ~3.7 s per scatterer.** Each scatterer makes
  73.7M `ComputeBound` calls (81,920 representatives × 2 conversions × 450
  bins). Each call costs about 50 ns and makes five libm calls.
  Everything except one `acos` and two `cos` calls depends only on
  (representative, k_in, k_out), yet it is recomputed in every bin.
- **Geometry: ~0.1 s per scatterer.** This covers the representative
  radii, the tessellations and the mirror match, and it is identical for
  every scatterer at a given take-off-angle degree.
- **Report statistics: ~0.43 s per axisymmetric scatterer (17 states) and
  ~0.23 s per general scatterer (9 states).** `GSATO` is recomputed for
  every state, though it depends only on the report direction.
- **Construction is serial and draws no random numbers.**

Lop Nor has 21 distinct parameter sets (20 axisymmetric) and crust pinch
has 7 (5 axisymmetric). Exact duplicates are already shared through
`GetScattererMatchingParams`, so caching across scatterers gains nothing.

The axisymmetric and triaxial benchmark rows describe **different media**.
At r=10, M is 100 for the axisymmetric medium against 31.6 for the triaxial
one, and the axisymmetric medium has 1.3–1.4× more real events.

## Decisions

1. **Restructure the bound loop without changing results.**
   - Per (representative, k_in, k_out), compute once: m̂, |m_rep|, ρ, the
     m² limits, and the numerator `1 + a_min²·m2_hi`.
   - Per (bin, conversion, cell), track the largest ratio and apply
     `min(M, M·pow(maxRatio, κ+1.5))` once.
   - With the same expressions this is bit-identical to `3db253a`,
     because `pow` and multiplication by M are monotone.
   - Cache the per-degree geometry (representative directions and radii,
     mirror permutations) in a static, keyed by take-off-angle count.
2. **Build in parallel.**
   - Build the 450 bins with a small `ParallelFor` on `std::thread`, using
     the run's thread count: `--threads`, or the automatic worker limit
     when it is 0.
   - Each bin is written by exactly one task, and the majorants are
     reduced in a fixed order after the join. The tables are therefore
     bit-identical for any thread count.
   - Report statistics compute `GSATO` once per report direction and
     integrate the states in parallel. Each state keeps its serial
     summation order, so the report values are also bit-identical.
3. **Trig-free cap kernel.**
   - With `c = m̂·u_c`, use cos(α_c ∓ δ) = c·cosδ ± √(1−c²)·sinδ, clamped
     at 0 and π through `c ≥ cosδ` and `c ≤ −cosδ`.
   - The maths is the same. Results differ from the reference by rounding
     only, within 1e-12 relative, which the 1e-6 bound margin covers.
4. **Faster axisymmetric sampling path.** The random-number call sequence
   and the bound-violation throw are unchanged.
   - Locate the incidence bin once per tentative collision (today `Fold`
     runs twice).
   - Locate bins by comparison against precomputed edge cosines and sines,
     with no `acos` or `atan2`.
   - Evaluate the acceptance ratio in the ray frame with a new
     `ScatterParams::AxiRatio(m2, mu2)`, where `mu = k_out·(t·u_rel) −
     k_in·u_rel.z` and `m2 = k_out² + k_in² − 2·k_out·k_in·t.z`. This
     drops the global `Express` and `ToFabric` transforms.
   - For P, `u_rel = (−√(1−cz²), 0, cz)` with `cz = up·d`. The canonical
     polarization is computed only for an accepted event.
   - Anisotropic outputs change at rounding level, so statistical tests
     replace byte comparisons from this step on.
5. **Same-medium performance gate.**
   - The success criterion is the axisymmetric sampler (`ax = ay`) against
     the general sampler on the same medium (`ay = ax·(1+δ)`, δ = 1e-5).
     Median-of-3 total time must satisfy axi ≤ 1.05 × general in every
     benchmark configuration.
   - The existing 2×-versus-isotropic gate is still reported.
   - If a configuration fails after decisions 1–4, the controller rules on
     a fallback, such as AUTO choosing GENERAL for axisymmetric specs below
     some M.
6. **Validation of the theory with ax ≈ ay.**
   - **Unit level, rates:** IntegratedRate with `ay = ax(1+δ)` stays
     within 5δ of `ay = ax` and scales linearly in δ.
   - **Unit level, samplers:** the axisymmetric and general samplers agree
     on mean distance to a real event and on the (conversion, cell)
     distribution.
   - **Unit level, reports:** report MFPs map AXI↔ANI within 1e-4.
   - **Unit level, cost:** tentative collisions and draws per real event
     are recorded for both samplers.
   - **Scenario level:** a manual script compares groups of replicate runs
     (`ax = ay` against `ay = ax(1+δ)`, plus a null control) with a
     permutation test on windowed trace energies.
   - **Seeds:** two seeds whose XOR is smaller than the phonon count reuse
     the same per-phonon streams (`MixRandomSeed(seed ^ phonon_id)`), so
     replicates use `596 + k·2^32`. This rule goes in the manual's
     `--seed` entry.
7. **Quadrature-bias documentation.**
   - Add a future-fix note, `docs/toa-quadrature-bias.md`, linked from the
     manual, the `ComputeMFPs` TODO and the test comment.
   - Correct the manual's claim that the equal-weight average is bias-free
     for isotropic regions. It gives a direction-independent offset, e.g.
     about +1.1% for Lop Nor crust S.

## Constraints

- Every Global Constraint of the 2026-10-16 plan still holds:
  - isotropic outputs stay byte-identical to `aaec540`;
  - random numbers come only from `RandomUniform()` and `RandomInteger()`;
  - bounds stay rigorous upper bounds, and a violated bound throws;
  - each axisymmetric scatterer stores at most about 9.2 MB;
  - the README, manual and handoff benchmark numbers stay unchanged;
  - the Lop Nor thread-equivalence scripts stay unchanged.
- Model construction draws no random numbers. Construction-time
  parallelism must never change a table, a report value or any output,
  whatever the thread count.
- Anisotropic outputs stay identical between `--threads 1` and
  `--threads N` at a fixed seed.
- Temporary build memory stays at about 32 MB or less per scatterer.
- Normal runs print no new stdout lines.

## Non-goals

- Area-weighted quadrature (documented only).
- Faster isotropic construction, e.g. caching `GSATO`'s `tgamma`
  prefactor. Noted as a separate follow-up.
- A command-line sampler override.
- Changing resolution (cells, bins, representatives).
- Lazy bin construction.

## Design calls

- δ = 1e-5 for the near-axisymmetric general-sampler comparison.
- Thresholds:
  - 1e-12 relative between the kernels;
  - 5δ for the rate check;
  - permutation p ≥ 0.01;
  - the 1.05 same-medium time gate.
- The parallel build uses the simulation's thread setting rather than a
  new option.
