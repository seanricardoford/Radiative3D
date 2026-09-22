# Axisymmetric Anisotropic Scattering Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace per-event full-sphere anisotropic scattering scans with a cached axisymmetric von Karman kernel while preserving 3D phonon trajectories, polarization, isotropic behavior, and seeded reproducibility.

**Architecture:** `ScatterParams` remains the source of PSD and Sato/Fehler radiation physics. A new `AxisymmetricScatteringKernel` builds immutable reduced quadrature/CDF tables over incoming `mu_in`, outgoing `cos(psi)`, and relative azimuth `zeta`; `Scatterer` uses those tables for directional MFPs and sampling. The sampled meridional direction and polarization are converted into the existing incoming-polarization relative frame before `Phonon::Transform()` applies the full 3D trajectory update.

**Tech Stack:** C++11, GCC/G++, `std::thread`, `std::unique_ptr`, `std::vector`, `std::mt19937_64`, GNU Make, and native assertion-based tests.

**Spec:** `docs/superpowers/specs/2026-09-22-axisymmetric-anisotropic-scattering-design.md`

## Global Constraints

- Keep the C++11 and `-pthread` build.
- Use explicit `RandomEngine&` instances for all simulation randomness.
- Build and freeze all shared quadrature/CDF tables before worker threads start.
- Use the local model vertical/radial unit vector as the anisotropy symmetry axis.
- Preserve the exact legacy path when `ah == av`.
- Preserve `--overridemfp` precedence over computed directional MFPs.
- Do not add user-facing quadrature-resolution options.
- Keep anisotropy global per model invocation; do not introduce per-cell CLI axes.
- Run `make -j2` and `make test` after simulation, scattering, or build changes.
- Run `make test-plotting` after the implementation because the complete recipe workflow must remain valid.

---

## File map

| File | Responsibility in this plan |
| --- | --- |
| `axisymmetric_scattering.hpp` | Reduced-kernel data types and read-only sampling/query interface. |
| `axisymmetric_scattering.cpp` | Quadrature, incoming-angle cache, CDF construction, interpolation, and sampling. |
| `scatparams.hpp/cpp` | Local-vertical meridional basis helper and corrected directional `GSATO` geometry. |
| `sources.hpp/cpp` | Compact `PhononSource` angular-storage constructor used only by anisotropic scatterers. |
| `scatterers.hpp/cpp` | Kernel ownership, anisotropic constructor branch, directional MFP/sampling delegation, frame conversion, and diagnostics. |
| `phonons.cpp` | Pass the incoming polarization angle into the geometry-aware scattering API. |
| `Makefile` | Build the new object and native regression executable. |
| `tests/test_anisotropic_scattering.cpp` | PSD and local-vertical geometry invariance tests. |
| `tests/test_axisymmetric_scattering.cpp` | Reduced quadrature reference, CDF, interpolation, and frame-conversion tests. |
| `tests/test_parallel_reproducibility.sh` | Add a small unequal-anisotropy one-worker/two-worker comparison. |
| `docs/DEVELOPMENT.md` | Document the cache, coordinates, freeze point, diagnostic summaries, and numerical limits. |
| `docs/MANUAL.md` | Record any user-visible clarification of existing anisotropic options. |
| `PLANS.md` | Mark this performance/capability work and remaining limitations. |

## Interfaces introduced by the plan

The following names and signatures are the contract between tasks:

~~~
// axisymmetric_scattering.hpp
class AxisymmetricScatteringKernel {
public:
  enum conversion_e { GPP, GPS, GSP, GSS, NUM_CONVERSIONS };

  struct Sample {
    Real psi;
    Real zeta;
    Real polarization;
    raytype output_type;
  };

  explicit AxisymmetricScatteringKernel(const ScatterParams & params);

  Real GetInverseMeanFreePath(raytype intype, Real mu_in) const;
  Real GetConversionWeight(conversion_e conversion, Real mu_in) const;
  Real GetAverageInverseMeanFreePath(raytype intype) const;
  Real GetAverageConversionWeight(conversion_e conversion) const;
  Real GetAverageDipole(raytype intype) const;
  Sample GetRandomSample(raytype intype, Real mu_in,
                         RandomEngine & rng) const;

  std::size_t IncomingBinCount() const;
  std::size_t DirectionNodeCount() const;
  bool CDFsAreValid(Real tolerance) const;
};

// scatterers.hpp
Phonon GetRandomScatteredRelativePhonon(
    raytype intype, const R3::XYZ & incoming,
    const R3::XYZ & vertical, Real incoming_polarization,
    RandomEngine & rng);
~~~

The existing geometry-aware overload without `incoming_polarization` remains a compatibility wrapper that passes `0.0`.

### Task 1: Correct local-vertical scattering geometry

**Files:**
- Modify: `scatparams.hpp` near the directional `GSATO` declarations.
- Modify: `scatparams.cpp` in the directional `GSATO` implementation.
- Test: `tests/test_anisotropic_scattering.cpp`.

**Interfaces:**
- Consumes: `incoming`, `vertical`, and relative `S2::ThetaPhi(psi, zeta)`.
- Produces: a shared `ScatterParams::MakeAxisymmetricBasis()` helper and a directional `GSATO` whose outgoing wavevector is expressed around `vertical`.

- [ ] **Step 1: Add the failing rotated-axis test.**

Add a deterministic rotation helper in the test and compare `GSATO` values for one canonical geometry and the same geometry after rotating `incoming` and `vertical` together. Use a nontrivial rotation that moves the symmetry axis away from global Z:

~~~
static R3::XYZ RotateY(const R3::XYZ & v, Real angle) {
  const Real c = std::cos(angle);
  const Real s = std::sin(angle);
  return R3::XYZ(c*v.x() + s*v.z(), v.y(), -s*v.x() + c*v.z());
}
~~~

Use `incoming=(0.6, 0.0, 0.8)`, `vertical=(0.0,0.0,1.0)`, `toa=S2::ThetaPhi(0.9,1.1)`, rotate both vectors by `0.7` radians, and assert all four G values and `spol` agree within `1e-11`.

- [ ] **Step 2: Run the focused test to verify it fails.**

Run:

~~~bash
make tests/test_anisotropic_scattering
./tests/test_anisotropic_scattering
~~~

Expected: failure in the new rotated-axis assertion because the current implementation constructs the outgoing direction using the global-Z `OrthoAxes` basis.

- [ ] **Step 3: Add the shared basis helper.**

Declare and implement:

~~~
static void MakeAxisymmetricBasis(const R3::XYZ & incoming,
                                  const R3::XYZ & vertical,
                                  R3::XYZ & meridian,
                                  R3::XYZ & azimuth);
~~~

The helper normalizes both inputs, projects `vertical` into the plane perpendicular to `incoming`, chooses the increasing-theta direction (`incoming*mu - vertical`), and sets `azimuth=incoming.Cross(meridian)`. When the projection magnitude is below `1e-12`, select the first stable global Cartesian axis not parallel to `incoming`, project it into the transverse plane, and normalize. Keep the helper in `ScatterParams` so the kernel's geometry and runtime frame conversion use identical conventions.

- [ ] **Step 4: Use the basis in directional `GSATO`.**

Replace the global-Z transform with:

~~~
MakeAxisymmetricBasis(incoming, vertical, meridian, azimuth);
const R3::XYZ out = in.ScaledBy(std::cos(psi))
                  + meridian.ScaledBy(std::sin(psi) * std::cos(zeta))
                  + azimuth.ScaledBy(std::sin(psi) * std::sin(zeta));
~~~

Keep the existing `XSATO` radiation factors and PSD calls unchanged. The `vertical` argument must continue into all four `PowerSpectralDensity()` calls.

- [ ] **Step 5: Run the focused test to verify it passes.**

Run:

~~~bash
make tests/test_anisotropic_scattering
./tests/test_anisotropic_scattering
~~~

Expected: PASS, including the existing equal/unequal PSD checks and the new rotated-axis GSATO check.

- [ ] **Step 6: Commit.**

~~~bash
git add scatparams.hpp scatparams.cpp tests/test_anisotropic_scattering.cpp
git commit -m "fix: express anisotropic scattering around local vertical"
~~~

### Task 2: Add compact angular storage for anisotropic scatterers

**Files:**
- Modify: `sources.hpp` constructor declarations and comments.
- Modify: `sources.cpp` constructor and `PrepareForSimulation()`.
- Modify: `scatterers.hpp` member declarations.
- Modify: `scatterers.cpp` constructor initialization.
- Test: `tests/test_axisymmetric_scattering.cpp` compile scaffold in Task 3.

**Interfaces:**
- Consumes: the existing `nTOA` static size.
- Produces: `PhononSource(int nraytypes_in, int nraytypes_out, int angular_size)` and a `Scatterer` anisotropic branch with no allocated full-sphere `mPDists`.

- [ ] **Step 1: Write the compact-constructor compile test.**

Create the test file with a small helper that constructs the new kernel through the eventual `Scatterer` path only after `PhononSource::Set_TOA_Array()` has been called. The initial assertion is that `IncomingBinCount()` is positive; the test will compile against the new compact storage seam once Task 3 adds the kernel.

- [ ] **Step 2: Run the new test target to verify it cannot build yet.**

Run:

~~~bash
make tests/test_axisymmetric_scattering
~~~

Expected: failure because `axisymmetric_scattering.hpp` and its implementation do not exist yet.

- [ ] **Step 3: Add the optional angular-size constructor.**

Keep the existing constructor source-compatible:

~~~
PhononSource(int nraytypes_in, int nraytypes_out)
  : PhononSource(nraytypes_in, nraytypes_out, nTOA) {}

PhononSource(int nraytypes_in, int nraytypes_out, int angular_size);
~~~

The three-argument constructor must allocate `mPDists` with `nraytypes_out` distributions only when `angular_size > 0`; `mWholeProbs` always retains its existing scalar output-condition storage. Make `PrepareForSimulation()` iterate only over distributions that exist, then prepare `mWholeProbs` as before.

- [ ] **Step 4: Select compact storage only for unequal anisotropy.**

Change the `Scatterer` base initializer to pass `0` angular size when `par.IsAnisotropic()` is true and `nTOA` otherwise. Do not change event-source construction or the equal-length path.

- [ ] **Step 5: Run the existing build and test targets.**

Run:

~~~bash
make -j2
make tests/test_anisotropic_scattering
./tests/test_anisotropic_scattering
~~~

Expected: the existing isotropic tests remain green and the main executable still links; anisotropic construction is not exercised until Task 4.

- [ ] **Step 6: Commit.**

~~~bash
git add sources.hpp sources.cpp scatterers.hpp
git commit -m "refactor: allow compact scatterer angular storage"
~~~

### Task 3: Build the immutable reduced quadrature kernel

**Files:**
- Create: `axisymmetric_scattering.hpp`.
- Create: `axisymmetric_scattering.cpp`.
- Modify: `Makefile` object list, header dependencies, and test target.
- Test: `tests/test_axisymmetric_scattering.cpp`.

**Interfaces:**
- Consumes: `ScatterParams::GSATO()`, `raytype`, `RandomEngine`, and the local-vertical basis convention from Task 1.
- Produces: the exact `AxisymmetricScatteringKernel` interface defined above.

- [ ] **Step 1: Write the failing reference and cache tests.**

Implement `tests/test_axisymmetric_scattering.cpp` with these checks:

~~~
static ReferenceTotals IntegrateReference(const ScatterParams & params,
                                          Real mu_in,
                                          int npsi, int nzeta);

int main() {
  Elastic::HSneak spectrum(0.8, 0.01, 1.0, 0.5);
  ScatterParams params(spectrum, 1.0, 1.7, 1.25, 0.625);
  AxisymmetricScatteringKernel kernel(params);

  assert(kernel.IncomingBinCount() >= 32);
  assert(kernel.DirectionNodeCount() >= 512);
  assert(kernel.CDFsAreValid(1.0e-12));

  for (Real mu : {-0.8, -0.2, 0.0, 0.4, 0.9}) {
    const ReferenceTotals ref = IntegrateReference(params, mu, 96, 192);
    const Real cached = kernel.GetInverseMeanFreePath(RAY_P, mu);
    assert(std::abs(cached - ref.inverse_p) / ref.inverse_p < 5.0e-3);
    assert(std::abs(kernel.GetConversionWeight(
      AxisymmetricScatteringKernel::GPP, mu) - ref.gpp) / ref.gpp < 5.0e-3);
  }
}
~~~

The reference uses the same normalized solid-angle measure as production but computes every `GSATO()` value directly on a dense midpoint grid. Add checks for both P and S totals, finite CDFs, and monotonic CDF entries.

- [ ] **Step 2: Run the focused target to verify it fails.**

Run:

~~~bash
make tests/test_axisymmetric_scattering
~~~

Expected: compilation failure for the missing kernel header/class.

- [ ] **Step 3: Define cache types and fixed private resolutions.**

Implement the public interface and private structures:

~~~
namespace {
const std::size_t kIncomingBins = 65;  // includes -1 and +1 endpoints
const int kPsiNodes = 48;
const int kZetaNodes = 64;
}

struct AxisymmetricScatteringKernel::Bin {
  Real mu;
  Real inverse_mfp[RAY_NBT];
  Real totals[NUM_CONVERSIONS];
  std::vector<Real> cdf[NUM_CONVERSIONS];
};
~~~

Store shared `S2::ThetaPhi` relative nodes and S/S canonical polarization angles once. Use `std::vector<Bin>` for incoming bins. Do not use global mutable tables.

- [ ] **Step 4: Implement Gauss-Legendre and periodic quadrature.**

Create a helper returning nodes and weights on `[-1,1]` using Newton iteration on Legendre roots. For each node pair use:

~~~
const Real solid_angle_average_weight = 0.5 * psi_weight / kZetaNodes;
~~~

because the `zeta` integral is periodic and the result is normalized by `4*pi`. Build relative nodes as `S2::ThetaPhi(std::acos(cos_psi), zeta)`.

- [ ] **Step 5: Populate each incoming bin.**

For each `mu` construct canonical `incoming=(sqrt(1-mu*mu),0,mu)` and `vertical=(0,0,1)`. For every reduced node call directional `GSATO`, multiply each G value by the normalized quadrature weight, and accumulate:

~~~
bin.totals[GPP] += gpp * weight;
bin.totals[GPS] += gps * weight;
bin.totals[GSP] += gsp * weight;
bin.totals[GSS] += gss * weight;
~~~

Store the weighted prefix sums in each conversion CDF. Set P/S inverse MFPs from the allowed conversion totals, rejecting neither tiny positive totals nor valid strongly anisotropic values. If a total is nonpositive, record inverse MFP zero and let the runtime safe fallback handle it.

- [ ] **Step 6: Implement bracket interpolation and CDF sampling.**

Clamp `mu_in` to `[-1,1]`, locate its two neighboring bins, and linearly interpolate inverse MFPs and conversion totals. For `GetRandomSample()`, first select an allowed conversion, then choose low/high CDF in proportion to its interpolated unnormalized contribution, and binary-search that CDF with the provided `RandomEngine`.

Map conversions to output types exactly as the existing code does:

~~~
const raytype output_types[NUM_CONVERSIONS] = {
  RAY_P, RAY_S, RAY_P, RAY_S
};
~~~

Return the cached `psi`, `zeta`, and canonical polarization. If both candidate bins have zero selected weight, return a deterministic forward P/S sample with zero polarization for the caller's input type.

- [ ] **Step 7: Implement diagnostics and validation.**

`CDFsAreValid(tolerance)` must check finite values, nondecreasing prefixes, positive final totals for every physically allowed conversion, and that every bin's normalized conversion probabilities sum to one for P and S inputs. The public node-count methods return the fixed table dimensions for focused tests. Implement `GetAverageInverseMeanFreePath()`, `GetAverageConversionWeight()`, and `GetAverageDipole()` as quadrature-weighted averages over the incoming bins; these are used only for legacy scalar diagnostics and compact `mWholeProbs` initialization.

- [ ] **Step 8: Add Makefile dependencies and run the target.**

Add `axisymmetric_scattering.o` to `objects`, define its header dependency on `scatparams.hpp`, `probability.hpp`, and `geom.hpp`, and link the new test with `axisymmetric_scattering.cpp`, `scatparams.cpp`, `geom_r3.cpp`, `geom_s2.cpp`, and `elastic.cpp`.

Run:

~~~bash
make tests/test_axisymmetric_scattering
./tests/test_axisymmetric_scattering
~~~

Expected: PASS for cache shape/CDF checks and reference agreement. If the dense reference exposes insufficient production resolution, increase only the private constants and rerun the test; do not add a CLI option.

- [ ] **Step 9: Commit.**

~~~bash
git add axisymmetric_scattering.hpp axisymmetric_scattering.cpp Makefile \
        tests/test_axisymmetric_scattering.cpp
git commit -m "feat: add cached axisymmetric scattering kernel"
~~~

### Task 4: Integrate MFPs, direction sampling, and polarization-frame conversion

**Files:**
- Modify: `scatterers.hpp`.
- Modify: `scatterers.cpp`.
- Modify: `phonons.cpp` scattering call.
- Modify: `Makefile` dependency graph to include `axisymmetric_scattering.hpp`.
- Test: `tests/test_axisymmetric_scattering.cpp`.

**Interfaces:**
- Consumes: immutable `AxisymmetricScatteringKernel`, `incoming`, `vertical`, `incoming_polarization`, and explicit RNG.
- Produces: anisotropic `GetRandomPathLength()` and geometry-aware sampling with a physically correct relative `Phonon` for the existing transform.

- [ ] **Step 1: Add the failing frame-conversion test.**

Add a deterministic helper test using identical RNG seeds for two incoming S
polarization angles. Apply the returned relative phonons to otherwise identical
current phonons and compare the physical results:

~~~
RandomEngine rng_a(0x1234);
RandomEngine rng_b(0x1234);
Phonon rel_a = scatterer.GetRandomScatteredRelativePhonon(
    RAY_S, incoming, vertical, 0.0, rng_a);
Phonon rel_b = scatterer.GetRandomScatteredRelativePhonon(
    RAY_S, incoming, vertical, Geometry::Pi90, rng_b);
Phonon out_a(S2::ThetaPhi(incoming), RAY_S);
Phonon out_b(S2::ThetaPhi(incoming), RAY_S);
out_a.SetPolarization(0.0);
out_b.SetPolarization(Geometry::Pi90);
out_a.Transform(rel_a);
out_b.Transform(rel_b);
assert(out_a.GetDirection().XYZ().DistFrom(out_b.GetDirection().XYZ()) < 1.0e-11);
assert(std::abs(out_a.DirectionOfMotion().Dot(out_a.GetDirection().XYZ())) < 1.0e-11);
assert(std::abs(out_b.DirectionOfMotion().Dot(out_b.GetDirection().XYZ())) < 1.0e-11);
~~~

The returned relative `Phi()` values should differ because they are expressed
in different incoming-polarization frames, while the transformed physical
outgoing directions should agree.

- [ ] **Step 2: Run the focused test to verify it fails.**

Run:

~~~bash
make tests/test_axisymmetric_scattering
./tests/test_axisymmetric_scattering
~~~

Expected: failure or missing-symbol compilation until the Scatterer integration and frame conversion exist.

- [ ] **Step 3: Add kernel ownership and compact constructor branching.**

Include `axisymmetric_scattering.hpp`, add:

~~~
std::unique_ptr<AxisymmetricScatteringKernel> mpAxisymmetric;
~Scatterer();
~~~

In the constructor, build `mpAxisymmetric` only when `mParams.IsAnisotropic()` is true. For that branch, set `mMeanFreeP` from the reciprocal of `GetAverageInverseMeanFreePath()`, set `mDipoles` from `GetAverageDipole()`, and populate `mWholeProbs[IN_P]` with the average GPP/GPS weights and `mWholeProbs[IN_S]` with the average GSP/GSS weights before calling `PrepareForSimulation()`. Keep the current `PopulateProbDists`, `PopulateWholeProbs`, `ComputeMFPs`, and `ComputeDipoles` sequence byte-for-byte in the isotropic branch where practical.

- [ ] **Step 4: Delegate directional MFP lookup.**

Replace the anisotropic body of `GetDirectionalMeanFreePath()` with:

~~~
const Real mu = incoming.UnitElse(vertical).Dot(vertical.UnitElse(
    R3::XYZ(0,0,1)));
const Real inverse = mpAxisymmetric->GetInverseMeanFreePath(intype, mu);
return inverse > 0.0 ? 1.0 / inverse : mMeanFreeP[intype];
~~~

Keep the explicit-MFP and isotropic branches in `GetRandomPathLength()` before this lookup. Use `-std::log(1.0 - rng.Uniform01())` exactly as the current code.

- [ ] **Step 5: Add the incoming-polarization overloads.**

Keep the current overload and delegate:

~~~
Phonon Scatterer::GetRandomScatteredRelativePhonon(
    raytype intype, const R3::XYZ & incoming,
    const R3::XYZ & vertical, RandomEngine & rng) {
  return GetRandomScatteredRelativePhonon(
      intype, incoming, vertical, 0.0, rng);
}
~~~

Change the propagation call in `phonons.cpp` to pass `mPol`. The no-geometry overload uses incoming `(0,0,1)`, vertical `(0,0,1)`, and polarization `0.0` only for compatibility diagnostics.

- [ ] **Step 6: Convert the cached meridional sample into Transform coordinates.**

Implement the following conversion in `scatterers.cpp`:

~~~
const R3::XYZ in = incoming.UnitElse(R3::XYZ(0,0,1));
const R3::XYZ up = vertical.UnitElse(R3::XYZ(0,0,1));
R3::XYZ meridian, azimuth;
ScatterParams::MakeAxisymmetricBasis(in, up, meridian, azimuth);

const Real sp = std::sin(sample.psi);
const Real cp = std::cos(sample.psi);
const Real cz = std::cos(sample.zeta);
const Real sz = std::sin(sample.zeta);
const R3::XYZ out = in.ScaledBy(cp)
                  + meridian.ScaledBy(sp * cz)
                  + azimuth.ScaledBy(sp * sz);
const R3::XYZ rtheta = meridian.ScaledBy(cp * cz)
                     + azimuth.ScaledBy(cp * sz)
                     - in.ScaledBy(sp);
const R3::XYZ rphi = meridian.ScaledBy(-sz) + azimuth.ScaledBy(cz);
const R3::XYZ output_pol = rtheta.ScaledBy(std::cos(sample.polarization))
                         + rphi.ScaledBy(std::sin(sample.polarization));
~~~

Build `AA=R3::OrthoAxes(in.Theta(), in.Phi(), incoming_polarization)`, express `out` and `output_pol` in AA components, derive relative `theta/phi`, then derive relative polarization with `BB=R3::OrthoAxes(rel_theta,rel_phi,0)`:

~~~
const Real rel_pol = std::atan2(
    pol_in_AA.Dot(BB.E2()), pol_in_AA.Dot(BB.E1()));
Phonon result(S2::ThetaPhi(rel_theta, rel_phi), sample.output_type);
result.SetPolarization(sample.output_type == RAY_P ? 0.0 : rel_pol);
~~~

For `GPS`, the cache's canonical polarization is zero; for `GSS`, use the cached Sato polarization; for P outputs polarization is ignored. Clamp the AA-frame longitudinal component to `[-1,1]` before `acos` and preserve the existing phonon singular-direction nudging.

- [ ] **Step 7: Update scatter-pattern diagnostics.**

Make `test_random_rayset()` call the geometry-aware sampler with incoming and vertical `(0,0,1)` when the kernel exists. Keep the existing isotropic output path for equal lengths. This prevents compact anisotropic scatterers from indexing empty `mPDists` distributions.

- [ ] **Step 8: Run focused tests and a small anisotropic process.**

Run:

~~~bash
make -j2
./tests/test_axisymmetric_scattering
./tests/test_anisotropic_scattering
./main --reports=ALL_OFF --num-phonons=32 --toa-degree=3 \
  --source=EXPL --source-loc=425.54,-169.53,-1.02 \
  --frequency=2.0 --timetolive=10 --binsize=1.0 --grid-compiled=1 \
  --range=1200 --flatten --model-args=0.8,0.01,0.5,0.2,50 \
  --seis-p2p=425.54,-169.53,0.98,-390.04,-167.18,1.457,1.0,2.0,2.0,4 \
  --scatter-horizontal=1.25 --scatter-vertical=0.625 \
  --workers=1 --seed=0x5eedc0de12345678
~~~

Expected: the process completes without invalidity and without a per-event 5.2-million-direction scan.

- [ ] **Step 9: Commit.**

~~~bash
git add scatterers.hpp scatterers.cpp phonons.cpp Makefile
git commit -m "feat: use cached axisymmetric scattering during propagation"
~~~

### Task 5: Add seeded process regression and complete numerical verification

**Files:**
- Modify: `tests/test_anisotropic_scattering.cpp`.
- Modify: `tests/test_axisymmetric_scattering.cpp`.
- Modify: `tests/test_parallel_reproducibility.sh`.
- Modify: `Makefile` test recipe.

**Interfaces:**
- Consumes: public kernel diagnostics and integrated `Scatterer` behavior from Tasks 1–4.
- Produces: reproducible small-run coverage for unequal anisotropy, overrides, no-deflection, and resolution convergence.

- [ ] **Step 1: Add direct axisymmetric invariance checks.**

For fixed `mu_in`, compare `GetInverseMeanFreePath()` and all conversion weights for canonical `(sqrt(1-mu*mu),0,mu)` and an azimuth-rotated incoming vector with the same `mu`. Require relative differences below `1e-12`; this verifies the table depends on polar angle, not absolute azimuth.

- [ ] **Step 2: Add resolution/reference assertions.**

Use `IntegrateReference(params, mu, 96, 192)` at `mu=-0.8,-0.2,0.0,0.4,0.9` for both P and S. Assert all positive totals agree within `5e-3`, then verify that the cache's interpolation at `mu=0.13` lies between the neighboring bin values and is finite.

- [ ] **Step 3: Add override and no-deflection checks.**

Construct a `Scatterer` with anisotropic parameters, call
`Scatterer::OverrideMFP(12.0, 34.0)`, and verify sampled path-length means are
consistent with those values over a fixed seeded sample. Run this override
case before the no-deflection case. Then enable `SetNoDeflect()` and assert the
relative direction is `theta=0` and output type matches input. Place both
static-override checks at the end of the native test process; the process exits
immediately afterward, so no later test depends on a reset API.

- [ ] **Step 4: Extend the process reproducibility shell test.**

Refactor `run_case()` to accept an output label and optional extra arguments. Keep the existing isotropic one-worker/two-worker comparison, then add:

~~~bash
run_case anisotropic 1 --scatter-horizontal=1.25 --scatter-vertical=0.625
run_case anisotropic 2 --scatter-horizontal=1.25 --scatter-vertical=0.625
~~~

Use 512 phonons, TOA degree 3, reports off, and the existing fixed seed. Diff only the `Loss surfaces:`, `Timeout:`, and `Invalidity:` summary lines between the two worker counts.

- [ ] **Step 5: Add all native targets to `make test`.**

Ensure the target list builds and runs both anisotropic executables before the process-level shell tests. Do not add a machine-specific time assertion.

- [ ] **Step 6: Run the complete focused suite.**

Run:

~~~bash
make -j2
make test
~~~

Expected: all native, process-level, anisotropic, and script tests pass.

- [ ] **Step 7: Commit.**

~~~bash
git add tests/test_anisotropic_scattering.cpp \
        tests/test_axisymmetric_scattering.cpp \
        tests/test_parallel_reproducibility.sh Makefile
git commit -m "test: verify axisymmetric anisotropic reproducibility"
~~~

### Task 6: Update documentation and roadmap

**Files:**
- Modify: `docs/DEVELOPMENT.md` anisotropic implementation section.
- Modify: `docs/MANUAL.md` anisotropic option/workflow wording to describe the
  unchanged CLI semantics and the new axisymmetric runtime behavior.
- Modify: `PLANS.md` current baseline and completed/follow-up sections.

**Interfaces:**
- Consumes: verified behavior from Tasks 1–5.
- Produces: durable coordinate, accuracy, performance, and limitation documentation without inventing new CLI options.

- [ ] **Step 1: Document the reduced coordinates and cache.**

Describe `mu_in`, `psi`, `zeta`, normalized solid-angle weights, local `up`, meridian-to-polarization-frame conversion, and the pre-worker freeze point. State that trajectories remain 3D and that equal lengths use the existing isotropic path.

- [ ] **Step 2: Document numerical and diagnostic limits.**

State the private fixed quadrature resolution, reference-validation tolerance, that scalar printed MFP/dipole values are incoming-direction averages, and that the axisymmetric approximation does not support three unequal axes or per-cell CLI axes.

- [ ] **Step 3: Update the roadmap.**

Move the anisotropic performance bottleneck from the active follow-up list to completed work, record the reduced-kernel validation evidence, and retain the future analytic-azimuth and fully anisotropic tensor possibilities as explicit follow-ups.

- [ ] **Step 4: Run documentation/script checks.**

Run:

~~~bash
git diff --check
make test-plotting
~~~

Expected: no whitespace errors and the existing plotting regression passes or reports only an unavailable optional toolchain.

- [ ] **Step 5: Commit.**

~~~bash
git add docs/DEVELOPMENT.md docs/MANUAL.md PLANS.md
git commit -m "docs: describe axisymmetric anisotropic scattering"
~~~

### Task 7: Final verification, integration, and branch cleanup

**Files:**
- Verify: all changed source, test, recipe, and documentation files.

- [ ] **Step 1: Run the required final verification commands.**

Run on the final feature-branch tree:

~~~bash
make -j2
make test
make test-plotting
git diff --check anisotropic-lopnor-script...HEAD
~~~

Also run the 100-phonon unequal-anisotropy probe with the repository's direct scenario command and record construction/simulation completion evidence. Do not claim success without checking each exit status.

- [ ] **Step 2: Review the complete diff and status.**

Run:

~~~bash
git status --short --branch
git diff --stat anisotropic-lopnor-script...HEAD
git log --oneline --decorate -8
~~~

Confirm no generated data, logs, binaries, or temporary PDF renderings are staged, and confirm the pre-existing recipe changes are present.

- [ ] **Step 3: Commit any final implementation corrections.**

Use a focused commit with an explanatory message if verification finds a real issue. Rerun the affected focused test and then the full required suite.

- [ ] **Step 4: Merge the feature branch back into `anisotropic-lopnor-script`.**

After all verification passes, switch to the original branch and merge without rewriting history:

~~~bash
git switch anisotropic-lopnor-script
git merge --no-ff axisymmetric-anisotropic-scattering \
  -m "merge: add cached axisymmetric anisotropic scattering"
~~~

Resolve only genuine merge conflicts, preserve the recipe changes, and rerun `git status --short --branch` after the merge.

- [ ] **Step 5: Delete the temporary feature branch.**

Once the merged branch contains the verified commits:

~~~bash
git branch -d axisymmetric-anisotropic-scattering
~~~

Do not delete the original `anisotropic-lopnor-script` branch or push to a remote.

- [ ] **Step 6: Report the final handoff.**

Provide the merged branch name, final commit/merge summary, test commands and results, anisotropic probe output directory, and any optional plotting or sanitizer limitation observed during verification.
