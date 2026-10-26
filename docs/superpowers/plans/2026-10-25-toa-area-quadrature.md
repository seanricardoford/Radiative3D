# Area-Weighted Take-off-Angle Quadrature Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in area-weighted take-off-angle quadrature, `--toa-quadrature=area`, to the scatterer tables, the event-source tables and the scatterer reports. Validate it, measure its scenario-level effect, and document it. The default `equal` mode stays byte-identical to `319c8dd`.

**Architecture:**
- `S2::TesselSphere` optionally computes each direction's triangle solid angle Ω_k.
- `PhononSource` holds the run's quadrature mode and supplies w_k = Ω_k·N/4π.
- The weights are applied once, where the tables are built. Every sampler, majorant, MFP and dipole inherits them unchanged.
- The report integrator takes the report set's solid angles explicitly.

**Tech Stack:** C++11 (GCC 16, `-O3`, default `-ffp-contract=fast` outside `aniso_axi.cpp`); POSIX sh and awk scripts; the `tests/anisotropy_tests.cpp` harness.

**Spec:** `docs/superpowers/specs/2026-10-25-toa-area-quadrature-design.md`. It implements the plan in `docs/toa-quadrature-bias.md`; read both.

## Global Constraints

- **Default mode `equal` behaves exactly as `319c8dd`.**
  - Traces, `.octv` files and normalized stdout are byte-identical, for isotropic and anisotropic runs.
  - Isotropic runs stay byte-identical to `aaec540` (`tests/test_aniso_isotropic_regression.sh`).
- **Equal-mode code paths keep today's arithmetic.**
  - Never insert a weight factor into an equal-mode expression that is later summed. GCC contracts multiply-adds in `scatterers.cpp`.
  - Branch on the mode instead.
- **Weights:** area mode uses w_k = Ω_k·N/(4π) everywhere. N is the size of the direction set being sampled or integrated, and Ω_k comes from `TesselSphere::SolidAngles()`.
- **Random numbers** come only from `RandomUniform()` and `RandomInteger()`. Model construction draws none. Both modes make the same draw calls; area mode changes table contents only.
- **Determinism, in both modes:**
  - construction parallelism never changes a table, a report value or an output;
  - outputs are identical between `--threads 1` and `--threads N` at a fixed seed.
- **Bounds:** every `AxiProposal` bound stays a rigorous upper bound. A violated bound throws `std::logic_error` and is never clamped.
- **Memory:**
  - area mode adds one `Real` per take-off direction, to the global set only;
  - `ReportDirs()`'s degree-6 set also carries 81,920 solid angles in both modes (about 0.66 MB, once per process);
  - no per-scatterer memory changes.
- **Output:**
  - equal mode prints no new stdout lines;
  - area mode prints exactly three new lines, listed in Task 4;
  - `out_mparams.octv` is unchanged in both modes.
- **Files not to touch:**
  - do not edit the benchmark numbers in `README.md`, `doc/MANUAL.md`, `docs/multithreaded-simulation-handoff.md` or `docs/anisotropy-performance.md`;
  - do not change `tests/run_lopnor_thread_case.sh`, `tests/test_lopnor_thread_equivalence.sh` or `tests/compare_axi_general.sh`.
- **Reference binaries.** Both are already built; do not rebuild them.
  - `BASELINE_EXE=/tmp/r3d-baseline/main`, built from `aaec540`.
  - `REF_EXE=/tmp/r3d-quad-ref/main`, built from `319c8dd`, the start of this plan.
  - If `/tmp` was swept, rebuild the REF binary with `mkdir -p /tmp/r3d-quad-ref-src /tmp/r3d-quad-ref && git archive 319c8dd | tar -x -C /tmp/r3d-quad-ref-src && make -C /tmp/r3d-quad-ref-src 'revision="\"319c8dd-ref\""' && cp /tmp/r3d-quad-ref-src/main /tmp/r3d-quad-ref/main`.
- **Full check**, run at the end of every code task (about 90 s at the start):
  `export BASELINE_EXE=/tmp/r3d-baseline/main; make && for t in tests/test_*.sh; do "$t" || exit 1; done`.
- **Default-mode reference check,** run at the end of Tasks 2, 3 and 4:
  `REF_EXE=/tmp/r3d-quad-ref/main tests/compare_aniso_reference.sh`. Every case must print PASS.
- **Model tiers:** implementers and per-task reviewers use sonnet; the final whole-branch review uses opus. Per this repository's history, do not use haiku for these tasks. Every number written into a document must trace to captured output.

## Review Focus

1. **Default-mode anisotropic byte-identity.** A forgotten branch, or a weight multiplied into a summed expression, would shift anisotropic report digits or traces in the default mode. Only the reference compare catches that.
   - Tests: `tests/compare_aniso_reference.sh` against `REF_EXE` (Tasks 2, 3, 4).
   - Test: explicit `--toa-quadrature=equal` strictly identical to no flag (Task 4).
2. **`--overridemfp` with area mode.** The redraw-until-accepted loop must still land on the weighted target distribution.
   - Test: `isolated_area_overridemfp` (Task 3).
3. **Curved mappings in area mode** (crust pinch, `RAE_CURVED`). The weights live in the ray frame, so the mapping must not matter. Runs must stay thread-invariant and complete.
   - Test: the crust-pinch area-mode checks in `tests/test_toa_quadrature.sh` (Task 4).
4. **High degrees in area mode.**
   - Degree 8: the solid-angle sum must hold. Test: `tessel_solid_angles` (Task 1).
   - Degree 9: the build must stay usable, with about 42 MB extra and a few tenths of a second. Test: the degree-9 build check in `tests/compare_toa_quadrature.sh` (Task 5).
5. **`--event-test --toa-quadrature=area`.** It must neither crash nor change the displayed pattern.
   - Test: in `tests/test_toa_quadrature.sh`, the output matches plain `--event-test` (Task 4).

---

### Task 1: Per-direction solid angles in `S2::TesselSphere`

**Files:**
- Modify: `geom_s2.hpp` (classes `Triangle`, `TesselSphere`), `geom_s2.cpp` (`Triangle` constructor, `populate_lists()`, `subdivide()`; `TesselSphere` constructor, `init_icosahedron()`, `grab_centerpoints_from_tri()`)
- Test: `tests/anisotropy_tests.cpp`. Add a new section banner, "TOA QUADRATURE TESTS", just above the TEST RUNNER section.

**Interfaces:**
- Produces:

```cpp
// geom_s2.hpp
class Triangle {   // existing class; additions only
public:
  Triangle(Node, Node, Node, int degree, bool solidAngles = false);
  const std::vector<Real> & SolidAngles() const;
      // One entry per leaf, in Centers order; empty unless solidAngles.
};
class TesselSphere : public ThetaPhi_Set {   // existing class; additions only
public:
  TesselSphere(TESS_Algorithm algo, int degree, bool withSolidAngles = false);
  bool HasSolidAngles() const;
  const std::vector<Real> & SolidAngles() const;
      // Omega_k (sr), index-aligned with *this; empty unless HasSolidAngles().
};
```

- [ ] **Step 1: Write the failing tests** in the new section.

```cpp
TEST(tessel_solid_angles)
  // writes build/toa_quadrature_solid_angles.txt:
  //   "degree N sum_rel_err min_over_mean max_over_mean rms", one line per degree
  for d in 0..8:
    S2::TesselSphere plain(S2::TESS_ICO, d), sa(S2::TESS_ICO, d, true);
    CHECK(!plain.HasSolidAngles() && plain.SolidAngles().empty());
    CHECK(sa.HasSolidAngles() && sa.size() == plain.size()
          && sa.SolidAngles().size() == sa.size());
    for each k:
      CHECK(sa[k].Theta() == plain[k].Theta() && sa[k].Phi() == plain[k].Phi());
      CHECK(sa.SolidAngles()[k] > 0);
    sum = plain loop over SolidAngles();  mean = sum / size;
    CHECK(std::fabs(sum - 4*Pi) <= 1e-12 * 4*Pi);
    if d == 0: CHECK(max/min - 1 <= 1e-14);          // 20 congruent faces
    if d >= 6: CHECK(min/mean in [0.9276, 0.9278]);
               CHECK(max/mean in [1.2064, 1.2067]);
               CHECK(rms(Omega/mean - 1) in [0.0864, 0.0865]);
TEST(tessel_solid_angles_mirror)
  // degree 3, with solid angles
  for each k, and for each of the images (-x, y, z) and (x, -y, z) of sa[k].XYZ():
    j = sa.GetBestIndexFromPoint(image as S2::ThetaPhi);
    CHECK(sa[j].XYZ().DistFrom(image) < 1e-9);
    CHECK(std::fabs(sa.SolidAngles()[j] - sa.SolidAngles()[k])
          <= 1e-12 * sa.SolidAngles()[k]);
```

  The spike measured these: at degree 6, min/mean = 0.9277266, max/mean = 1.2065026, rms = 0.0864447. The sum's relative error stays ≤ 2e-13 up to degree 8.
- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh tessel_solid_angles`.** Expect a compile failure on the three-argument constructor.
- [ ] **Step 3: Implement.**
  - When `solidAngles` is set, each degree-0 `Triangle` pushes Ω = 2·atan2(|n1·(n2×n3)|, 1 + n1·n2 + n2·n3 + n3·n1), computed from its own three vertex nodes.
  - `subdivide()` concatenates the children's lists exactly as it concatenates `_CenterPoints`, so the two stay index-aligned.
  - `TesselSphere` passes the flag down, appends each face's list in `grab_centerpoints_from_tri()` (or the equivalent), and keeps a `std::vector<Real>`.
  - Do not touch any `Node` arithmetic that produces the centers.
- [ ] **Step 4: Run `tests/test_anisotropy_unit.sh tessel_solid_angles`.** Both tests pass. Then run the full check.
- [ ] **Step 5: Commit.**

```bash
git add geom_s2.hpp geom_s2.cpp tests/anisotropy_tests.cpp
git commit -m "Compute per-direction solid angles in TesselSphere on request"
```

---

### Task 2: Quadrature mode in `PhononSource` and area-weighted report integration

**Files:**
- Modify:
  - `sources.hpp`, `sources.cpp` (mode state and API);
  - `model.hpp` (`mpTOA_Array` becomes `S2::TesselSphere *`);
  - `main.cpp` (the `--event-test` block's `pTOAA` becomes `S2::TesselSphere *`; it stays in equal mode);
  - `scatterers.hpp`, `scatterers.cpp` (`IntegrateRows`, `IntegrateRateAndDipole`, `IntegratedRate`, `ComputeReportStats`, `ReportDirs`).
- Test: `tests/anisotropy_tests.cpp`.

**Interfaces:**
- Consumes: Task 1's `TesselSphere(algo, degree, true)`, `HasSolidAngles()`, `SolidAngles()`.
- Produces:

```cpp
// sources.hpp, class PhononSource
public:
  enum toa_quadrature_e {TOA_EQUAL, TOA_AREA};
  static void Set_TOA_Array(S2::TesselSphere * toa, toa_quadrature_e quad = TOA_EQUAL);
      // Replaces Set_TOA_Array(S2::S2Set*). TOA_AREA requires toa->HasSolidAngles();
      // otherwise throws std::logic_error and leaves the previous array and mode unchanged.
  static toa_quadrature_e Get_TOA_Quadrature();
protected:
  static Real TOAWeight(Index k);
      // TOA_EQUAL: 1.0.  TOA_AREA: SolidAngles()[k] * (nTOA / (4*Geometry::Pi)),
      // the scale computed once in Set_TOA_Array().

// scatterers.hpp, class Scatterer
public:
  Real IntegratedRate(raytype intype, const R3::XYZ & dir, const R3::XYZ & pol,
                      const S2::S2Set & dirs,
                      const std::vector<Real> * solidAngles = nullptr) const;
  void IntegrateRateAndDipole(raytype intype, const R3::XYZ & dir, const R3::XYZ & pol,
                              const S2::S2Set & dirs, Real & rate, Real & dipole,
                              const std::vector<Real> * solidAngles = nullptr) const;
protected:
  void IntegrateRows(raytype intype, const R3::XYZ & dir, const R3::XYZ & pol,
                     const S2::S2Set & dirs, const std::vector<GRow> & rows,
                     const std::vector<Real> * solidAngles,
                     Real & rate, Real & dipole) const;
      // solidAngles == nullptr: today's loop, unchanged.
      // Otherwise its size must equal dirs.size() (std::logic_error if not),
      // with w_k = (*solidAngles)[k] * n/(4*Pi):
      //   rate   = Sum_k w_k*Sum_conv(G*R) / n
      //   dipole = Sum w_k*G*R*cos(psi) / Sum w_k*G*R   (0 when the denominator is 0)
```

- Test-side helpers, produced in this task and reused by Task 3:

```cpp
static const S2::TesselSphere & UseTOADegree(int degree,
    PhononSource::toa_quadrature_e quad = PhononSource::TOA_EQUAL);
    // Cached spheres are now built with solid angles; calls Set_TOA_Array(slot, quad).
struct AreaQuadratureScope {
  explicit AreaQuadratureScope(int degree);  // UseTOADegree(degree, TOA_AREA)
  ~AreaQuadratureScope();                    // UseTOADegree(degree): back to TOA_EQUAL
  const S2::TesselSphere & toas;
};
static Real ExactRate(const ScatterParams & sp, raytype intype, const R3::XYZ & dir,
                      const R3::XYZ & pol, int nTheta, int nPhi);
    // Independent reference for IntegratedRate(), with the same integrand: GSATO() times
    // AnisoRatio(k_out*e_out - k_in*dir), summed over intype's two conversions, where
    // e_out = AA.Express(t), AA = OrthoAxes(dir.Theta(), dir.Phi(),
    // RotationForS1(dirTP, pol)), and t runs over the ray-frame sphere.
    // Integrated as Gauss-Legendre in cos(psi) (nTheta nodes, Newton iteration to
    // 1e-15) times the midpoint rule in phi (nPhi nodes), then divided by 4*Pi.
```

- [ ] **Step 1: Write the failing tests.**

```cpp
TEST(quadrature_mode_switch)
  S2::TesselSphere plain(S2::TESS_ICO, 3), withSA(S2::TESS_ICO, 3, true);
  struct WeightProbe : public PhononSource {static Real W(Index k) {return TOAWeight(k);}};
  PhononSource::Set_TOA_Array(&withSA, PhononSource::TOA_AREA);
  CHECK(PhononSource::Get_TOA_Quadrature() == PhononSource::TOA_AREA);
  sum over k of W(k): CHECK(NearRel(sum, (Real) withSA.size(), 1e-12));
  each k: CHECK(NearRel(W(k), withSA.SolidAngles()[k]*(withSA.size()/(4*Pi)), 1e-14));
  bool threw = false; try {PhononSource::Set_TOA_Array(&plain, PhononSource::TOA_AREA);}
                      catch (const std::logic_error &) {threw = true;}
  CHECK(threw && PhononSource::Get_TOA_Quadrature() == PhononSource::TOA_AREA);  // unchanged
  PhononSource::Set_TOA_Array(&withSA);
  CHECK(Get_TOA_Quadrature() == TOA_EQUAL);  each k: CHECK(W(k) == 1.0);
  UseTOADegree(5);   // leave the usual state for later tests
TEST(integrated_rate_solid_angles)
  const S2::TesselSphere & d6 = UseTOADegree(6);
  ScatterParams sp = LopNorCrust(2, 2, 2);   // isotropic, a = 2 km
  Scatterer s(sp, Scatterer::SAMPLER_GENERAL);
  U = (0,0,1), E = (1,0,0);
  std::vector<Real> bad(d6.size()-1, 1.0);
  CHECK(IntegratedRate(RAY_P, U, E, d6, &bad) throws std::logic_error);
  CHECK(NearRel(ExactRate(sp, RAY_P, U, E, 400, 64),
                ExactRate(sp, RAY_P, U, E, 800, 128), 1e-12));   // reference self-check
  for in in {RAY_P, RAY_S}:   // MFP relative error = exact/rate - 1
    x = ExactRate(sp, in, U, E, 400, 64);
    CHECK(std::fabs(x/s.IntegratedRate(in, U, E, d6, &d6.SolidAngles()) - 1) <= 1e-4);
    CHECK(std::fabs(x/s.IntegratedRate(in, U, E, d6) - 1) >= 3e-3);  // discriminates
TEST(area_report_directional_mfps)
  // writes build/toa_quadrature_report.txt: "state equal_mfp area_mfp" for all 8 states
  AreaQuadratureScope q(5);
  Scatterer s(LopNorCrust(2, 2, 0.5), Scatterer::SAMPLER_GENERAL);
  // Order P0 P45 P90 S0 S45SV S45SH S90SV S90SH; docs/toa-quadrature-bias.md's area column:
  expected[8] = {1602.7, 0, 1203.2, 754.7, 0, 0, 442.8, 374.7};
  for i in {0, 2, 3, 6, 7}: CHECK(NearRel(s.GetReportMFPs()[i], expected[i], 2e-4));
  // The equal_mfp column comes from a second Scatterer built after leaving the scope
  // (or under UseTOADegree(5)).
TEST(area_report_thread_invariant)
  AreaQuadratureScope q(5);
  for spec in {LopNorCrust(2,2,0.5), LopNorCrust(2,1,0.5)}:
    SetBuildThreads(1) -> Scatterer a; SetBuildThreads(7) -> Scatterer b; SetBuildThreads(1);
    CHECK(a.GetReportMFPs() == b.GetReportMFPs());   // element-wise ==
TEST(area_p_frame_residual)
  // writes build/toa_quadrature_pframe.txt: "layer rot_deg equal_mfp area_mfp",
  // plus one "layer exact exact_mfp" line per layer
  const S2::TesselSphere & d6 = UseTOADegree(6);
  layers: crust  {vp 6.40, vs 3.63, ax = ay = 5.0, az 0.5, kappa 0.3, theta 60},
          mantle {vp 8.07, vs 4.63, ax = ay = 7.0, az 0.7, kappa 0.5, theta 89};   phi 20 deg
  per layer: sp = LopNorCrust(a, a, az, kappa, 0.01, vp, vs);
             Scatterer s(sp, SAMPLER_GENERAL);  dir = Dir(theta, 20);
    for i in 0..23: pol = R3::OrthoAxes(theta*DtoR, 20*DtoR, 2*Pi*i/24).S1();   // radians
                    mfpA[i] = 1/s.IntegratedRate(RAY_P, dir, pol, d6, &d6.SolidAngles());
                    mfpE[i] = 1/s.IntegratedRate(RAY_P, dir, pol, d6);
    CHECK(max(mfpA)/min(mfpA) - 1 < 1e-3);   // spike: 1.3e-4 crust, 2.8e-4 mantle
    CHECK(max(mfpE)/min(mfpE) - 1 > 4e-2);   // spike: 5.3% crust, 11.2% mantle
    exact = 1/ExactRate(sp, RAY_P, dir, pol(i=0), 1000, 1000);
    each i: CHECK(std::fabs(mfpA[i]/exact - 1) < 1e-3);
```

  The `area_report_directional_mfps` checks fail before the implementation, because area mode does not exist yet.
- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.** Expect compile failures on `toa_quadrature_e`, the new `Set_TOA_Array` and `IntegratedRate`'s fifth argument.
- [ ] **Step 3: Implement the `PhononSource` API.**
  - `Set_TOA_Array` validates first, then sets `pTOA`, `nTOA`, the mode, the solid-angle pointer and the scale `nTOA/(4*Pi)`.
  - Change `mpTOA_Array`'s type in `model.hpp`. `model.cpp`'s existing `Set_TOA_Array(mpTOA_Array)` call then compiles unchanged, in equal mode, until Task 4.
  - Change the `--event-test` pointer type in `main.cpp`.
- [ ] **Step 4: Implement the weighted `IntegrateRows`** as a separate loop.
  - The `nullptr` path stays today's code, byte for byte, in the same summation order.
  - `IntegrateRateAndDipole` and `IntegratedRate` forward `solidAngles`.
  - `ReportDirs()` builds `TesselSphere(S2::TESS_ICO, 6, true)`.
  - `ComputeReportStats()` passes `&dirs.SolidAngles()` when `Get_TOA_Quadrature() == TOA_AREA`, and `nullptr` otherwise.
  - Update the header comments of the functions touched.
- [ ] **Step 5: Update `tests/anisotropy_tests.cpp` helpers.**
  - `UseTOADegree` gains the mode parameter and builds spheres with solid angles.
  - Add `AreaQuadratureScope` and `ExactRate`.
  - Existing call sites need no change.
- [ ] **Step 6: Run `tests/test_anisotropy_unit.sh`.** Everything passes, including every pre-existing test. Then run the full check, then the default-mode reference check: every case prints PASS.
- [ ] **Step 7: Commit.**

```bash
git add sources.hpp sources.cpp model.hpp main.cpp scatterers.hpp scatterers.cpp tests/anisotropy_tests.cpp
git commit -m "Add the take-off-angle quadrature mode and area-weighted report integration"
```

---

### Task 3: Area-weighted scatterer and event-source tables

**Files:**
- Modify:
  - `scatterers.cpp`: `PopulateProbDists()`; header comments of `PopulateProbDists()`, `ComputeMFPs()` (replace the "TesselSphere needs to code area elements" TODO) and `ComputeDipoles()`;
  - `events.cpp`: the `ShearDislocation` constructor;
  - `aniso_axi.hpp`, `aniso_axi.cpp`: one comment each, stating that the per-cell sums and in-cell draws inherit the weights and the bounds are unchanged. No code change.
- Test: `tests/anisotropy_tests.cpp`.

**Interfaces:**
- Consumes: Task 2's `Get_TOA_Quadrature()`, `TOAWeight(k)`, `IntegratedRate(..., solidAngles)`, `UseTOADegree(d, quad)`, `AreaQuadratureScope`, `ExactRate`.
- Produces:
  - in area mode, every `mPDists[c]` entry k equals `TOAWeight(k)` times its equal-mode value (`Scatterer` and `ShearDislocation`);
  - `m_spol` is unchanged in both modes;
  - test helpers `CheckDistributionMatchesTarget(..., const std::vector<Real> * solidAngles = nullptr)` and `BuildIsoPDists(sp, toas, pdists, const std::vector<Real> * solidAngles = nullptr)`, each multiplying its weight or entry by w_k = Ω_k·n/(4π) when the pointer is non-null;
  - `CheckAxiBoundsValid(..., const std::vector<Real> * solidAngles = nullptr)`, which forwards to `BuildIsoPDists`;
  - `static void CheckAxiMatchesGeneral(const ScatterParams & sp)`, the factored body of `axi_matches_general`, which keeps that test's behavior unchanged.

- [ ] **Step 1: Write the failing tests.**

```cpp
TEST(area_tables_weighted)
  // degree 3; ScatterProbe : Scatterer exposes
  //   Real Entry(int c, Index k) const {return mPDists[c].GetRangeMagnitude(k, k+1);}
  //   Real Spol(Index k) const {return m_spol[k];}
  for spec in {LopNorCrust(0.5,0.5,0.5), LopNorCrust(2,2,0.5)}:
    UseTOADegree(3) -> ScatterProbe e(spec);
    UseTOADegree(3, TOA_AREA) -> ScatterProbe a(spec), with w_k from the sphere's solid angles;
    UseTOADegree(3);
    each c in 0..3, k with e.Entry(c,k) > 0:
      CHECK(NearRel(a.Entry(c,k)/e.Entry(c,k), w_k, 1e-9));
    each k: CHECK(a.Spol(k) == e.Spol(k));
TEST(area_event_source_tables)
  // SourceProbe : ShearDislocation(mt, R3::XYZ(0,0,0)) exposes Entry(rt,k) as above, and
  //   Real PFrac() {return mWholeProbs[0].GetDiffProb(RAY_P);}
  for mt in {Tensor::SDR(125,40,90), Tensor::SDR(22.5,90,0)}:
    equal and area probes at degree 3, as above;
    rt in {RAY_P, RAY_SH, RAY_SV}, k with e.Entry > 1e-12*max: ratio == w_k within 1e-9;
    CHECK(std::fabs(e.PFrac() - 0.4) <= 1e-12 && std::fabs(a.PFrac() - 0.4) <= 1e-12);
    // icosahedral symmetry integrates these degree-4 patterns exactly
TEST(area_iso_mfp_exact)
  // writes build/toa_quadrature_iso.txt:
  //   "a kappa degree P_equal_relerr P_area_relerr S_equal_relerr S_area_relerr"
  //   (MFP rel err = exact/rate - 1, where 1/rate = GetMajorantMeanFreePath)
  for (a, kappa, checkDeg) in {(0.5, 0.3, 5), (2.0, 0.5, 6)}, degrees 5..7:
    sp = LopNorCrust(a, a, a, kappa); xP = ExactRate(sp, RAY_P, U, E, 400, 64); xS likewise;
    equal: UseTOADegree(d) -> Scatterer e(sp); area: AreaQuadratureScope(d) -> Scatterer s(sp);
    at d == checkDeg:
      CHECK(std::fabs(xP*s.GetMajorantMeanFreePath(RAY_P) - 1) <= 1e-4);   // and S
      if (a == 2) CHECK(std::fabs(xP*e.GetMajorantMeanFreePath(RAY_P) - 1) >= 3e-3);   // and S
TEST(area_iso_forced)
  // general_iso_forced's checks, inside AreaQuadratureScope q(5)
TEST(area_tables_thread_invariant)
  // iso_tables_thread_invariant's checks inside AreaQuadratureScope q(5). Factor its body
  // into a helper with no mode parameter; the scope selects the mode.
TEST(area_general_rate_matches_integration)
  // general_rate_matches_integration's 10 cases and 4-sigma tolerance, inside
  // AreaQuadratureScope q(5), with expected = 1/IntegratedRate(..., toas, &toas.SolidAngles())
TEST(area_general_distribution_matches_target)
  // general_distribution_matches_target's two states, inside AreaQuadratureScope q(5),
  // passing &toas.SolidAngles() to CheckDistributionMatchesTarget
TEST(area_axi_matches_general)
  // AreaQuadratureScope q(5); CheckAxiMatchesGeneral(LopNorCrust(2,2,0.5));
TEST(area_axi_bounds_valid)
  // AreaQuadratureScope q(5); axi_bounds_valid's P/S direction lists, specs {2,2,0.5} and
  // {0.5,0.5,2}, kappa 0.5, frequency 2 Hz;
  // CheckAxiBoundsValid(sp, q.toas, thetasP, thetaChiS, 5, &q.toas.SolidAngles())
TEST(isolated_area_overridemfp)
  // isolated_overridemfp's checks inside AreaQuadratureScope q(5), passing
  // &q.toas.SolidAngles() to CheckDistributionMatchesTarget
```

- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.**
  - Expect failures in `area_tables_weighted`, `area_event_source_tables`, `area_iso_mfp_exact`, `area_general_rate_matches_integration` and `area_general_distribution_matches_target`.
  - `area_iso_forced`, `area_tables_thread_invariant`, `area_axi_matches_general` and `area_axi_bounds_valid` may already pass. They guard invariants, so record which ones did.
- [ ] **Step 3: Implement.**
  - `PopulateProbDists()`: read the mode once, before the `ParallelFor`. In area mode, multiply each of the four GSATO values by `TOAWeight(k)` before `SetRelativeProb`. The equal branch stores exactly today's values. `m_spol[k]` is never weighted.
  - `ShearDislocation`: in area mode, multiply each squared amplitude by `TOAWeight(ctr)` before `SetRelativeProb`.
  - Update the comments listed under Files. The `ComputeMFPs()` comment must say that the ÷nTOA stays correct in both modes because Σw_k = N.
- [ ] **Step 4: Run `tests/test_anisotropy_unit.sh`** (it runs `isolated_*` tests in their own processes). Everything passes. Then run the full check and the default-mode reference check: every case prints PASS.
- [ ] **Step 5: Commit.**

```bash
git add scatterers.cpp events.cpp aniso_axi.hpp aniso_axi.cpp tests/anisotropy_tests.cpp
git commit -m "Weight scatterer and event-source tables by solid angle in area mode"
```

---

### Task 4: `--toa-quadrature` option, provenance lines, scenario tests

**Files:**
- Modify:
  - `cmdline.hpp`: add `OPT_TOAQUAD` after `OPT_TOA`, and `map["--toa-quadrature"] = OPT_TOAQUAD;` under "Universe";
  - `main.cpp`: a `process_option` case;
  - `model.hpp`: a `ModelParams` member and its default;
  - `model.cpp`: `ModelParams::Output()` and the `Model` constructor;
  - `scatterers.cpp`: the `PrintAllScatteringStats()` header.
- Create: `tests/test_toa_quadrature.sh`

**Interfaces:**
- Consumes: Task 2's `Set_TOA_Array(toa, quad)` and `Get_TOA_Quadrature()`; Task 1's three-argument `TesselSphere`.
- Produces:
  - `ModelParams::TOA_Quadrature`, of type `PhononSource::toa_quadrature_e`, default `TOA_EQUAL`;
  - the option `--toa-quadrature=equal|area`;
  - the three area-only stdout lines below.

- [ ] **Step 1: Write `tests/test_toa_quadrature.sh`.** Model it on `tests/test_aniso_scenarios.sh`: source `aniso_lib.sh`, use a `mktemp` work directory with a trap, and run `tests/run_aniso_case.sh` with `EXTRA_ARGS`. It fails, printing what failed, unless:
  - **Options:**
    - `./main --toa-quadrature=area --help` and `./main --toa-quadrature equal --help` exit 0;
    - `./main --toa-quadrature=bogus --help` and `./main --toa-quadrature --help` exit nonzero, and their output contains `Error processing command-line option`.
  - **Explicit equal equals the default:** Lop Nor with default args, `TOA_DEGREE=5`, 20000 phonons, 1 thread. A run with no flag and a run with `EXTRA_ARGS=--toa-quadrature=equal` pass `compare_runs ... strict`.
  - **Provenance:**
    - the area run's `stdout.txt` contains each of these lines exactly: `TOA Quadrature: area`, `|          (Area-weighted take-off-angle quadrature.)` and `#  TOA quadrature: area`;
    - `grep -ci quadrature` on the default run's `stdout.txt` prints 0.
  - **Effect:** the area run's normalized traces differ from the default run's.
  - **Determinism in area mode.** For Lop Nor default args, the anisotropic `lopnor_21`, and crust pinch `crustpinch_42` (layouts copied from `tests/test_aniso_scenarios.sh`):
    - two `THREADS=1` runs pass `compare_runs ... strict`;
    - a `THREADS=12` run's address-masked traces equal the `THREADS=1` traces.
    - Use 20000 phonons for every run.
  - **Seven- vs five-value equivalence in area mode:** the Lop Nor 15-value and 21-value isotropic layouts (copied from `tests/test_aniso_inputs.sh`), `TOA_DEGREE=4`, 20000 phonons, pass `compare_runs ... noargs`.
  - **Event test:** `./main --event-test --source=SDR,125,40,90` and the same with `--toa-quadrature=area` both exit 0. Their stdout is identical from the `@@ __EVENT_SOURCE_TEST__` line onward; the parameter echo above it legitimately gains the area line.
- [ ] **Step 2: Run it.** It fails at the first option check.
- [ ] **Step 3: Implement.**
  - **Option parsing:** `opt.PopValue_Text()` is `equal` or `area`. Anything else throws `Runtime("take-off-angle quadrature must be 'equal' or 'area'")`.
  - **`ModelParams::Output()`:** directly after the `TOA_Degree:` line, area mode only, print `TOA Quadrature: area`.
  - **`Model` constructor:**
    - build `new S2::TesselSphere(S2::TESS_ICO, par.TOA_Degree, par.TOA_Quadrature == PhononSource::TOA_AREA)`;
    - in area mode only, print `|          (Area-weighted take-off-angle quadrature.)` directly after the `(TesselSphere of degree N.)` line;
    - call `PhononSource::Set_TOA_Array(mpTOA_Array, par.TOA_Quadrature)`.
  - **`PrintAllScatteringStats()`:** in area mode only, print `#  TOA quadrature: area` directly after the `#  Overrides:` line.
- [ ] **Step 4: Run `tests/test_toa_quadrature.sh`,** then the full check (which now includes it), then the default-mode reference check: every case prints PASS.
- [ ] **Step 5: Commit.**

```bash
git add cmdline.hpp main.cpp model.hpp model.cpp scatterers.cpp tests/test_toa_quadrature.sh
git commit -m "Add the --toa-quadrature option with area-mode provenance output"
```

---

### Task 5: Scenario-level report and degree-9 build check

**Files:**
- Create:
  - `tests/compare_toa_quadrature.sh` (manual; not named `test_*`);
  - `docs/toa-quadrature-runs-2026-10/` (captured output: `compare-toa-quadrature.out`, `.err`, and a `README.md` holding the exact command, date, host, `git rev-parse --short HEAD` and wall time).

**Interfaces:**
- Consumes: `--toa-quadrature`; `tests/run_aniso_case.sh` (`EXTRA_ARGS`, `SEED`, `TOA_DEGREE`); `tests/aniso_lib.sh` (`scenario_model`, `auto_threads`).
- Produces: the capture that Task 6 quotes.

- [ ] **Step 1: Write the script.**
  - **Header comment:** purpose, environment variables, statistic and output, in the style of `tests/compare_axi_general.sh`.
  - **Environment, with defaults:**
    - `TQ_DEGREE=7`, `TQ_REPLICATES=4`, `TQ_THREADS=auto`, `TQ_SEED_BASE=596`;
    - `TQ_LOPNOR_PHONONS=200000`, `TQ_CRUSTPINCH_PHONONS=50000`;
    - `TQ_CONFIGS="lopnor:iso lopnor:axi4 crustpinch:iso crustpinch:axi4"`, where `axiR` means `scenario_model SCENARIO 0.01 axi R` and `iso` means `scenario_model SCENARIO 0.01 iso`;
    - `TQ_BUILD_CHECK=1`.
  - **Runs per config:**
    - replicates k = 1..K run in equal mode, with seed `TQ_SEED_BASE + k*2^32` (computed with `awk 'BEGIN{printf "%.0f", ...}'`, exact below 2^53);
    - replicates k = K+1..2K run in area mode (`EXTRA_ARGS=--toa-quadrature=area`).
  - **Energy table, per config:**
    - E_s is a station's total energy: the sum of the X+Y+Z columns over every trace bin. Parse it as `extract_energy` in `tests/compare_axi_general.sh` does, but without windows.
    - For each station, compute the per-mode mean and standard error, the ratio r_s = mean_area/mean_equal, and SE(r_s) by error propagation.
    - Print one Markdown row per config: `config | stations | min r | median r | max r | stations with |r-1| > 3 SE`.
  - **MFP table, per config:** from the first run of each mode, print one row per scatterer-dump data row. Rows lie between `BEGIN SCATTERER DUMP` and `END SCATTERERS`, excluding those that start with `#`; match rows by order. Columns: `a ay az mode | MFP(P) equal area %change | MFP(S) equal area %change`. Isotropic rows have no ay, az or mode; print `-`.
  - **Build check** (`TQ_BUILD_CHECK=1`):
    - degree 9, `--num-phonons=0`, Lop Nor `iso` and `axi10`, equal against area, 3 runs each;
    - measure median wall seconds and maximum resident set size with `/usr/bin/time -l` on macOS (parse `maximum resident set size`, in bytes), falling back to `-v` on Linux (parse `Maximum resident set size (kbytes)`);
    - print a Markdown table.
  - **Streams:** progress goes to stderr. The script exits nonzero if any run fails or lacks `__SIMULATION_COMPLETE__`.
- [ ] **Step 2: Smoke-test** with `TQ_CONFIGS=lopnor:iso TQ_REPLICATES=2 TQ_LOPNOR_PHONONS=20000 TQ_BUILD_CHECK=0 TQ_DEGREE=5 tests/compare_toa_quadrature.sh`. It prints both tables and exits 0.
- [ ] **Step 3: Run the default configuration** with stdout to `docs/toa-quadrature-runs-2026-10/compare-toa-quadrature.out` and stderr to `.err`. Write `README.md` with the provenance fields.
- [ ] **Step 4: Copy the measurement files** written by Tasks 1–3's tests (`build/toa_quadrature_solid_angles.txt`, `build/toa_quadrature_report.txt`, `build/toa_quadrature_pframe.txt`, `build/toa_quadrature_iso.txt`) into `docs/toa-quadrature-runs-2026-10/`, after a fresh `tests/test_anisotropy_unit.sh` run. List them in the README.
- [ ] **Step 5: Commit.**

```bash
git add tests/compare_toa_quadrature.sh docs/toa-quadrature-runs-2026-10
git commit -m "Measure how area-weighted quadrature moves scenario outputs"
```

---

### Task 6: Documentation

**Files:**
- Modify: `doc/MANUAL.md`, `docs/toa-quadrature-bias.md`, `tests/anisotropy_tests.cpp` (the S90SH reference comment in `report_directional_mfps` only: point to `area_report_directional_mfps`)

**Interfaces:**
- Consumes: `docs/toa-quadrature-runs-2026-10/*` (every number must come from there, from a test assertion, or from the spec's spike table, labelled as such); the spec's decisions 7–9.

- [ ] **Step 1: `doc/MANUAL.md`.**
  - **Option entry.** Add `* <tt>--toa-quadrature=_MODE_</tt>` directly after the `--toa-degree` entry, covering:
    - `equal` (default): every direction carries the same weight, as in all earlier versions;
    - `area`: each direction is weighted by its triangle's solid angle, w = Ω·N/4π;
    - area mode changes every output, isotropic regions included, and prints the three provenance lines;
    - `--event-test` ignores the option.
  - **Validity limits.** Rewrite the quadrature paragraph so it keeps every existing equal-mode fact and number, and adds:
    - area mode removes the fixed offset (cite `toa_quadrature_iso.txt`);
    - anisotropic report values converge to the exact integral, e.g. S90SH 374.7 km;
    - the P frame dependence becomes a residual that shrinks about 4× per degree (cite `toa_quadrature_pframe.txt` at degree 6);
    - equal weights also distort the sampled direction density by the factor mean(Ω)/Ω_k, between 0.829 and 1.078 (from `toa_quadrature_solid_angles.txt`), and area mode removes this;
    - replace "it has not been implemented" with the option name.
  - **Scatterer dump.** Add one sentence: the dump's values follow the selected quadrature, and area mode adds the `#  TOA quadrature: area` header line.
- [ ] **Step 2: `docs/toa-quadrature-bias.md`.** Restructure into these sections:
  - Status (opt-in since branch `quadrature`; default unchanged; ruling R9 history);
  - Mechanism (existing text, plus the direction-density effect);
  - Isotropic effect (equal against area against exact, per degree);
  - Anisotropic effect (keep the table, and name `area_report_directional_mfps` as its in-repo reproduction);
  - Frame orientation for P incidence (add the measured area residual);
  - Event source (totals exact in both modes, direction density changes);
  - Implementation (replaces "Affected code": what changed in each file, and what inherits the weights unchanged);
  - Consequences (opt-in; the forced-isotropic and seven- vs five-value equivalences hold in both modes; provenance lines);
  - Validation (each list item of the old section, with the test or script that checks it and its result);
  - Scenario-level effect (tables from `compare-toa-quadrature.out`, including the degree-9 build cost);
  - Provenance.
- [ ] **Step 3: Verify.** `make && tests/test_anisotropy_unit.sh` passes. `git diff --stat` shows only the three files. Every number in the new text appears in a file under `docs/toa-quadrature-runs-2026-10/`, in a test assertion, or in the spec's spike table; list each number and its source in the task report.
- [ ] **Step 4: Commit.**

```bash
git add doc/MANUAL.md docs/toa-quadrature-bias.md tests/anisotropy_tests.cpp
git commit -m "Document the area-weighted take-off-angle quadrature option"
```
