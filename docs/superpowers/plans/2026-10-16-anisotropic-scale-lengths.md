# Anisotropic Heterogeneity Scale Lengths Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let every region's heterogeneity spectrum carry three von Kármán scale lengths `ax, ay, az`
(isotropic, axisymmetric, or anisotropic). Scattering is sampled exactly from the anisotropic
spectrum, and isotropic runs stay byte-identical.

**Architecture:**
- `HetSpec` and `ScatterParams` carry `ax, ay, az`. An anisotropic scatterer builds today's
  isotropic tables at `a_min` and uses them as the proposal for delta tracking (null
  collisions) in `Phonon::Propagate()`.
- Triaxial scatterers use the general sampler: an isotropic proposal plus an exact
  accept/reject step.
- Axisymmetric scatterers add an incidence-tailored coarse proposal (`AxiProposal`). It relies
  on the hierarchical ordering of the take-off-angle set and on rigorous bounds.
- Isotropic scatterers never enter the new code.

**Tech Stack:** C++11, GNU Make, POSIX shell tests. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-10-16-anisotropic-scale-lengths-design.md`

## Global Constraints

- Isotropic runs keep today's arithmetic and random-draw order. Fixed-seed non-image outputs
  (`seis_traces_asc.dat`, `seis_*.octv`, `out_mparams.octv`, stdout with addresses and the
  `BUILD STATS` line masked) are byte-identical to the baseline build of `aaec540`.
- Random numbers come only from `RandomUniform()` and `RandomInteger()` in `random.hpp`.
  Propagation writes no shared state. Scatterers are read-only once constructed.
- Output changes only append columns or header text. Runs whose regions are all isotropic
  print exactly what they print today.
- The fabric frame is local (East, North, Up) at the scattering point, from the global `ECS`
  (`GetEast`, `GetNorth`, `GetUp`). That's internal `+x, +y, +z` in `ENU_ORTHO`/`RAE_ORTHO`.
- No tabulation or interpolation of anisotropic rates or patterns. Every bound is a rigorous
  upper bound. A violated bound throws `std::logic_error` and is never clamped.
- Memory: isotropic scatterers are unchanged. Each axisymmetric scatterer adds at most about
  9.2 MB, for its bin tables stored as `double`.
- Do not edit the benchmark numbers in `README.md`, `doc/MANUAL.md`, or
  `docs/multithreaded-simulation-handoff.md`.
- Do not change `tests/run_lopnor_thread_case.sh` or
  `tests/test_lopnor_thread_equivalence.sh`.
- **Full check**, run at the end of every task:
  `make && for t in tests/test_*.sh; do "$t" || exit 1; done`. To skip rebuilding the
  baseline each time, build it once and export its path:
  `git worktree add --detach /tmp/r3d-baseline aaec540 && make -C /tmp/r3d-baseline && export BASELINE_EXE=/tmp/r3d-baseline/main`.

## Review Focus

1. **Vertical incidence** (incoming ray parallel to local Up). The local θ̂ and the
   polarization angle χ are undefined there. Results must be finite (no NaN) and use a
   consistent frame. Tests: `general_vertical_incidence` (Task 4), `axi_vertical_incidence`
   (Task 7).
2. **Curved mappings** (`RAE_CURVED`, used by crust pinch). The fabric must follow the local
   vertical, and the P canonical frame must use the local θ̂, not the lab θ̂. Tests:
   `general_curved_frame` (Task 4), `axi_curved_frame` (Task 7).
3. **`--overridemfp` and `--nodeflect` with anisotropic scatterers.** The override fixes the
   path length but the pattern stays anisotropic. No-deflect keeps the anisotropic event
   rate but doesn't deflect. Tests: `isolated_overridemfp`, `isolated_nodeflect` (Task 4).
4. **`eps = 0` anisotropic regions.** There's no scattering, so the mean free path must be
   infinite, with no NaN from zero majorants and no division by zero in thinning. Tests:
   `general_eps_zero` (Task 4), `axi_eps_zero` (Task 7).
5. **Low take-off-angle degrees** (below 3, common in quick runs). Coarse cells become single
   fine directions, and bounds and sampling must still work. Test: `axi_low_degree` (Task 7).

---

### Task 1: Isotropic regression harness

**Files:**
- Create: `tests/aniso_lib.sh`, `tests/run_aniso_case.sh`,
  `tests/test_aniso_isotropic_regression.sh`

**Interfaces:**
- Produces:
  - `tests/run_aniso_case.sh SCENARIO THREADS OUTDIR [PHONONS] [MODEL_ARGS]`
    - `SCENARIO` is `lopnor` or `crustpinch`. It runs the simulation-only part of
      `do-lopnor.sh` or `do-crustpinch.sh`.
    - `PHONONS` defaults to 200000 (`lopnor`) or 50000 (`crustpinch`).
    - `MODEL_ARGS` defaults to the do-script values.
    - Environment variables: `R3D_EXE` (default `$ROOT/main`), `TOA_DEGREE` (5), `SEED` (596),
      `REPORTS` (`ALL_OFF`), and `EXTRA_ARGS` (appended verbatim, e.g. `--dump-grid`).
    - Writes `stdout.txt`, `stderr.txt`, `reports.dat`, `out_mparams.octv`,
      `seis_traces_asc.dat`, and `seis_*.octv` into `OUTDIR`.
  - Functions in `tests/aniso_lib.sh`:
    - `normalize_stdout FILE`
    - `normalize_stdout_noargs FILE`
    - `compare_runs DIR_A DIR_B MODE`, where `MODE` is `strict` or `noargs`
    - `trace_has_signal FILE`
    - `require_sct_bands SCENARIO REPORT BAND...`

- [ ] **Step 1: Write `tests/run_aniso_case.sh`.** Use the arguments from the do-scripts:

| | `lopnor` | `crustpinch` |
|---|---|---|
| `--grid-compiled` | `1` | `5` |
| `--source` | `SDR,125,40,90,0.0` | `SDR,22.5,90,0` |
| `--source-loc` | `425.54,-169.53,-31.02` | `0,0,-10` |
| other options | `--range=1200 --flatten` | (none) |
| model args | `0.8,0.01,0.5,0.2,50,0.8,0.01,0.5,0.3,1000,0.8,0.01,0.7,0.5,300` | `SCAT1..SCAT5,LAYERS,PINCH` from `do-crustpinch.sh` (32 values) |
| seismometers | `--seis-p2p="$XINXYZ,$WUSXYZ,1.0,2.0,40.0,8"` and `--seis-p2p="$XINXYZ,$MAKXYZ,1.0,2.0,40.0,8"` (coordinates as in `run_lopnor_thread_case.sh`) | three `--seis-p2p=0,AZ,0,950,AZ,0,1.0,2.0,40.0,8` for `AZ` = 67.5, 112.5, 90 |

Both scenarios also get `--frequency=2.0 --timetolive=600 --binsize=2.0`, `--seed`, and
`--threads`. Run from inside `OUTDIR`, as `run_lopnor_thread_case.sh` does.

- [ ] **Step 2: Write `tests/aniso_lib.sh`.**
  - `normalize_stdout`: rewrite `0x[0-9a-f]+` to `0xADDR` and delete the `BUILD STATS` line.
  - `normalize_stdout_noargs`: as above, then also delete the line after
    `Args provided for "compiled" model are:`, and rewrite `with [0-9]+ args` to
    `with N args`.
  - `compare_runs`:
    - always `cmp` `seis_traces_asc.dat` (addresses masked as in
      `test_lopnor_thread_equivalence.sh`) and every `seis_*.octv`;
    - `strict` also `cmp`s `out_mparams.octv` and the `normalize_stdout` output;
    - `noargs` skips `out_mparams.octv`, which echoes the arguments, and compares
      `normalize_stdout_noargs` output.
  - `trace_has_signal`: copy the awk from `test_lopnor_thread_equivalence.sh`.
  - `require_sct_bands`: count `SCT:` lines whose location falls in each named band, and fail
    with the band name if any count is 0.
    - `lopnor` uses internal `z` (the 3rd number after `xyz:(`).
    - `crustpinch` (curved) uses `depth = 6371 − |xyz − (0,0,−6371)|` and
      `range = 6371·atan2(√(x²+y²), z+6371)`.
    - Starting bands, to be confirmed inside each region against `--dump-grid`:

| Scenario | Band | Location range |
|---|---|---|
| `lopnor` | `sediments` | z > 1.0 |
| `lopnor` | `crust` | −30 < z < −15 |
| `lopnor` | `mantle` | z < −60 |
| `crustpinch` | `sediments` | depth < 1.5 |
| `crustpinch` | `crust` | 12 < depth < 20, range < 300 or range > 550 |
| `crustpinch` | `pinched` | 12 < depth < 20, 380 < range < 460 |
| `crustpinch` | `moho` | 33 < depth < 36 |
| `crustpinch` | `mantle` | depth > 50 |

- [ ] **Step 3: Write `tests/test_aniso_isotropic_regression.sh`.**
  - Build the current tree with `make`.
  - If `BASELINE_EXE` is unset, run
    `git worktree add --detach "$work/baseline" "${BASELINE_REV:-aaec540}"` and
    `make -C "$work/baseline"`, then remove the worktree in the `trap`.
  - For each scenario, run the baseline and the current build with `THREADS=8` and default
    phonons, then `compare_runs … strict`, and require `trace_has_signal`.
  - Then run the current build with `REPORTS=SCT` and 20000 phonons on 1 thread, and call
    `require_sct_bands` with every band for the scenario.
  - Finish with `echo "isotropic regression checks passed"`.

- [ ] **Step 4: Run it.** `tests/test_aniso_isotropic_regression.sh` should print
  `isotropic regression checks passed`. The code is unchanged, so any difference points to a
  harness bug. If a band count is 0, adjust that band inside its region using
  `EXTRA_ARGS=--dump-grid`, or raise the coverage phonon count. Never widen a band across a
  region boundary.

- [ ] **Step 5: Run the full check, then commit.**

```bash
git add tests/aniso_lib.sh tests/run_aniso_case.sh tests/test_aniso_isotropic_regression.sh
git commit -m "Add isotropic regression harness for anisotropy work"
```

---

### Task 2: Three scale lengths in `HetSpec` and `ScatterParams`, plus the unit-check harness

**Files:**
- Modify: `elastic.hpp` (class `HetSpec`, new schema `HSnexyzk`), `scatparams.hpp`,
  `scatparams.cpp`, `Makefile`
- Create: `tests/anisotropy_tests.cpp`, `tests/test_anisotropy_unit.sh`

**Interfaces:**
- Produces, in `elastic.hpp`, namespace `Elastic`:
  - `HetSpec`: protected `HetSpec(Real nu, Real eps, Real ax, Real ay, Real az, Real k)` throws
    `std::invalid_argument("HetSpec: scale lengths ax, ay, az must be positive")` unless all
    three are greater than 0. The existing four-argument constructor sets
    `ax = ay = az = a`.
  - Public accessors: `Real a()` (equals `ax()`), `Real ax()`, `Real ay()`, `Real az()`,
    `bool IsIsotropic()` (`ax==ay && ay==az`), `bool IsAxisymmetric()` (`ax==ay && ay!=az`).
  - `class HSnexyzk : public HetSpec { HSnexyzk(Real nu, Real eps, Real ax, Real ay, Real az, Real k); }`
- Produces, in `ScatterParams`:
  - New members `ax, ay, az`; the existing member `a` becomes `min(ax, ay, az)`.
  - Accessors: `Real GetAx()`, `GetAy()`, `GetAz()`, `bool IsIsotropic()`,
    `bool IsAxisymmetric()`, `Real GetKP()` (`el/gam0`), `Real GetKS()` (`el`).
  - `Real AnisoMajorantFactor() const` returns `M = ax·ay·az/a³`.
  - `Real PSATO_Aniso(const R3::XYZ & m) const` implements Sato19 Eq. 3, with `m` in fabric
    coordinates.
  - `Real AnisoRatio(const R3::XYZ & m) const` returns `PSATO_Aniso(m)/PSATO(|m|)`.
  - `Real AnisoAcceptance(const R3::XYZ & m) const` returns
    `[(1+a²|m|²)/(1+ax²mx²+ay²my²+az²mz²)]^(κ+3/2)`.
  - `CompareRoughly` adds `dax²+day²+daz²`.
- Produces, in the harness: `build/anisotropy_tests [FILTER | --list]`.
  - `--list` prints every test name, one per line.
  - With no argument it runs every test whose name doesn't start with `isolated_`.
  - With an argument it runs tests whose name contains `FILTER`; `isolated_*` tests run only
    when named exactly.
  - It prints `PASS name` or `FAIL name: message`, and exits 1 on any failure.
  - Helpers:
    - `UseTOADegree(int)` caches one `S2::TesselSphere` per degree and calls
      `PhononSource::Set_TOA_Array`.
    - `ScatterParams LopNorCrust(Real ax, Real ay, Real az, Real kappa=0.5, Real eps=0.01)`
      calls `ScatterParams::SetFrequencyHertz(2.0)` and returns
      `ScatterParams(Elastic::VpVs(6.2,3.6), Elastic::HSnexyzk(0.8,eps,ax,ay,az,kappa))`.
    - `R3::XYZ Dir(Real theta_deg, Real phi_deg)`.

- [ ] **Step 1: Add the harness and the failing tests.**
  - Add the `Makefile` target:

```make
test_objects := $(filter-out $(objdir)/main.o,$(objects))
$(objdir)/anisotropy_tests : tests/anisotropy_tests.cpp $(test_objects) $(scatterers_hpp) $(ecs_hpp) $(phonons_hpp)
	$(CPP) tests/anisotropy_tests.cpp $(test_objects) -I. $(FLAGS) -o $@
.PHONY : anisotropy-tests
anisotropy-tests : directories $(objdir)/anisotropy_tests
```

  - `tests/test_anisotropy_unit.sh` runs `make`, `make anisotropy-tests`, then
    `build/anisotropy_tests`, then each `isolated_*` name from `--list` in its own process,
    and prints `anisotropy unit checks passed`. Separate processes are needed because the
    override flags are class-static and can't be reset.
  - Tests:

```cpp
test_hetspec_modes:   HSneak(0.8,0.01,0.5,0.3): IsIsotropic, !IsAxisymmetric, ax()==ay()==az()==a()==0.5
                      HSnexyzk(0.8,0.01,2,2,0.5,0.3): IsAxisymmetric; HSnexyzk(...,2,1,0.5,...): neither
                      HSnexyzk(...,0.5,0.5,0.5,...): IsIsotropic and a()==0.5
test_hetspec_validation: HSnexyzk with ax=0, ay=-1, az=0 each throw std::invalid_argument
test_scatparams_iso_unchanged: ScatterParams from HSneak(0.8,0.01,0.5,0.3) vs HSnexyzk(...,0.5,0.5,0.5,...):
                      CompareRoughly()==0, GetA()==0.5, and GSATO outputs bitwise equal at 20 TOAs (degree 2)
test_scatparams_aniso_reference: LopNorCrust(2,2,0.5): GetA()==0.5, GetAx()==2, GetAz()==0.5, AnisoMajorantFactor()==16
test_psato_aniso_iso_limit: isotropic a=0.5: PSATO_Aniso(m)==PSATO(|m|) within 1e-12 rel for 50 random m, |m|<=8
test_aniso_ratio_bounds: LopNorCrust(2,2,0.5,0.3): AnisoRatio(0)==16 (1e-12); AnisoRatio((0,0,m))==16 for m in {0.5,1,4};
                      AnisoRatio(m)<=16*(1+1e-12) and AnisoAcceptance(m)==AnisoRatio(m)/16 (1e-12) for 1000 random m
test_matching_distinguishes: CompareRoughly(2,2,0.5 vs 2,0.5,2)>0, (2,2,0.5 vs iso 0.5)>0, identical aniso ==0
```

- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.** Expect a compile failure: `HSnexyzk`,
  `GetAx`, and the other new names are undefined.

- [ ] **Step 3: Implement the interfaces above** in `elastic.hpp` (add `#include <stdexcept>`),
  `scatparams.hpp`, and `scatparams.cpp`.
  - The two `ScatterParams` constructors and the three special test constructors set
    `ax, ay, az` from the spec (or all to `a`), with `a = min`.
  - `PSATO()` is untouched. For isotropic specs `a` is unchanged, so `GSATO` is bitwise
    identical.

- [ ] **Step 4: Run `tests/test_anisotropy_unit.sh`.** Expect every test to `PASS`, followed by
  `anisotropy unit checks passed`.

- [ ] **Step 5: Run the full check, then commit.** The isotropic regression must pass.

```bash
git add elastic.hpp scatparams.hpp scatparams.cpp Makefile tests/anisotropy_tests.cpp tests/test_anisotropy_unit.sh
git commit -m "Carry three heterogeneity scale lengths through HetSpec and ScatterParams"
```

---

### Task 3: Seven-value input layouts, grid-dump columns, and the sampling guard

**Files:**
- Modify: `user.cpp` (model-index 1 dispatch), `user_LopNorCyl_inc.cpp`,
  `user_LopNorCylMoho_inc.cpp`, `user_NSCP_inc.cpp`, `grid.hpp`, `grid.cpp`,
  `scatterers.cpp` (`GetRandomScatteredRelativePhonon`), `doc/MANUAL.md`, `do-lopnor.sh`,
  `do-crustpinch.sh`
- Create: `tests/test_aniso_inputs.sh`
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Consumes: `HSnexyzk` and `HetSpec::IsIsotropic` (Task 2), and `run_aniso_case.sh` and
  `aniso_lib.sh` (Task 1).
- Produces:
  - `void GridNode::OutputAsAscii(std::ostream & out, const std::string & prefix, bool aniso_columns = false) const`.
  - Seven-value layouts, each region as `nu, eps, ax, ay, az, kappa, Q`:

| Model | Count | Regions |
|---|---:|---|
| Lop Nor baseline | 21 | sediments, crust, mantle |
| Lop Nor Moho | 28 | sediments, crust, mantle, Moho |
| Crust pinch | 35 | sediments, crust, pinched crust, Moho, mantle |
| Crust pinch | 38 | 35, then `SediThick, CrustThick, MohoThick` |
| Crust pinch | 42 | 38, then `PinchFrac, DepthFrac, SediFrac, MohoFrac` |

  - Model index 1 selects the Moho model only for counts 20 and 28.
  - `Scatterer::GetRandomScatteredRelativePhonon` throws
    `std::logic_error("anisotropic scatterers must be sampled through AttemptAnisotropicScatter")`
    when `!mParams.IsIsotropic()`.

- [ ] **Step 1: Write the failing tests.**
  - `tests/test_aniso_inputs.sh`, using phonon count 20000 (`lopnor`) or 5000 (`crustpinch`),
    `TOA_DEGREE=4`, and `compare_runs … noargs`:
    - **Equivalence.** Each five-value layout and its seven-value twin with
      `ax=ay=az=a` must be identical, for counts 15→21, 20→28, 25→35, 28→38, and 32→42.
    - **Dispatch.** Stdout must contain `Selected Lop Nor Baseline Model` for 15 and 21, and
      `Selected Lop Nor Moho Model` for 20 and 28.
    - **Grid dump.** Run the anisotropic Lop Nor layout with 0 phonons and
      `EXTRA_ARGS=--dump-grid`:
      `0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,2.0,2.0,0.5,0.3,1000,0.8,0.01,1.4,0.7,0.35,0.5,300`.
      - The header contains `nu eps a k ay az`.
      - Every node row has 17 fields.
      - Rows with field 15 (`kappa`) `0.30000` (crust) have fields 14, 16, 17 =
        `2.00000 2.00000 0.50000`.
      - Rows with field 15 `0.50000` (mantle) have `1.40000 0.70000 0.35000`.
      - The default layout's dump has 15 fields per row.
    - **Bad counts.** 22 values exits nonzero with `wrong number of model args`.
    - **Bad scale length.** The 21-value layout with crust `ay=0` exits nonzero, and
      `stdout.txt` or `stderr.txt` contains `scale lengths ax, ay, az must be positive`.
  - In the unit harness, `test_scatterer_guard`: `Scatterer(LopNorCrust(2,2,0.5))` at degree 2;
    `GetRandomScatteredRelativePhonon(RAY_P)` throws `std::logic_error`.

- [ ] **Step 2: Run `tests/test_aniso_inputs.sh`.** Expect it to fail at the 21-value
  equivalence check (wrong number of model args).

- [ ] **Step 3: Implement.**
  - In each model, existing layouts set `Ay = Az = A` for every region. All `HSneak` calls in
    the three model files become `HSnexyzk(Nu, Eps, A, Ay, Az, K)`.
  - The new cases parse `Q` through `QmQk` as today. Crust pinch 42 falls through to 38 and
    then 35, matching the 32 → 28 → 25 pattern.
  - `Grid::DumpGridToAscii` scans node data for any `!getHS().IsIsotropic()`. If it finds one,
    the header line becomes `#  i j k    x y x    vp vs rho qp qs    nu eps a k ay az`. That
    is the same line count, so the do-scripts' `tail -n +9` still works. It then passes
    `aniso_columns=true`, so every row appends `ay az` (`***` rows append two more `***`).
  - Document in `doc/MANUAL.md`, in the heterogeneity section:
    - the three modes and their auto-detection;
    - the fabric frame: x=East, y=North, z=Up at the scattering point, following Earth
      curvature in curved mappings;
    - the layout table above;
    - that every region uses the same layout width;
    - the validity notes `k²·a_max²·ε² ≪ 1` and Margerin's very-flat-scatterer caution.
  - Add one comment line under the `SCAT` lines in each do-script showing the seven-value form.

- [ ] **Step 4: Run `tests/test_aniso_inputs.sh` and `tests/test_anisotropy_unit.sh`.** Both
  must pass.

- [ ] **Step 5: Run the full check, then commit.**

```bash
git add user.cpp user_LopNorCyl_inc.cpp user_LopNorCylMoho_inc.cpp user_NSCP_inc.cpp grid.hpp grid.cpp scatterers.cpp doc/MANUAL.md do-lopnor.sh do-crustpinch.sh tests/test_aniso_inputs.sh tests/anisotropy_tests.cpp
git commit -m "Add seven-value anisotropic model-args layouts and grid-dump columns"
```

---

### Task 4: General anisotropic sampler

**Files:**
- Modify: `scatterers.hpp`, `scatterers.cpp`
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Consumes: the `ScatterParams` anisotropic methods (Task 2), the global `ECS`, and
  `PhononSource::pTOA`, `mPDists`, `mWholeProbs`.
- Produces, in class `Scatterer`:

```cpp
enum sampler_e { SAMPLER_AUTO, SAMPLER_GENERAL, SAMPLER_AXISYMMETRIC };
Scatterer(ScatterParams par, sampler_e choice = SAMPLER_AUTO);
    // AUTO: isotropic spec -> today's path; axisymmetric -> AXISYMMETRIC (Task 7; GENERAL until then); else GENERAL.
bool IsAnisotropic() const;                       // true unless the isotropic path is active
Real GetMajorantMeanFreePath(raytype rt) const;   // mMeanFreeP[rt]
struct AnisoAttempt {
  bool IsReal;  Real InPol;  Phonon Rel;  int Conv;  Index Toa;  int Draws;
  AnisoAttempt();   // IsReal=false, InPol=0, Rel(S2::ThetaPhi(0,0),RAY_P), Conv=-1, Toa=0, Draws=0
};
AnisoAttempt AttemptAnisotropicScatter(raytype intype, const R3::XYZ & loc,
                                       const S2::ThetaPhi & dir, Real pol);
Real IntegratedRate(raytype intype, const R3::XYZ & dir, const R3::XYZ & pol,
                    const S2::S2Set & dirs) const;   // fabric-coordinate dir/pol; mean over relative dirs of
                                                     // the in-type's two conversions, GSATO(a_min) * AnisoRatio
static R3::XYZ ToFabric(const R3::XYZ & loc, const R3::XYZ & v);          // (v·E, v·N, v·U) at loc
static Real RotationForS1(const S2::ThetaPhi & dir, const R3::XYZ & s1);  // rot so OrthoAxes(dir,rot).S1()
                                                     // is along s1 projected normal to dir; 0 if that projection is ~0
```

  - Under the general sampler, `mMeanFreeP[in]` is the isotropic mean free path at `a_min`
    divided by `M`.
  - Conversion wavenumbers `(k_in, k_out)`: GPP `(kP,kP)`, GPS `(kP,kS)`, GSP `(kS,kP)`,
    GSS `(kS,kS)`.

- [ ] **Step 1: Write the failing tests.** Use degree 5 unless noted. A helper
  `MeanDistanceToRealEvent(s, in, dir, pol, n)` sums `GetRandomPathLength(in)` draws until
  `AttemptAnisotropicScatter(in, origin, dir, RotationForS1(dir,pol))` returns `IsReal`, and
  averages over `n` events.

```cpp
test_general_iso_forced: Scatterer(LopNorCrust(.5,.5,.5), SAMPLER_GENERAL): IsAnisotropic();
    GetMajorantMeanFreePath equals Scatterer(same,AUTO) MFP (1e-12); 2000 attempts all IsReal
test_general_rate_matches_integration (4a): LopNorCrust(2,2,0.5); in-states: P along E,N,U and Dir(60,30);
    S along E (pol N, pol U), N (pol E, pol U), U (pol E, pol N). n=20000:
    |mean - 1/g0| < 4*(1/g0)/sqrt(n), g0 = IntegratedRate(in, dir, pol, own TOA set)
test_general_distribution_matches_target (4b): P Dir(60,30) and S along E pol U; 200000 accepted draws binned by
    (Conv, Toa / 4^(5-2)); expected from GSATO*AnisoRatio summed per bin: chi2/dof < 1.5 over bins with expected >= 20
test_general_symmetry (4d): IntegratedRate equal within 1e-9 rel for: axisymmetric (2,2,0.5) P along E vs N;
    P Dir(60,30) vs its x-mirror; triaxial (2,1,0.5) P along E vs -E
test_general_vertical_incidence: S along U, pol E: 1000 attempts give finite fields, IntegratedRate finite and > 0
test_general_curved_frame: ECS.SetMapping(EarthCoords::RAE_CURVED); at loc 950 km along x on the surface,
    ToFabric(loc, ECS.GetUp(loc)) == (0,0,1) within 1e-12 and ECS.GetUp(loc).z() < cos(8 deg); restore ENU_ORTHO
test_general_eps_zero: LopNorCrust(2,2,0.5,0.5,eps=0): GetMajorantMeanFreePath(RAY_P) is +inf; construction does not throw
isolated_overridemfp: Scatterer::OverrideMFP(100,50) then LopNorCrust(2,2,0.5): mean GetRandomPathLength(RAY_P)
    over 20000 = 100 within 3%; every AttemptAnisotropicScatter IsReal; repeat 4b for P Dir(60,30): chi2/dof < 1.5
isolated_nodeflect: Scatterer::SetNoDeflect() then LopNorCrust(2,2,0.5): real attempts return Rel direction theta==0,
    Rel raytype == intype, InPol == incoming pol
```

- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.** Expect a compile failure on
  `AttemptAnisotropicScatter`, `IntegratedRate`, `SAMPLER_GENERAL`, and the other new names.

- [ ] **Step 3: Implement.** In the constructor, an anisotropic sampler runs `ComputeMFPs()`
  first and then divides by `M`, unless the override is set. `AttemptAnisotropicScatter` for
  the general sampler:
  1. Build `AA = R3::OrthoAxes(dir.Theta(), dir.Phi(), pol)`.
  2. Draw `Conv` from `mWholeProbs[intype].GetRandomIndex()`, then `Toa` from
     `mPDists[Conv].GetRandomIndex()`.
  3. Compute `ê_out = AA.Express((*pTOA)[Toa].XYZ())` and
     `m = k_out·ê_out − k_in·dir.XYZ()`.
  4. Accept if `RandomUniform() < mParams.AnisoAcceptance(ToFabric(loc, m))`.
  5. Build `Rel` and its polarization by today's rules (`m_spol` for GSS, otherwise 0), and
     set `InPol = pol` and `Draws = 1`.

  With `cm_MFPOverride_b`, loop steps 2–4 until accepted and count `Draws`. With
  `cm_NoDeflect_b`, an accepted `Rel` is `Phonon(S2::ThetaPhi(0,0), intype)` with polarization
  0. `ToFabric` calls `ECS.GetEast`, `GetNorth`, and `GetUp` at `loc`.

- [ ] **Step 4: Run `tests/test_anisotropy_unit.sh`.** Every test must pass, including the
  `isolated_*` invocations.

- [ ] **Step 5: Run the full check, then commit.**

```bash
git add scatterers.hpp scatterers.cpp tests/anisotropy_tests.cpp
git commit -m "Add exact general sampler for anisotropic scatterers"
```

---

### Task 5: Null-collision propagation and end-to-end anisotropic scenarios

**Files:**
- Modify: `phonons.cpp` (`Phonon::Propagate`)
- Create: `tests/test_aniso_scenarios.sh`

**Interfaces:**
- Consumes: `Scatterer::IsAnisotropic`, `GetRandomPathLength`, `AttemptAnisotropicScatter`,
  and `AnisoAttempt` (Task 4).
- Produces: anisotropic runs end to end. Later tasks don't call new code here.

- [ ] **Step 1: Write `tests/test_aniso_scenarios.sh`.** Use the anisotropic layouts below,
  default phonons, and `TOA_DEGREE=5`.

| Scenario | Layout |
|---|---|
| `lopnor` (21 values) | `0.8,0.01,0.5,0.5,0.5,0.2,50,0.8,0.01,2.0,2.0,0.5,0.3,1000,0.8,0.01,1.4,0.7,0.35,0.5,300` |
| `crustpinch` (42 values) | `0.8,0.01,0.2,0.2,0.2,0.2,200,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.8,0.8,0.2,0.3,1500,0.8,0.01,0.6,0.3,0.15,0.4,1500,0.8,0.01,0.2,0.2,0.2,0.5,900,2.0,30.0,5.0,.3666667,.4736842,1,1` |

  For each scenario:
  - stdout contains `__SIMULATION_COMPLETE__`;
  - two `THREADS=1` runs pass `compare_runs … strict`;
  - the `THREADS=12` traces are identical to the `THREADS=1` traces, after address masking;
  - `trace_has_signal` holds;
  - the anisotropic traces differ from the default isotropic run's (`! cmp -s`);
  - `REPORTS=SCT` with 20000 phonons passes `require_sct_bands` for the anisotropic regions
    (`lopnor`: `crust mantle`; `crustpinch`: `crust pinched moho`).

  The script ends with `echo "anisotropic scenario checks passed"`.

- [ ] **Step 2: Run it.** Expect failure: the first scattering event throws the Task 3 guard
  (`std::logic_error`), which surfaces as `** Error during simulation execution`.

- [ ] **Step 3: Implement the anisotropic branch in `Propagate()`.** Keep the isotropic
  branch's statements and their order exactly as they are now.

```cpp
TravelRec travel = mpCell->GetPathToBoundary(mType, mLoc, mDir);   // unchanged, incl. the inf check
if (!mpScat->IsAnisotropic()) {
  /* existing scatlen / edgelen block, verbatim */
} else {
  Real edgelen = travel.PathLength, s = 0;  bool scattered = false;
  while (true) {
    s += mpScat->GetRandomPathLength(mType);
    if (!(s < edgelen)) break;                                     // boundary first
    TravelRec tent = mpCell->AdvanceLength(mType, s, mLoc, mDir);
    Scatterer::AnisoAttempt att =
        mpScat->AttemptAnisotropicScatter(mType, tent.NewLoc, tent.NewDir, mPol);
    if (att.IsReal) { Move(tent); mPol = att.InPol; Transform(att.Rel);
                      dataout.ReportScatterEvent(*this); scattered = true; break; }
  }                                                                // null collision: no move, no report
  if (scattered) continue;
}
this->Move(travel);   // existing boundary handling continues unchanged
```

- [ ] **Step 4: Run `tests/test_aniso_scenarios.sh`.** Expect
  `anisotropic scenario checks passed`.

- [ ] **Step 5: Run the full check, then commit.** The isotropic regression proves the
  isotropic branch is untouched.

```bash
git add phonons.cpp tests/test_aniso_scenarios.sh
git commit -m "Propagate anisotropic scattering with null collisions"
```

---

### Task 6: Anisotropic scatterer reporting

**Files:**
- Modify: `scatterers.hpp`, `scatterers.cpp` (`PrintStats`, `PrintAllScatteringStats`,
  constructor), `doc/MANUAL.md`, `tests/test_aniso_inputs.sh`
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Consumes: `IntegratedRate` (Task 4).
- Produces:
  - `void Scatterer::IntegrateRateAndDipole(raytype in, const R3::XYZ & dir, const R3::XYZ & pol, const S2::S2Set & dirs, Real & rate, Real & dipole) const`.
    `IntegratedRate` becomes a wrapper around it. The dipole is `Σ cosψ·G / Σ G`, as in
    `ComputeDipoles`.
  - `const std::vector<Real> & Scatterer::GetReportMFPs() const` returns the directional mean
    free paths in dump order.
  - The report direction set is a static `S2::TesselSphere(S2::TESS_ICO, 6)`, built on the
    first anisotropic construction during the single-threaded model build.

- [ ] **Step 1: Write the failing tests.**

```cpp
test_report_directional_mfps: LopNorCrust(2,2,0.5) at degree 5; GetReportMFPs() (AXI order:
    P0 P45 P90 S0 S45SV S45SH S90SV S90SH): P0=1603, P90=1203, S0=755, S90SV=443, S90SH=375 each within 3%
test_report_ani_order: LopNorCrust(2,1,0.5): 9 values ordered P_E P_N P_U S_E/N S_E/U S_N/E S_N/U S_U/E S_U/N;
    P_E < P_N (propagation along the longer axis scatters more, design table: 1426 vs 1724 for kappa 0.5)
```

  Extend `tests/test_aniso_inputs.sh` with a zero-phonon run of the anisotropic Lop Nor
  layout. In the scatterer dump:
  - isotropic rows have 11 fields;
  - `AXI` rows have 22 fields, with field 14 `AXI`;
  - `ANI` rows have 23 fields, with field 14 `ANI`;
  - field 3 equals `ax` (`2` for crust rows);
  - `awk '{print $5*$9}'` is positive on every row;
  - the default layout's dump header has no `ay`.

- [ ] **Step 2: Run both test scripts.** Expect failures: `GetReportMFPs` is undefined, and the
  field counts are wrong.

- [ ] **Step 3: Implement.** At construction, an anisotropic scatterer evaluates 9 axis states
  (P along E, N, U; S along each axis with pol along each other axis) and, for `AXI`, the 8
  angle states above, over the degree-6 set.
  - Legacy `MFP (P)` is `3/Σ` of the 3 P rates, and `MFP (S)` is `6/Σ` of the 6 S rates.
    `DM` is the mean of the corresponding dipoles.
  - `PrintStats` prints `GetAx()` in the `a` column; for isotropic scatterers that is
    identical to today. It then appends `ay az MODE` and the directional values (8 for `AXI`,
    9 for `ANI`) at the existing widths.
  - `PrintAllScatteringStats` appends ` ay az mode directional-MFPs` to its column-header line
    only when some scatterer `IsAnisotropic()`.
  - Document the appended columns, in order, in `doc/MANUAL.md`.

- [ ] **Step 4: Run `tests/test_anisotropy_unit.sh` and `tests/test_aniso_inputs.sh`.** Both
  must pass.

- [ ] **Step 5: Run the full check, then commit.**

```bash
git add scatterers.hpp scatterers.cpp doc/MANUAL.md tests/test_aniso_inputs.sh tests/anisotropy_tests.cpp
git commit -m "Report anisotropic scale lengths and directional mean free paths"
```

---

### Task 7: Axisymmetric sampler

**Files:**
- Create: `aniso_axi.hpp`, `aniso_axi.cpp`
- Modify: `probability.hpp`, `probability.cpp`, `scatterers.hpp`, `scatterers.cpp`, `Makefile`
  (add `aniso_axi.o` to `objects`, the
  `aniso_axi_hpp = aniso_axi.hpp $(scatparams_hpp) $(probability_hpp) $(geom_hpp)` variable,
  its rule, and add it to the `scatterers.o` and test-target prerequisites)
- Test: `tests/anisotropy_tests.cpp`

**Interfaces:**
- Consumes: `ScatterParams` (Task 2), `Scatterer::RotationForS1` and `ToFabric` (Task 4).
- Produces, in `ProbDist`, all `const` and requiring integrated form:
  - `bool IsIntegrated() const`
  - `Real GetCumulative(Index idx) const`
  - `Real GetRangeMagnitude(Index begin, Index end) const`
  - `Index GetRandomIndexInRange(Index begin, Index end) const`, which needs a range
    magnitude greater than 0.
- Produces, in `aniso_axi.hpp`:

```cpp
class AxiProposal {   // read-only after construction; conv index order is Scatterer's GPP, GPS, GSP, GSS
public:
  static const int kCoarseDegree = 3, kRepDegree = 6, kThetaBins = 45, kChiBins = 9;  // 2 deg, 10 deg, folded to [0,90]
  static constexpr Real kBoundMargin = 1e-6, kWeightFloor = 1e-9;
  AxiProposal(const ScatterParams & par, const S2::S2Set & toa, const std::vector<ProbDist> & pdists);
  Real MajorantRate(raytype intype) const;                       // g*_max (1/length)
  Real BinRate(raytype intype, const R3::XYZ & u_rel) const;     // g*_bin of u_rel's folded bin
  struct Draw { int Conv; Index Toa; Real Weight; Real CellSum; };  // acceptance = AnisoRatio*CellSum/Weight
  Draw DrawProposal(raytype intype, const R3::XYZ & u_rel) const;
  Real AcceptanceBound(raytype intype, const R3::XYZ & u_rel, int conv, Index cell) const; // Weight/CellSum
  Index CellOf(Index toa) const;  Index NumCells() const;
};
```

- Produces, in `Scatterer`:
  - `static Real CanonicalPPolarization(const S2::ThetaPhi & dir, const R3::XYZ & up)` returns
    `RotationForS1(dir, (up·d)d − up)`.
  - `SAMPLER_AUTO` selects `SAMPLER_AXISYMMETRIC` for `IsAxisymmetric()` specs.
  - Under the axisymmetric sampler, `mMeanFreeP[in] = 1/MajorantRate(in)`, which is `+inf`
    when the rate is 0.

- [ ] **Step 1: Write the failing tests.**

```cpp
test_probdist_range: weights 1..10: GetCumulative(3)==10, GetRangeMagnitude(2,5)==12; 100000 GetRandomIndexInRange(2,5)
    frequencies within 4 sigma of 3:4:5/12 and always in [2,5)
test_axi_bounds_valid (4g): specs (2,2,0.5), (5,5,0.5), (0.5,0.5,2) x kappa {0.3,0.5} x f {1,4} Hz at degree 5:
    for u_rel at every P bin's corners+center (theta to 180 deg) and S bins' corners+centers in all four chi
    quadrants (a 3x3 subset of bins incl. edges): R(k)*CellSum(actual cell) <= Weight*(1+1e-12) at every fine k
test_axi_matches_general (4f): (2,2,0.5): AXISYMMETRIC vs GENERAL on the 4a in-states: mean distances within
    4*sqrt(sigma_a^2+sigma_g^2); two-sample chi2/dof < 1.5 on (Conv, degree-2 cell) bins
test_axi_draws_per_event (4h): degree 6; mean AnisoAttempt.Draws per real event for P at theta 0,30,60,89 and
    S at 30,89 with chi 0,45,90; write all values (both samplers) to build/aniso_draws.txt;
    AXI <= 2.2 (P) / 2.6 (S) for (2,2,0.5); <= 3.6 / 6.0 for (5,5,0.5); AXI < GENERAL in every state
test_axi_mfp_curves (4e): (2,2,0.5): MFP vs theta every 5 deg (P, SV, SH) via IntegratedRate on the degree-6 set,
    written to build/aniso_mfp_curves.txt; all finite > 0; P(90) < P(0)
test_axi_vertical_incidence: theta exactly 0, P and S: 1000 attempts finite, no throw
test_axi_curved_frame: RAE_CURVED, loc 950 km along x: OrthoAxes(dir, CanonicalPPolarization(dir, up)).S1()
    parallel to (up·d)d-up within 1e-12 for dir Dir(70,20) rotated into the local frame; restore ENU_ORTHO
test_axi_low_degree: degree 2: NumCells()==320 (one fine TOA per cell); bounds valid as in 4g; 1000 real events
test_axi_eps_zero: eps=0: MajorantRate==0, GetMajorantMeanFreePath is +inf, construction does not throw
test_axi_extreme_ratio: (10,10,0.5) at degree 4: bounds valid; mean Draws per real event finite (< 50)
isolated_axi_overridemfp: Scatterer::OverrideMFP(100,50) then axisymmetric LopNorCrust(2,2,0.5): mean
    GetRandomPathLength(RAY_P) over 20000 = 100 within 3%; every attempt IsReal; (Conv, degree-2 cell) distribution
    matches a SAMPLER_GENERAL scatterer under the same override: two-sample chi2/dof < 1.5
```

- [ ] **Step 2: Run `tests/test_anisotropy_unit.sh`.** Expect a compile failure on
  `AxiProposal`, `GetRangeMagnitude`, `SAMPLER_AXISYMMETRIC`, and the other new names.

- [ ] **Step 3: Add the `ProbDist` range methods.** Do binary search over `mDist` within
  `[begin, end)`, drawing with `RandomUniform()`. `GetRandomIndex` stays untouched.

- [ ] **Step 4: Implement `AxiProposal`.**
  - **Cells.** With `D = log4(toa.size()/20)`, a cell owns fine indices
    `[j·4^(D−3), (j+1)·4^(D−3))`; below degree 3 each fine index is a cell. Build the x- and
    y-mirror cell permutations from a degree-`min(D,3)` `TesselSphere`, by matching mirrored
    centers within 1e-9, and throw if any cell has no match.
  - **Representatives.** Use `TesselSphere(degree min(D,6))`. Representative `i` owns fine block
    `[i·4^(D−6), …)`. Its radius `r_i` is the maximum angle to its fine points, or 0 when
    `D ≤ 6`.
  - **Bins.** θ bin `b` covers `[2b°, 2b+2°]`; χ bin `c` covers `[10c°, 10c+10°]`. The cap radius
    is `δ = Δθ/2` for P and `δ = Δθ/2 + sin(θ_hi)·Δχ/2` for S.
  - **Bounds** for each representative `t` and conversion `(k_in, k_out)`, with
    `u_c(θ,χ) = (−sinθ cosχ, sinθ sinχ, cosθ)`:

```text
m_rep = k_out*t − k_in*ẑ;  mr = |m_rep|;  rho = k_out*2*sin(r_i/2)
α_c = angle(m_rep, u_c) (0 if mr == 0);  α ∈ [max(0, α_c−δ), min(π, α_c+δ)]
s_lo = mr*cos(α_hi) − rho;  s_hi = mr*cos(α_lo) + rho
s2_hi = max(s_lo², s_hi²);  s2_lo = (s_lo ≤ 0 ≤ s_hi) ? 0 : min(s_lo², s_hi²)
m2_lo = max(mr−rho, 0)²;  m2_hi = (mr+rho)²
Q_lo = ah ≥ az ? max(ah²·m2_lo − (ah²−az²)·min(s2_hi, m2_hi), az²·m2_lo) : ah²·m2_lo + (az²−ah²)·s2_lo
bound = min(M, M·((1 + a_min²·m2_hi)/(1 + Q_lo))^(κ+3/2))
```

  - **Weights.** A cell's bound is the maximum over its representatives. The weight is
    `w = max(bound·(1+kBoundMargin)·S, kWeightFloor·Σ_bin bound·S)` when `S > 0`, and 0
    otherwise, where `S = pdists[conv].GetRangeMagnitude(cell range)`.
  - **Storage.** Store only the per-bin cumulative `double` array over `(conv, cell)` for the
    in-type's two conversions; recover `w` as a difference. Then
    `g*_bin = Σw/toa.size()`, and `g*_max` is the maximum over bins.
  - **Fold**, given `u = u_rel`:

```text
if (u.z < 0) u = −u                            // R depends on (m·u)²
mx = (u.x > 0); if (mx) u.x = −u.x
P: ignore u.y (canonical frame, χ bin 0).   S: my = (u.y < 0); if (my) u.y = −u.y
θ = acos(clamp(u.z)); χ = sinθ > 1e-12 ? atan2(u.y, −u.x) : 0;  bins clamped to [0,44], [0,8]
actual cell = Mx^mx(My^my(table cell))         // mirrors are involutions and commute
```

  - **Draw.** Pick a table entry `(conv, table cell)` from the bin's cumulative array with
    `RandomUniform()`. Map it to the actual cell. Draw
    `Toa = pdists[conv].GetRandomIndexInRange(actual cell range)`. Return
    `{conv, Toa, w, S(actual cell)}`.

- [ ] **Step 5: Wire the axisymmetric path into `AttemptAnisotropicScatter`.**
  1. Compute `up = ECS.GetUp(loc)`.
  2. Use `pol_used = (intype==RAY_P) ? CanonicalPPolarization(dir, up) : pol`.
  3. Build `AA` from `pol_used`, and set `u_rel = (up·S1, up·S2, up·S3)`.
  4. With `Draws = 0`: if `RandomUniform()·MajorantRate ≥ BinRate(u_rel)`, return a null
     attempt.
  5. Otherwise call `DrawProposal`, set `Draws = 1`, and compute `m` as in Task 4.
  6. Set `p = AnisoRatio(ToFabric(loc,m))·CellSum/Weight`. If `p > 1 + kBoundMargin`, throw
     `std::logic_error` with the spec, `u_rel`, `conv`, and `Toa`.
  7. Accept if `RandomUniform() < p`, and set `InPol = pol_used`.

  With the override, skip thinning and loop steps 5–7 until accepted. Build the `AxiProposal`
  in the constructor when the sampler is `SAMPLER_AXISYMMETRIC`.

- [ ] **Step 6: Run `tests/test_anisotropy_unit.sh` and `tests/test_aniso_scenarios.sh`.** Both
  must pass. Axisymmetric regions in the scenarios now use this sampler.

- [ ] **Step 7: Run the full check, then commit.**

```bash
git add aniso_axi.hpp aniso_axi.cpp probability.hpp probability.cpp scatterers.hpp scatterers.cpp Makefile tests/anisotropy_tests.cpp
git commit -m "Add incidence-tailored axisymmetric sampler"
```

---

### Task 8: Performance measurement and final verification

**Files:**
- Create: `tests/benchmark_aniso.sh`, `docs/anisotropy-performance.md`

**Interfaces:**
- Consumes: `run_aniso_case.sh` (Task 1), and `build/aniso_draws.txt` written by
  `test_axi_draws_per_event` (Task 7).

- [ ] **Step 1: Write `tests/benchmark_aniso.sh`.** It is manual and not named `test_*`, so the
  full check skips it.
  - Environment: `BENCH_DEGREE` (9), `BENCH_LOPNOR_PHONONS` (1000000),
    `BENCH_CRUSTPINCH_PHONONS` (200000), `BENCH_THREADS` (`auto`).
  - For each scenario and ε in {do-script value, 0.05}, time the isotropic run and the runs
    for aspect ratio r ∈ {2, 4, 10}:
    - axisymmetric: `ax=ay=r·a`, `az=a`;
    - triaxial: `ax=r·a`, `ay=√r·a`, `az=a`;
    - use the do-script `a` per region, and apply the anisotropy to the crust and mantle
      regions.
  - Print a Markdown table with these columns: scenario, ε, mode, r, seconds, slowdown versus
    isotropic, and gate (`ok` if the slowdown is at most 2×, otherwise `GATE`).
  - Isotropic no-slowdown check: `BASELINE_EXE` is required. Time each scenario's default
    isotropic case three times with `BASELINE_EXE` and three times with the current build.
    Print the median current/baseline ratio, and exit nonzero if it exceeds 1.05.

- [ ] **Step 2: Run the measurements and write them up.** Run
  `tests/test_anisotropy_unit.sh axi_draws` and `tests/benchmark_aniso.sh`. Record in
  `docs/anisotropy-performance.md`:
  - date, host, compiler, and workload, in the README's table style;
  - the benchmark table;
  - the isotropic current/baseline ratio;
  - draws per real event from `build/aniso_draws.txt` for both samplers;
  - the gate verdict.

  If any row reads `GATE`, say in the doc which follow-up the spec's measurement gate calls
  for. Do not implement it.

- [ ] **Step 3: Run the full check one last time.** Every `tests/test_*.sh` passes, including
  `test_lopnor_thread_equivalence.sh` unchanged.

- [ ] **Step 4: Commit.**

```bash
git add tests/benchmark_aniso.sh docs/anisotropy-performance.md
git commit -m "Measure anisotropic scattering performance against the gate"
```
