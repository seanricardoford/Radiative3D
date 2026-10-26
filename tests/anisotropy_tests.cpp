// tests/anisotropy_tests.cpp
//
// Minimal, dependency-free unit-check harness for the anisotropic
// heterogeneity scale-length work (ax, ay, az in Elastic::HetSpec and
// ScatterParams). Built by `make anisotropy-tests` into
// build/anisotropy_tests and driven by tests/test_anisotropy_unit.sh.
// This file is expected to grow as later tasks of the
// anisotropic-scale-lengths plan add their own checks; see
// docs/superpowers/plans/2026-10-16-anisotropic-scale-lengths.md and
// docs/superpowers/specs/2026-10-16-anisotropic-scale-lengths-design.md.
//
// USAGE:
//
//   build/anisotropy_tests            Runs every test whose name does
//                                      not start with "isolated_".
//   build/anisotropy_tests --list     Prints every registered test
//                                      name, one per line.
//   build/anisotropy_tests FILTER     Runs tests whose name contains
//                                      FILTER. A test named
//                                      "isolated_*" only runs when
//                                      FILTER names it exactly -- such
//                                      tests need a fresh process
//                                      because the overrides they
//                                      exercise are class-static and
//                                      can't be reset within one run.
//                                      (See tests/test_anisotropy_unit.sh,
//                                      which drives each isolated_*
//                                      test in its own process.)
//
// Each run test prints "PASS name" or "FAIL name: message"; the program
// exits 1 if any run test failed, 0 otherwise.
//
// This is deliberately not a general-purpose framework: no fixtures, no
// suites, no external dependencies. Just a registry of named functions
// that signal failure by throwing TestFailure (directly, or via the
// CHECK()/CHECK_MSG() macros below).
//
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "elastic.hpp"
#include "geom.hpp"
#include "geom_s2.hpp"
#include "parallel.hpp"
#include "scatparams.hpp"
#include "sources.hpp"
#include "phonons.hpp"
#include "scatterers.hpp"
#include "aniso_axi.hpp"
#include "ecs.hpp"


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  MINIMAL TEST FRAMEWORK                                      ****
// ****                                                              ****
//

// ::: TestFailure :::
//
//   Thrown by CHECK()/CHECK_MSG() (or directly by a test body) to
//   report a failed check. Any other uncaught exception escaping a
//   test, or none at all, is also handled by the runner in main().
//
struct TestFailure {
  std::string message;
  explicit TestFailure(const std::string & msg) : message(msg) {}
};

typedef void (*TestFunc)();

struct TestCase {
  std::string name;
  TestFunc func;
  TestCase(const std::string & n, TestFunc f) : name(n), func(f) {}
};

// Registry of all tests, populated by the static Registrar objects that
// each TEST(...) block below creates. A function-local static avoids
// any dependency on static-initialization order (moot here, since
// everything lives in this one translation unit, but it's the robust
// idiom regardless).
//
static std::vector<TestCase> & Registry() {
  static std::vector<TestCase> reg;
  return reg;
}

struct Registrar {
  Registrar(const std::string & name, TestFunc func) {
    Registry().push_back(TestCase(name, func));
  }
};

// TEST(name) { ... } declares and registers a test function named
// test_<name> under the string "<name>", then opens the function body
// for the caller to supply.
//
#define TEST(NAME)                                                \
  static void test_##NAME();                                      \
  static Registrar registrar_##NAME(#NAME, test_##NAME);          \
  static void test_##NAME()

#define CHECK(COND)                                                \
  do {                                                              \
    if (!(COND)) {                                                  \
      std::ostringstream _check_oss;                                \
      _check_oss << "check failed: " #COND                          \
                 << " (" << __FILE__ << ":" << __LINE__ << ")";     \
      throw TestFailure(_check_oss.str());                          \
    }                                                                \
  } while (0)

#define CHECK_MSG(COND, MSG)                                       \
  do {                                                               \
    if (!(COND)) {                                                   \
      std::ostringstream _check_oss;                                 \
      _check_oss << MSG                                              \
                 << " (" << __FILE__ << ":" << __LINE__ << ")";      \
      throw TestFailure(_check_oss.str());                           \
    }                                                                 \
  } while (0)

// NearRel(actual, expected, reltol):
//
//   True if actual is within relative tolerance reltol of expected.
//   Falls back to an absolute comparison against reltol itself when
//   expected is exactly zero.
//
static bool NearRel(Real actual, Real expected, Real reltol) {
  if (actual == expected) return true;
  Real scale = std::fabs(expected);
  if (scale == 0.0) return std::fabs(actual) <= reltol;
  return std::fabs(actual - expected) <= reltol * scale;
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  SHARED FIXTURE HELPERS                                      ****
// ****                                                              ****
//

// ::: UseTOADegree(degree) :::
//
//   Ensures a take-off-angle set (S2::TesselSphere) of the given
//   icosahedral tessellation degree exists, building and caching it on
//   first request (one TesselSphere per distinct degree, kept for the
//   life of the program -- rebuilding a high-degree sphere repeatedly
//   would be wasteful across tests that share a degree), points
//   PhononSource at it via Set_TOA_Array(), and returns it by reference
//   so callers can also iterate its directions directly.
//
static const S2::TesselSphere & UseTOADegree(int degree) {
  static std::map<int, S2::TesselSphere *> cache;

  S2::TesselSphere *& slot = cache[degree];
  if (slot == nullptr) {
    slot = new S2::TesselSphere(S2::TESS_ICO, degree);
  }
  PhononSource::Set_TOA_Array(slot);
  return *slot;
}

// ::: LopNorCrust(ax,ay,az,kappa,eps,vp,vs) :::
//
//   A ScatterParams reference fixture at Lop Nor crustal velocities
//   (Vp 6.2, Vs 3.6 km/s by default) and 2 Hz, with caller-chosen
//   heterogeneity scale lengths. nu is fixed at 0.8, matching the Lop
//   Nor crust layer (see user_LopNorCylMoho_inc.cpp). The trailing
//   vp/vs parameters default to the crust values so every existing call
//   site is unaffected; test_p_frame_consistency passes the Lop Nor
//   mantle velocities through them.
//
static ScatterParams LopNorCrust(Real ax, Real ay, Real az,
                                  Real kappa=0.5, Real eps=0.01,
                                  Real vp=6.2, Real vs=3.6) {
  ScatterParams::SetFrequencyHertz(2.0);
  return ScatterParams(Elastic::VpVs(vp,vs),
                        Elastic::HSnexyzk(0.8,eps,ax,ay,az,kappa));
}

// ::: Dir(theta_deg, phi_deg) :::
//
//   A unit vector at the given take-off-angle-style spherical
//   coordinates (theta measured from +z, phi from +x towards +y), in
//   degrees.
//
static R3::XYZ Dir(Real theta_deg, Real phi_deg) {
  return S2::ThetaPhi(theta_deg * Geometry::DtoR,
                       phi_deg   * Geometry::DtoR).XYZ();
}

// ::: RandomXYZ(rng, halfwidth) :::
//
//   An R3::XYZ with each component drawn independently and uniformly
//   from [-halfwidth,+halfwidth]. Used to spot-check the anisotropic
//   spectrum identities at arbitrary wavevectors. Callers pass a
//   locally-seeded generator so tests stay deterministic.
//
static R3::XYZ RandomXYZ(std::mt19937 & rng, Real halfwidth) {
  std::uniform_real_distribution<Real> d(-halfwidth, halfwidth);
  return R3::XYZ(d(rng), d(rng), d(rng));
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 2 TESTS: HetSpec and ScatterParams scale lengths       ****
// ****                                                              ****
//

TEST(hetspec_modes) {

  Elastic::HSneak iso(0.8, 0.01, 0.5, 0.3);
  CHECK(iso.IsIsotropic());
  CHECK(!iso.IsAxisymmetric());
  CHECK(iso.ax() == 0.5);
  CHECK(iso.ay() == 0.5);
  CHECK(iso.az() == 0.5);
  CHECK(iso.a()  == 0.5);

  Elastic::HSnexyzk axisym(0.8, 0.01, 2, 2, 0.5, 0.3);
  CHECK(axisym.IsAxisymmetric());
  CHECK(!axisym.IsIsotropic());

  Elastic::HSnexyzk triaxial(0.8, 0.01, 2, 1, 0.5, 0.3);
  CHECK(!triaxial.IsIsotropic());
  CHECK(!triaxial.IsAxisymmetric());

  Elastic::HSnexyzk isoxyzk(0.8, 0.01, 0.5, 0.5, 0.5, 0.3);
  CHECK(isoxyzk.IsIsotropic());
  CHECK(isoxyzk.a() == 0.5);

}


TEST(hetspec_validation) {

  bool threw;

  threw = false;
  try { Elastic::HSnexyzk bad(0.8, 0.01, 0.0, 1.0, 1.0, 0.5); }
  catch (const std::invalid_argument &) { threw = true; }
  CHECK_MSG(threw, "ax=0 should throw std::invalid_argument");

  threw = false;
  try { Elastic::HSnexyzk bad(0.8, 0.01, 1.0, -1.0, 1.0, 0.5); }
  catch (const std::invalid_argument &) { threw = true; }
  CHECK_MSG(threw, "ay=-1 should throw std::invalid_argument");

  threw = false;
  try { Elastic::HSnexyzk bad(0.8, 0.01, 1.0, 1.0, 0.0, 0.5); }
  catch (const std::invalid_argument &) { threw = true; }
  CHECK_MSG(threw, "az=0 should throw std::invalid_argument");

}


TEST(scatparams_iso_unchanged) {

  ScatterParams::SetFrequencyHertz(2.0);
  Elastic::Velocity vpvs = Elastic::VpVs(6.2,3.6);
  ScatterParams spFour(vpvs, Elastic::HSneak(0.8,0.01,0.5,0.3));
  ScatterParams spSix (vpvs, Elastic::HSnexyzk(0.8,0.01,0.5,0.5,0.5,0.3));

  CHECK(spFour.CompareRoughly(spSix) == 0.0);
  CHECK(spFour.GetA() == 0.5);
  CHECK(spSix.GetA()  == 0.5);

  // Bitwise GSATO agreement at a spread of take-off angles: with ax ==
  // ay == az == a for both specs, GSATO/XSATO/PSATO never touch ax, ay,
  // or az, so routing isotropic specs through the new triaxial fields
  // must not perturb a single bit of the result.
  const S2::TesselSphere & toas = UseTOADegree(2);
  const std::size_t n = toas.size();
  const std::size_t stride = n / 20;
  for (std::size_t i = 0; i < 20; ++i) {
    S2::S2Point toa = toas[i*stride];
    Real gppA,gpsA,gspA,gssA,spolA;
    Real gppB,gpsB,gspB,gssB,spolB;
    spFour.GSATO(toa, gppA,gpsA,gspA,gssA,spolA);
    spSix.GSATO (toa, gppB,gpsB,gspB,gssB,spolB);
    CHECK(gppA == gppB);
    CHECK(gpsA == gpsB);
    CHECK(gspA == gspB);
    CHECK(gssA == gssB);
    CHECK(spolA == spolB);
  }

}


TEST(scatparams_aniso_reference) {

  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  CHECK(sp.GetA()  == 0.5);
  CHECK(sp.GetAx() == 2.0);
  CHECK(sp.GetAz() == 0.5);
  CHECK(sp.AnisoMajorantFactor() == 16.0);

}


TEST(psato_aniso_iso_limit) {

  ScatterParams sp = LopNorCrust(0.5, 0.5, 0.5);  // isotropic spec
  std::mt19937 rng(20261016u);

  for (int i = 0; i < 50; ++i) {
    R3::XYZ m = RandomXYZ(rng, 8.0);
    Real paniso = sp.PSATO_Aniso(m);
    Real piso   = sp.PSATO(m.Mag());
    CHECK(NearRel(paniso, piso, 1e-12));
  }

}


TEST(aniso_ratio_bounds) {

  ScatterParams sp = LopNorCrust(2, 2, 0.5, 0.3);
  const Real M = 16.0;  // ax*ay*az/a^3 = 2*2*0.5 / 0.5^3

  CHECK(NearRel(sp.AnisoRatio(R3::XYZ(0,0,0)), M, 1e-12));

  // Along the shortest fabric axis (here z, since az==a==0.5), the
  // denominator terms of PSATO_Aniso and PSATO happen to match exactly
  // for every magnitude, so the ratio stays pinned at M.
  Real zvals[] = {0.5, 1.0, 4.0};
  for (int i = 0; i < 3; ++i) {
    R3::XYZ m = Dir(0.0, 0.0).ScaledBy(zvals[i]);  // (0,0,zvals[i])
    CHECK(NearRel(sp.AnisoRatio(m), M, 1e-12));
  }

  std::mt19937 rng(777u);
  for (int i = 0; i < 1000; ++i) {
    R3::XYZ m = RandomXYZ(rng, 8.0);
    Real ratio = sp.AnisoRatio(m);
    CHECK(ratio <= M * (1.0 + 1e-12));
    Real acceptance = sp.AnisoAcceptance(m);
    CHECK(NearRel(acceptance, ratio/M, 1e-12));
  }

}


TEST(matching_distinguishes) {

  ScatterParams spA    = LopNorCrust(2, 2, 0.5);
  ScatterParams spPerm = LopNorCrust(2, 0.5, 2);      // axes permuted
  ScatterParams spIso  = LopNorCrust(0.5, 0.5, 0.5);
  ScatterParams spA2   = LopNorCrust(2, 2, 0.5);      // identical to spA

  CHECK(spA.CompareRoughly(spPerm) > 0.0);
  CHECK(spA.CompareRoughly(spIso)  > 0.0);
  CHECK(spA.CompareRoughly(spA2) == 0.0);

}


TEST(scatterer_guard) {

  // When a Scatterer is constructed with anisotropic HetSpec parameters,
  // GetRandomScatteredRelativePhonon() must throw std::logic_error rather
  // than silently proceeding with isotropic sampling.

  ScatterParams aniso = LopNorCrust(2, 2, 0.5);  // anisotropic: ax != az
  CHECK(!aniso.IsIsotropic());

  UseTOADegree(2);  // Set up TOA array for Scatterer
  Scatterer scatt(aniso);

  bool threw = false;
  try {
    scatt.GetRandomScatteredRelativePhonon(RAY_P);
  } catch (const std::logic_error & e) {
    const std::string msg(e.what());
    if (msg.find("anisotropic scatterers must be sampled through AttemptAnisotropicScatter") != std::string::npos) {
      threw = true;
    }
  }
  CHECK_MSG(threw, "GetRandomScatteredRelativePhonon should throw std::logic_error for anisotropic params");

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 4 TESTS: General anisotropic sampler                   ****
// ****                                                              ****
//

// ::: MeanDistanceToRealEvent(s, in, dir, pol, n) :::
//
//   Sums GetRandomPathLength(in) draws until AttemptAnisotropicScatter()
//   returns a real (accepted) event, and averages the per-event total
//   distance over n events. Thinning preserves memorylessness, so this
//   mean should converge to 1/IntegratedRate(in, dir.XYZ(), pol, ...).
//
static Real MeanDistanceToRealEvent(Scatterer & s, raytype in,
                                    const S2::ThetaPhi & dir,
                                    const R3::XYZ & pol, int n) {
  const Real rot = Scatterer::RotationForS1(dir, pol);
  const R3::XYZ origin(0, 0, 0);
  Real total = 0;
  for (int i = 0; i < n; ++i) {
    Real dist = 0;
    Scatterer::AnisoAttempt att;
    do {
      dist += s.GetRandomPathLength(in);
      att = s.AttemptAnisotropicScatter(in, origin, dir, rot);
    } while (!att.IsReal);
    total += dist;
  }
  return total / n;
}

// ::: MeanDrawsPerRealEvent(s, in, dir, pol, n) :::
//
//   Sums AnisoAttempt.Draws over every AttemptAnisotropicScatter() call
//   (thinned nulls contribute 0, a drawn-but-rejected or accepted
//   proposal contributes its own Draws) made while chasing down n real
//   (accepted) events, and returns the mean total per real event. Used
//   by test_axi_draws_per_event and test_axi_extreme_ratio.
//
static Real MeanDrawsPerRealEvent(Scatterer & s, raytype in,
                                  const S2::ThetaPhi & dir,
                                  const R3::XYZ & pol, int n) {
  const Real rot = Scatterer::RotationForS1(dir, pol);
  const R3::XYZ origin(0, 0, 0);
  long long totalDraws = 0;
  for (int i = 0; i < n; ++i) {
    Scatterer::AnisoAttempt att;
    do {
      att = s.AttemptAnisotropicScatter(in, origin, dir, rot);
      totalDraws += att.Draws;
    } while (!att.IsReal);
  }
  return (Real) totalDraws / n;
}

// ::: CheckDistributionMatchesTarget(...) :::
//
//   Shared body for test 4b ("test_general_distribution_matches_target")
//   and the matching repeat under --overridemfp ("isolated_overridemfp").
//   Draws 'ndraws' *accepted* events from 'scat' for incoming state
//   (intype, dir, pol), bins them by (conversion, coarse take-off-angle),
//   and compares against the expected distribution computed directly
//   from GSATO() (isotropic, a=a_min) times AnisoRatio() at the
//   corresponding exchange wavevector -- the same quantity
//   IntegratedRate() sums, just binned rather than totalled. 'divisor'
//   coarsens the fine (degree-5) take-off-angle index down to a coarser
//   grid (e.g. 4^(5-2)=64 maps degree 5 down to degree 2) so that bins
//   collect enough counts for a meaningful chi-square test.
//
//   Conversions are identified by the raw out_types_e integer codes
//   Scatterer::AnisoAttempt::Conv uses (GPP=0, GPS=1, GSP=2, GSS=3),
//   since out_types_e itself is a protected implementation detail of
//   Scatterer.
//
static void CheckDistributionMatchesTarget(Scatterer & scat,
                                           const ScatterParams & sp,
                                           raytype intype,
                                           const S2::ThetaPhi & dir,
                                           const R3::XYZ & pol,
                                           const S2::TesselSphere & toas,
                                           int divisor, int ndraws) {
  const Real rot = Scatterer::RotationForS1(dir, pol);
  const R3::OrthoAxes AA(dir.Theta(), dir.Phi(), rot);
  const R3::XYZ dirXYZ = dir.XYZ();
  const Real k_in = (intype == RAY_P) ? sp.GetKP() : sp.GetKS();

  int convs[2];
  if (intype == RAY_P) { convs[0] = 0; convs[1] = 1; }  // GPP, GPS
  else                 { convs[0] = 2; convs[1] = 3; }  // GSP, GSS
  const raytype outrt[4] = {RAY_P, RAY_S, RAY_P, RAY_S};

  // Expected weight (proportional to the probability of landing in
  // that bin, given acceptance) summed over the fine take-off-angle
  // grid:
  std::map<std::pair<int,int>, Real> expectedWeight;
  Real totalWeight = 0;
  const std::size_t n = toas.size();
  for (std::size_t t = 0; t < n; ++t) {
    Real gpp, gps, gsp, gss, spol;
    sp.GSATO(toas[t], gpp, gps, gsp, gss, spol);
    const Real gvals[4] = {gpp, gps, gsp, gss};
    const R3::XYZ e_out = AA.Express(toas[t].XYZ());
    const int coarse = (int)(t / divisor);
    for (int ci = 0; ci < 2; ++ci) {
      const int conv = convs[ci];
      const Real k_out = (outrt[conv] == RAY_P) ? sp.GetKP() : sp.GetKS();
      const R3::XYZ m = e_out.ScaledBy(k_out) + dirXYZ.ScaledBy(-k_in);
      const Real w = gvals[conv] * sp.AnisoRatio(m);
      expectedWeight[std::make_pair(ci, coarse)] += w;
      totalWeight += w;
    }
  }

  // Observed counts, drawn from the actual sampler:
  std::map<std::pair<int,int>, int> observed;
  const R3::XYZ origin(0, 0, 0);
  int collected = 0;
  while (collected < ndraws) {
    Scatterer::AnisoAttempt att =
      scat.AttemptAnisotropicScatter(intype, origin, dir, rot);
    if (!att.IsReal) continue;
    const int ci = (att.Conv == convs[0]) ? 0 : 1;
    const int coarse = (int)(att.Toa / divisor);
    observed[std::make_pair(ci, coarse)]++;
    ++collected;
  }

  // Chi-square over bins with expected count >= 20:
  Real chi2 = 0;
  int dof = 0;
  for (std::map<std::pair<int,int>, Real>::const_iterator it = expectedWeight.begin();
       it != expectedWeight.end(); ++it) {
    const Real expectedCount = ndraws * (it->second / totalWeight);
    if (expectedCount < 20.0) continue;
    std::map<std::pair<int,int>, int>::const_iterator oit = observed.find(it->first);
    const int obs = (oit == observed.end()) ? 0 : oit->second;
    const Real diff = obs - expectedCount;
    chi2 += diff * diff / expectedCount;
    ++dof;
  }
  CHECK_MSG(dof > 0, "no bins had expected count >= 20");
  const Real chi2_per_dof = chi2 / dof;
  CHECK_MSG(chi2_per_dof < 1.5,
           "chi2/dof too large: " << chi2_per_dof << " (dof=" << dof << ")");
}


TEST(general_iso_forced) {

  UseTOADegree(5);
  ScatterParams sp = LopNorCrust(0.5, 0.5, 0.5);  // isotropic spec
  Scatterer general(sp, Scatterer::SAMPLER_GENERAL);
  Scatterer autoSc(sp, Scatterer::SAMPLER_AUTO);

  CHECK(general.IsAnisotropic());
  CHECK(!autoSc.IsAnisotropic());
  CHECK(NearRel(general.GetMajorantMeanFreePath(RAY_P),
               autoSc.GetMajorantMeanFreePath(RAY_P), 1e-12));
  CHECK(NearRel(general.GetMajorantMeanFreePath(RAY_S),
               autoSc.GetMajorantMeanFreePath(RAY_S), 1e-12));

  const R3::XYZ origin(0, 0, 0);
  const S2::ThetaPhi dir(60.0 * Geometry::DtoR, 30.0 * Geometry::DtoR);
  for (int i = 0; i < 2000; ++i) {
    Scatterer::AnisoAttempt att =
      general.AttemptAnisotropicScatter(RAY_P, origin, dir, 0.0);
    CHECK(att.IsReal);
  }

}


TEST(general_rate_matches_integration) {

  const S2::TesselSphere & toas = UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);

  const R3::XYZ E(1,0,0), N(0,1,0), U(0,0,1), D(Dir(60,30));

  struct Case { raytype in; R3::XYZ dir; R3::XYZ pol; };
  const Case cases[] = {
    {RAY_P, E, E}, {RAY_P, N, N}, {RAY_P, U, U}, {RAY_P, D, D},
    {RAY_S, E, N}, {RAY_S, E, U},
    {RAY_S, N, E}, {RAY_S, N, U},
    {RAY_S, U, E}, {RAY_S, U, N},
  };
  const int n = 20000;

  for (std::size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i) {
    const S2::ThetaPhi dirTP(cases[i].dir.Theta(), cases[i].dir.Phi());
    const Real g0 = scat.IntegratedRate(cases[i].in, cases[i].dir,
                                        cases[i].pol, toas);
    const Real expected = 1.0 / g0;
    const Real mean = MeanDistanceToRealEvent(scat, cases[i].in, dirTP,
                                              cases[i].pol, n);
    const Real tol = 4.0 * expected / std::sqrt((Real) n);
    CHECK_MSG(std::fabs(mean - expected) < tol,
             "case " << i << ": mean=" << mean << " expected=" << expected
             << " tol=" << tol);
  }

}


TEST(general_distribution_matches_target) {

  const S2::TesselSphere & toas = UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);
  int divisor = 1;
  for (int i = 0; i < 5-2; ++i) divisor *= 4;  // 4^(5-2) = 64

  const S2::ThetaPhi dirP(60.0 * Geometry::DtoR, 30.0 * Geometry::DtoR);
  CheckDistributionMatchesTarget(scat, sp, RAY_P, dirP, dirP.XYZ(),
                                 toas, divisor, 200000);

  const S2::ThetaPhi dirS(Geometry::Pi90, 0.0);  // E
  const R3::XYZ polU(0, 0, 1);                   // U
  CheckDistributionMatchesTarget(scat, sp, RAY_S, dirS, polU,
                                 toas, divisor, 200000);

}


TEST(general_symmetry) {

  const S2::TesselSphere & toas = UseTOADegree(5);

  {
    ScatterParams sp = LopNorCrust(2, 2, 0.5);  // axisymmetric
    Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);
    const R3::XYZ E(1,0,0), N(0,1,0);
    const Real rE = scat.IntegratedRate(RAY_P, E, E, toas);
    const Real rN = scat.IntegratedRate(RAY_P, N, N, toas);
    CHECK(NearRel(rE, rN, 1e-9));
  }

  {
    ScatterParams sp = LopNorCrust(2, 2, 0.5);  // axisymmetric
    Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);
    const R3::XYZ d1 = Dir(60, 30);
    const R3::XYZ d2 = Dir(60, -30);  // x-mirror
    const Real r1 = scat.IntegratedRate(RAY_P, d1, d1, toas);
    const Real r2 = scat.IntegratedRate(RAY_P, d2, d2, toas);
    CHECK(NearRel(r1, r2, 1e-9));
  }

  {
    ScatterParams sp = LopNorCrust(2, 1, 0.5);  // triaxial
    Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);
    const R3::XYZ E(1,0,0), negE(-1,0,0);
    const Real rE = scat.IntegratedRate(RAY_P, E, E, toas);
    const Real rNegE = scat.IntegratedRate(RAY_P, negE, negE, toas);
    CHECK(NearRel(rE, rNegE, 1e-9));
  }

}


TEST(general_vertical_incidence) {

  const S2::TesselSphere & toas = UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);

  const S2::ThetaPhi dir(0.0, 0.0);  // U
  const R3::XYZ pol(1, 0, 0);        // E
  const Real rot = Scatterer::RotationForS1(dir, pol);
  const R3::XYZ origin(0, 0, 0);

  for (int i = 0; i < 1000; ++i) {
    Scatterer::AnisoAttempt att =
      scat.AttemptAnisotropicScatter(RAY_S, origin, dir, rot);
    CHECK(std::isfinite(att.InPol));
    if (att.IsReal) {
      CHECK(std::isfinite(att.Rel.GetDirection().Theta()));
      CHECK(std::isfinite(att.Rel.GetDirection().Phi()));
      CHECK(std::isfinite(att.Rel.GetPolarization()));
    }
  }

  const Real rate = scat.IntegratedRate(RAY_S, R3::XYZ(0,0,1), pol, toas);
  CHECK(std::isfinite(rate));
  CHECK(rate > 0.0);

}


TEST(general_curved_frame) {

  ECS.SetMapping(EarthCoords::RAE_CURVED);

  const R3::XYZ loc(950.0, 0.0, 0.0);
  const R3::XYZ up = ECS.GetUp(loc);
  const R3::XYZ fab = Scatterer::ToFabric(loc, up);

  CHECK(std::fabs(fab.x()) < 1e-12);
  CHECK(std::fabs(fab.y()) < 1e-12);
  CHECK(std::fabs(fab.z() - 1.0) < 1e-12);
  CHECK(up.z() < std::cos(8.0 * Geometry::DtoR));

  ECS.SetMapping(EarthCoords::ENU_ORTHO);  // restore for later tests

}


TEST(general_eps_zero) {

  UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 2, 0.5, 0.5, 0.0);  // eps=0
  bool threw = false;
  Scatterer * scat = 0;
  try {
    scat = new Scatterer(sp, Scatterer::SAMPLER_GENERAL);
  } catch (...) {
    threw = true;
  }
  CHECK_MSG(!threw, "construction with eps=0 must not throw");
  const Real mfp = scat->GetMajorantMeanFreePath(RAY_P);
  CHECK(std::isinf(mfp));
  CHECK(mfp > 0.0);

  // I3: eps=0 must not leak NaN into the reported (axis-averaged and
  // directional) mean free paths -- AnisoRatio() used to compute
  // PSATO_Aniso(m)/PSATO(|m|), a 0/0 NaN at eps=0. +inf is the correct,
  // isotropic-convention value for a zero scattering rate.
  const std::vector<Real> & mfps = scat->GetReportMFPs();
  CHECK_MSG(mfps.size() == 8, "expected 8 AXI directional MFPs, got "
           << mfps.size());
  for (std::size_t i = 0; i < mfps.size(); ++i) {
    CHECK_MSG(!std::isnan(mfps[i]), "GetReportMFPs()[" << i << "] is NaN");
    CHECK_MSG(mfps[i] > 0.0, "GetReportMFPs()[" << i << "] is not positive");
  }
  delete scat;

}


TEST(general_eps_zero_triaxial_ani) {

  // I3 (triaxial/ANI path): the review's reproduction was specifically a
  // triaxial (ax != ay != az) eps=0 crust row, which exercises the ANI
  // (9-column) branch of ComputeReportStats()/GetReportMFPs(), not just
  // the axisymmetric 8-column branch general_eps_zero already covers.
  UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 1, 0.5, 0.5, 0.0);  // eps=0, triaxial
  CHECK(!sp.IsAxisymmetric());
  CHECK(!sp.IsIsotropic());
  Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);

  const std::vector<Real> & mfps = scat.GetReportMFPs();
  CHECK_MSG(mfps.size() == 9, "expected 9 ANI directional MFPs, got "
           << mfps.size());
  for (std::size_t i = 0; i < mfps.size(); ++i) {
    CHECK_MSG(!std::isnan(mfps[i]), "GetReportMFPs()[" << i << "] is NaN");
    CHECK_MSG(mfps[i] > 0.0, "GetReportMFPs()[" << i << "] is not positive");
  }

}


TEST(isolated_overridemfp) {

  const S2::TesselSphere & toas = UseTOADegree(5);
  Scatterer::OverrideMFP(100, 50);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);

  const int n = 20000;
  Real total = 0;
  for (int i = 0; i < n; ++i) total += scat.GetRandomPathLength(RAY_P);
  const Real mean = total / n;
  CHECK(NearRel(mean, 100.0, 0.03));

  const R3::XYZ origin(0, 0, 0);
  const S2::ThetaPhi dirP(60.0 * Geometry::DtoR, 30.0 * Geometry::DtoR);
  for (int i = 0; i < 2000; ++i) {
    Scatterer::AnisoAttempt att =
      scat.AttemptAnisotropicScatter(RAY_P, origin, dirP, 0.0);
    CHECK(att.IsReal);
  }

  int divisor = 1;
  for (int i = 0; i < 5-2; ++i) divisor *= 4;  // 4^(5-2) = 64
  CheckDistributionMatchesTarget(scat, sp, RAY_P, dirP, dirP.XYZ(),
                                 toas, divisor, 200000);

}


TEST(isolated_nodeflect) {

  UseTOADegree(5);
  Scatterer::SetNoDeflect();
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);

  const R3::XYZ origin(0, 0, 0);
  const S2::ThetaPhi dirP(60.0 * Geometry::DtoR, 30.0 * Geometry::DtoR);
  const Real pol = 0.37;

  // Task 6b: SetNoDeflect() only swaps BuildRelativePhonon()'s result for
  // a non-deflected dummy; it doesn't touch the frame AttemptAnisotropic-
  // Scatter() samples from or the InPol it reports. For P incidence that
  // frame -- and a real event's InPol -- is the canonical one, not the
  // carried 'pol' (see docs/toa-quadrature-bias.md).
  const Real canonicalPolP =
    Scatterer::CanonicalPPolarization(dirP, ECS.GetUp(origin));

  int realCountP = 0;
  for (int i = 0; i < 2000; ++i) {
    Scatterer::AnisoAttempt att =
      scat.AttemptAnisotropicScatter(RAY_P, origin, dirP, pol);
    if (att.IsReal) {
      ++realCountP;
      CHECK(att.Rel.GetDirection().Theta() < 1.e-6);  // nudge_if_singular() keeps
                                                      // theta slightly off zero
      CHECK(att.Rel.GetRaytype() == RAY_P);
      CHECK(att.InPol == canonicalPolP);
    }
  }
  CHECK_MSG(realCountP > 0, "no real P attempts observed");

  int realCountS = 0;
  for (int i = 0; i < 2000; ++i) {
    Scatterer::AnisoAttempt att =
      scat.AttemptAnisotropicScatter(RAY_S, origin, dirP, pol);
    if (att.IsReal) {
      ++realCountS;
      CHECK(att.Rel.GetDirection().Theta() < 1.e-6);  // nudge_if_singular() keeps
                                                      // theta slightly off zero
      CHECK(att.Rel.GetRaytype() == RAY_S);
      CHECK(att.InPol == pol);
    }
  }
  CHECK_MSG(realCountS > 0, "no real S attempts observed");

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 6 TESTS: Anisotropic scatterer reporting               ****
// ****                                                              ****
//

TEST(report_directional_mfps) {

  UseTOADegree(5);  // Degree for the scatterer's own (general-sampler)
                    // tables; GetReportMFPs() always integrates over
                    // its own fixed degree-6 set regardless.
  ScatterParams sp = LopNorCrust(2, 2, 0.5);  // axisymmetric
  CHECK(sp.IsAxisymmetric());
  Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);

  const std::vector<Real> & mfps = scat.GetReportMFPs();
  CHECK_MSG(mfps.size() == 8, "expected 8 AXI directional MFPs, got "
           << mfps.size());

  // Order: P0 P45 P90 S0 S45SV S45SH S90SV S90SH
  //
  // Reference values are the design doc's "Approximate directional mean
  // free paths for 2, 2, 0.5" table, except S90SH: the exact integral
  // there is ~375 (confirmed by an independent 200,000-point
  // Fibonacci-sphere cross-check, which reproduces the design table to
  // <0.1% on every entry including S90SH: 374.7), but this Scatterer's
  // mandated reporting grid -- the equal-weight-node S2::TesselSphere
  // (TESS_ICO, 6) that ComputeMFPs() and the general sampler's mPDists
  // also use -- carries a systematic node-density bias for this
  // particular state that does not shrink with tessellation degree
  // (confirmed at degree 7 and 8: ~397.8 at all three). Per design
  // ruling R9, the reference below (398) is this implementation's own
  // converged value, not the exact integral, so the report states the
  // rate the simulation actually realizes under its own quadrature. See
  // docs/toa-quadrature-bias.md for the bias mechanism and a future-fix
  // plan.
  const Real expected[8] = {1603, 0, 1203, 755, 0, 0, 443, 398};
  const int checkIdx[]   = {0, 2, 3, 6, 7};  // P0, P90, S0, S90SV, S90SH
                                             // -- the reference values
                                             // the design doc gives
                                             // (P45/S45SV/S45SH are not
                                             // independently given
                                             // there).
  for (int i = 0; i < 5; ++i) {
    int idx = checkIdx[i];
    CHECK_MSG(NearRel(mfps[idx], expected[idx], 0.03),
             "directional MFP[" << idx << "]=" << mfps[idx]
             << " expected " << expected[idx] << " within 3%");
  }

}


TEST(report_ani_order) {

  UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 1, 0.5);  // triaxial: ANI mode
  CHECK(!sp.IsAxisymmetric());
  CHECK(!sp.IsIsotropic());
  Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);

  const std::vector<Real> & mfps = scat.GetReportMFPs();
  CHECK_MSG(mfps.size() == 9, "expected 9 ANI directional MFPs, got "
           << mfps.size());

  // Order: P_E P_N P_U S_E/N S_E/U S_N/E S_N/U S_U/E S_U/N
  const Real P_E = mfps[0];
  const Real P_N = mfps[1];
  CHECK_MSG(P_E < P_N,
           "expected P_E < P_N (propagation along the longer fabric "
           "axis scatters more): P_E=" << P_E << " P_N=" << P_N);
  for (std::size_t i = 0; i < mfps.size(); ++i) {
    CHECK(mfps[i] > 0.0);
  }

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 7 TESTS: Axisymmetric sampler                          ****
// ****                                                              ****
//

TEST(probdist_range) {

  ProbDist pd(10);
  for (Index i = 0; i < 10; i++) {
    pd.SetRelativeProb(i, (Real)(i+1));   // weights 1..10
  }
  pd.GetMagnitude();   // Force differential->cumulative conversion
                       // (GetCumulative()/GetRangeMagnitude()/
                       // GetRandomIndexInRange() require it, same as
                       // Scatterer::PopulateWholeProbs() does for
                       // mPDists before AxiProposal is built).

  CHECK(pd.GetCumulative(3) == 10.0);      // 1+2+3+4
  CHECK(pd.GetRangeMagnitude(2,5) == 12.0); // 3+4+5 (indices 2,3,4)

  const int n = 100000;
  int counts[3] = {0,0,0};
  for (int i = 0; i < n; ++i) {
    Index idx = pd.GetRandomIndexInRange(2,5);
    CHECK(idx >= 2 && idx < 5);
    counts[idx-2]++;
  }

  const Real expectedFrac[3] = {3.0/12.0, 4.0/12.0, 5.0/12.0};
  for (int i = 0; i < 3; ++i) {
    const Real p = expectedFrac[i];
    const Real sigma = std::sqrt(p*(1.0-p)/n);
    const Real frac = (Real) counts[i] / n;
    CHECK_MSG(std::fabs(frac - p) < 4.0*sigma,
             "index " << (i+2) << ": frac=" << frac << " expected=" << p
             << " sigma=" << sigma);
  }

}


// ::: BuildIsoPDists(sp, toas, pdists) :::
//
//   Populates 'pdists' (4 ProbDist objects, GPP,GPS,GSP,GSS) the same
//   way Scatterer::PopulateProbDists() does, and forces integration
//   (as Scatterer::PopulateWholeProbs() does via GetMagnitude()) so the
//   result is suitable for passing straight into an AxiProposal
//   constructor, exactly as Scatterer does with its own mPDists.
//
static void BuildIsoPDists(const ScatterParams & sp,
                           const S2::S2Set & toas,
                           std::vector<ProbDist> & pdists) {
  pdists.clear();
  for (int c = 0; c < 4; ++c) {pdists.push_back(ProbDist((int)toas.size()));}
  for (Index k = 0; k < toas.size(); ++k) {
    Real gpp,gps,gsp,gss,spol;
    sp.GSATO(toas[k], gpp,gps,gsp,gss,spol);
    pdists[0].SetRelativeProb(k, gpp);
    pdists[1].SetRelativeProb(k, gps);
    pdists[2].SetRelativeProb(k, gsp);
    pdists[3].SetRelativeProb(k, gss);
  }
  for (int c = 0; c < 4; ++c) {pdists[c].GetMagnitude();} // force Integrate()
}


// ::: RatioForAxisymmetric(sp, m, u_rel) :::
//
//   sp.AnisoRatio() expects 'm' expressed in the fabric (East,North,Up)
//   frame. 'm' here is expressed in the local ray frame instead (as
//   AxiProposal's bound construction works entirely in that frame --
//   see aniso_axi.cpp); 'u_rel' is the local vertical Û, also in the
//   ray frame. Since the dot product m.u_rel (= m's true fabric-Up
//   component) and |m| are both frame-invariant, and an axisymmetric
//   spec's denominator only depends on |m| and the vertical component
//   (any horizontal split works, as ax==ay), reconstructing a fabric-
//   frame stand-in with z = m.u_rel and x = the remaining horizontal
//   magnitude gives the exact same AnisoRatio() value a real
//   ToFabric()-based 'm' would.
//
static Real RatioForAxisymmetric(const ScatterParams & sp, const R3::XYZ & m,
                                 const R3::XYZ & u_rel) {
  const Real s = m.Dot(u_rel);
  const Real h2 = std::max(m.MagSquared() - s*s, (Real)0.0);
  const R3::XYZ mFabric(std::sqrt(h2), 0.0, s);
  return sp.AnisoRatio(mFabric);
}


// ::: CheckAxiBoundsValid(sp, toas, thetaDegsP, thetaChiPairsS, stride) :::
//
//   For each u_rel direction named by 'thetaDegsP' (P incidence angles,
//   degrees, chi implicitly 0) and 'thetaChiPairsS' (S (theta,chi)
//   incidence pairs, degrees), and for every stride'th fine take-off
//   angle index, checks that the exact ratio R(k)=AnisoRatio(m) never
//   exceeds the stored bound by more than the 1e-12 relative margin
//   test 4g calls for: R(k)*CellSum(actual cell) <= Weight*(1+1e-12),
//   equivalently (dividing by CellSum, skipped when CellSum==0 since
//   both sides are then exactly 0) R(k) <= AcceptanceBound(..)*(1+1e-12).
//
static void CheckAxiBoundsValid(const ScatterParams & sp,
                                const S2::TesselSphere & toas,
                                const std::vector<Real> & thetaDegsP,
                                const std::vector<std::pair<Real,Real> > & thetaChiPairsS,
                                unsigned stride) {

  std::vector<ProbDist> pdists;
  BuildIsoPDists(sp, toas, pdists);
  AxiProposal axi(sp, toas, pdists);

  const Real kP = sp.GetKP();
  const Real kS = sp.GetKS();
  const Real kOut[2] = {kP, kS};  // conv-local 0 -> P out, 1 -> S out

  // P states (intype RAY_P: conv 0=GPP out P, conv1=GPS out S):
  for (std::size_t ti = 0; ti < thetaDegsP.size(); ++ti) {
    const Real theta = thetaDegsP[ti] * Geometry::DtoR;
    const R3::XYZ u_rel(-std::sin(theta), 0.0, std::cos(theta));
    for (Index k = 0; k < toas.size(); k += stride) {
      const Index cell = axi.CellOf(k);
      const R3::XYZ t = toas[k].XYZ();
      for (int li = 0; li < 2; ++li) {
        const int conv = li;  // GPP=0, GPS=1
        const R3::XYZ m = t.ScaledBy(kOut[li]) + R3::XYZ(0,0,-kP);
        const Real ratio = RatioForAxisymmetric(sp, m, u_rel);
        const Real bound = axi.AcceptanceBound(RAY_P, u_rel, conv, cell);
        const Real cellSum = pdists[conv].GetRangeMagnitude(
          cell*(toas.size()/axi.NumCells()), (cell+1)*(toas.size()/axi.NumCells()));
        if (cellSum <= 0.0) continue;
        CHECK_MSG(ratio <= bound*(1.0+1e-12),
                 "P bound violated: theta=" << thetaDegsP[ti]
                 << " k=" << k << " conv=" << conv
                 << " ratio=" << ratio << " bound=" << bound);
      }
    }
  }

  // S states (intype RAY_S: conv 2=GSP out P, conv3=GSS out S):
  for (std::size_t si = 0; si < thetaChiPairsS.size(); ++si) {
    const Real theta = thetaChiPairsS[si].first  * Geometry::DtoR;
    const Real chi   = thetaChiPairsS[si].second * Geometry::DtoR;
    const R3::XYZ u_rel(-std::sin(theta)*std::cos(chi),
                        std::sin(theta)*std::sin(chi),
                        std::cos(theta));
    for (Index k = 0; k < toas.size(); k += stride) {
      const Index cell = axi.CellOf(k);
      const R3::XYZ t = toas[k].XYZ();
      for (int li = 0; li < 2; ++li) {
        const int conv = 2 + li;  // GSP=2, GSS=3
        const R3::XYZ m = t.ScaledBy(kOut[li]) + R3::XYZ(0,0,-kS);
        const Real ratio = RatioForAxisymmetric(sp, m, u_rel);
        const Real bound = axi.AcceptanceBound(RAY_S, u_rel, conv, cell);
        const Real cellSum = pdists[conv].GetRangeMagnitude(
          cell*(toas.size()/axi.NumCells()), (cell+1)*(toas.size()/axi.NumCells()));
        if (cellSum <= 0.0) continue;
        CHECK_MSG(ratio <= bound*(1.0+1e-12),
                 "S bound violated: theta=" << thetaChiPairsS[si].first
                 << " chi=" << thetaChiPairsS[si].second
                 << " k=" << k << " conv=" << conv
                 << " ratio=" << ratio << " bound=" << bound);
      }
    }
  }

}


TEST(axi_bounds_valid) {

  const S2::TesselSphere & toas = UseTOADegree(5);

  struct SpecCase {Real ax,ay,az;};
  const SpecCase specs[3] = { {2,2,0.5}, {5,5,0.5}, {0.5,0.5,2} };
  const Real kappas[2] = {0.3, 0.5};
  const Real freqs[2]  = {1.0, 4.0};

  // P: every bin corner (0,2,...,90) plus a handful past 90 (up to
  // 180) to exercise the z<0 fold.
  std::vector<Real> thetaDegsP;
  for (int b = 0; b <= AxiProposal::kThetaBins; ++b) {thetaDegsP.push_back(2.0*b);}
  const Real extra[] = {100, 135, 150, 170, 180};
  for (std::size_t i = 0; i < sizeof(extra)/sizeof(extra[0]); ++i) {
    thetaDegsP.push_back(extra[i]);
  }

  // S: a 3x3 subset of (theta,chi) bin corners, repeated into all four
  // chi quadrants via the chi value itself (0-360).
  std::vector<std::pair<Real,Real> > thetaChiS;
  const Real thetaCorners[3] = {0.0, 44.0, 90.0};
  const Real chiBase[3] = {0.0, 40.0, 90.0};
  for (int qi = 0; qi < 4; ++qi) {
    const Real qoff = 90.0*qi;
    for (int ti = 0; ti < 3; ++ti) {
      for (int ci = 0; ci < 3; ++ci) {
        Real chi = qoff + ((qi%2==0) ? chiBase[ci] : (90.0-chiBase[ci]));
        thetaChiS.push_back(std::make_pair(thetaCorners[ti], chi));
      }
    }
  }

  for (int si = 0; si < 3; ++si) {
    for (int ki = 0; ki < 2; ++ki) {
      for (int fi = 0; fi < 2; ++fi) {
        ScatterParams::SetFrequencyHertz(freqs[fi]);
        ScatterParams sp(Elastic::VpVs(6.2,3.6),
                        Elastic::HSnexyzk(0.8,0.01,specs[si].ax,
                                          specs[si].ay,specs[si].az,
                                          kappas[ki]));
        CheckAxiBoundsValid(sp, toas, thetaDegsP, thetaChiS, 5);
      }
    }
  }

}


TEST(axi_sampler_rejects_ax_ne_ay) {

  // M4: forcing SAMPLER_AXISYMMETRIC onto a non-axisymmetric (ax!=ay)
  // spec must throw at construction rather than silently build
  // AxiProposal with ah=GetAx() (ignoring ay), which gives an invalid
  // bound for the actual spectrum.
  UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 1, 0.5);  // ax=2 != ay=1
  CHECK(!sp.IsAxisymmetric());
  bool threw = false;
  try {
    Scatterer scat(sp, Scatterer::SAMPLER_AXISYMMETRIC);
    (void) scat;
  } catch (const std::logic_error & e) {
    threw = true;
  }
  CHECK_MSG(threw, "SAMPLER_AXISYMMETRIC on ax!=ay spec should throw "
           "std::logic_error");

  // Sanity: the same spec through SAMPLER_GENERAL (or AUTO) must not be
  // affected by the new guard.
  bool threwGeneral = false;
  try {
    Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);
    (void) scat;
  } catch (...) {
    threwGeneral = true;
  }
  CHECK_MSG(!threwGeneral, "SAMPLER_GENERAL on ax!=ay spec must not throw");

}


TEST(axi_proposal_rejects_non_icosahedral_toa) {

  // M5: AxiProposal's cell/block arithmetic assumes an icosahedral
  // (TESS_ICO, size = 20*4^D) take-off-angle set. A non-icosahedral set
  // must not be silently accepted -- it would make mCellBlock truncate
  // (integer division), dropping the tail directions from every cell's
  // range, a silent bias rather than a loud failure. S2::S2Set is just
  // a std::vector<ThetaPhi> (see geom.hpp's typedef), so build one by
  // hand at a size (100) that isn't a multiple of 20.
  S2::S2Set oddToas;
  for (int i = 0; i < 100; ++i) {
    const Real theta = Geometry::Pi * (i + 0.5) / 100.0;
    const Real phi    = 2.0 * Geometry::Pi * 0.61803398875 * i; // golden-angle spiral
    oddToas.push_back(S2::ThetaPhi(theta, phi));
  }
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  std::vector<ProbDist> pdists;
  BuildIsoPDists(sp, oddToas, pdists);

  bool threw = false;
  try {
    AxiProposal axi(sp, oddToas, pdists);
    (void) axi;
  } catch (const std::logic_error & e) {
    threw = true;
  }
  CHECK_MSG(threw, "AxiProposal on a non-icosahedral TOA set should throw "
           "std::logic_error");

}


TEST(axi_matches_general) {

  UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  Scatterer axiScat(sp, Scatterer::SAMPLER_AXISYMMETRIC);
  Scatterer genScat(sp, Scatterer::SAMPLER_GENERAL);

  const R3::XYZ E(1,0,0), N(0,1,0), U(0,0,1), D(Dir(60,30));
  struct Case {raytype in; R3::XYZ dir; R3::XYZ pol;};
  const Case cases[] = {
    {RAY_P, E, E}, {RAY_P, N, N}, {RAY_P, U, U}, {RAY_P, D, D},
    {RAY_S, E, N}, {RAY_S, E, U},
    {RAY_S, N, E}, {RAY_S, N, U},
    {RAY_S, U, E}, {RAY_S, U, N},
  };
  const int n = 20000;
  int divisor = 1;
  for (int i = 0; i < 5-2; ++i) divisor *= 4;  // 4^(5-2)=64

  std::map<std::pair<int,int>, int> countsA, countsG; // (state*8+conv, coarse)

  for (std::size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i) {
    const S2::ThetaPhi dirTP(cases[i].dir.Theta(), cases[i].dir.Phi());

    const Real meanA = MeanDistanceToRealEvent(axiScat, cases[i].in, dirTP,
                                               cases[i].pol, n);
    const Real meanG = MeanDistanceToRealEvent(genScat, cases[i].in, dirTP,
                                               cases[i].pol, n);
    const Real sigA = meanA/std::sqrt((Real)n);
    const Real sigG = meanG/std::sqrt((Real)n);
    const Real tol = 4.0*std::sqrt(sigA*sigA+sigG*sigG);
    CHECK_MSG(std::fabs(meanA-meanG) < tol,
             "case " << i << ": meanA=" << meanA << " meanG=" << meanG
             << " tol=" << tol);

    const Real rot = Scatterer::RotationForS1(dirTP, cases[i].pol);
    const R3::XYZ origin(0,0,0);
    const int targetReal = 50000;
    int collectedA = 0, collectedG = 0;
    while (collectedA < targetReal || collectedG < targetReal) {
      if (collectedA < targetReal) {
        Scatterer::AnisoAttempt a = axiScat.AttemptAnisotropicScatter(
          cases[i].in, origin, dirTP, rot);
        if (a.IsReal) {
          countsA[std::make_pair((int)i*8+a.Conv, (int)(a.Toa/divisor))]++;
          ++collectedA;
        }
      }
      if (collectedG < targetReal) {
        Scatterer::AnisoAttempt g = genScat.AttemptAnisotropicScatter(
          cases[i].in, origin, dirTP, rot);
        if (g.IsReal) {
          countsG[std::make_pair((int)i*8+g.Conv, (int)(g.Toa/divisor))]++;
          ++collectedG;
        }
      }
    }
  }

  // Two-sample chi-square test for homogeneity between the two
  // samplers' (Conv, degree-2 cell) distributions. The two samplers
  // collect different total real-event counts out of the same number
  // of attempts (different acceptance rates), so the standard
  // unequal-N two-sample statistic is used:
  //   chi2 = sum_i (a_i*Ng - g_i*Na)^2 / (Na*Ng*(a_i+g_i))
  Real Na = 0, Ng = 0;
  for (std::map<std::pair<int,int>,int>::const_iterator it = countsA.begin();
       it != countsA.end(); ++it) {Na += it->second;}
  for (std::map<std::pair<int,int>,int>::const_iterator it = countsG.begin();
       it != countsG.end(); ++it) {Ng += it->second;}

  Real chi2 = 0; int dof = 0;
  std::map<std::pair<int,int>,int> keys = countsA;
  for (std::map<std::pair<int,int>,int>::const_iterator it = countsG.begin();
       it != countsG.end(); ++it) {keys[it->first] = 0;}
  for (std::map<std::pair<int,int>,int>::const_iterator it = keys.begin();
       it != keys.end(); ++it) {
    const Real a = countsA.count(it->first) ? countsA[it->first] : 0;
    const Real g = countsG.count(it->first) ? countsG[it->first] : 0;
    if (a+g < 20.0) continue;
    const Real num = a*Ng - g*Na;
    chi2 += num*num/(Na*Ng*(a+g));
    ++dof;
  }
  CHECK_MSG(dof > 0, "no bins had combined count >= 20");
  const Real chi2_per_dof = chi2/dof;
  CHECK_MSG(chi2_per_dof < 1.5,
           "chi2/dof too large: " << chi2_per_dof << " (dof=" << dof << ")");

}


TEST(axi_draws_per_event) {

  UseTOADegree(6);

  struct SpecLimits {Real ax,ay,az; Real limitP, limitS;};
  const SpecLimits specs[2] = {
    {2,2,0.5, 2.2, 2.6},
    {5,5,0.5, 3.6, 6.0},
  };

  std::ofstream out("build/aniso_draws.txt");
  out << "# spec theta chi raytype sampler mean_draws\n";

  const Real pThetas[4] = {0,30,60,89};
  const Real sThetas[2] = {30,89};
  const Real sChis[3]   = {0,45,90};
  const int n = 3000;

  for (int si = 0; si < 2; ++si) {
    ScatterParams sp = LopNorCrust(specs[si].ax, specs[si].ay, specs[si].az);
    Scatterer axiScat(sp, Scatterer::SAMPLER_AXISYMMETRIC);
    Scatterer genScat(sp, Scatterer::SAMPLER_GENERAL);

    Real maxP_axi = 0, maxS_axi = 0;

    for (int ti = 0; ti < 4; ++ti) {
      const S2::ThetaPhi dir(pThetas[ti]*Geometry::DtoR, 0.0);
      const R3::XYZ pol = dir.XYZ();
      const Real dAxi = MeanDrawsPerRealEvent(axiScat, RAY_P, dir, pol, n);
      const Real dGen = MeanDrawsPerRealEvent(genScat, RAY_P, dir, pol, n);
      out << "(" << specs[si].ax << "," << specs[si].ay << "," << specs[si].az
          << ") " << pThetas[ti] << " - P AXI " << dAxi << "\n";
      out << "(" << specs[si].ax << "," << specs[si].ay << "," << specs[si].az
          << ") " << pThetas[ti] << " - P GENERAL " << dGen << "\n";
      maxP_axi = std::max(maxP_axi, dAxi);
      CHECK_MSG(dAxi < dGen, "P theta=" << pThetas[ti] << ": axi=" << dAxi
               << " general=" << dGen);
    }

    for (int ti = 0; ti < 2; ++ti) {
      for (int ci = 0; ci < 3; ++ci) {
        const S2::ThetaPhi dir(sThetas[ti]*Geometry::DtoR, 0.0);
        const R3::OrthoAxes rotAA(dir.Theta(), dir.Phi(),
                                  sChis[ci]*Geometry::DtoR);
        const R3::XYZ pol = rotAA.S1();
        const Real dAxi = MeanDrawsPerRealEvent(axiScat, RAY_S, dir, pol, n);
        const Real dGen = MeanDrawsPerRealEvent(genScat, RAY_S, dir, pol, n);
        out << "(" << specs[si].ax << "," << specs[si].ay << "," << specs[si].az
            << ") " << sThetas[ti] << " " << sChis[ci] << " S AXI " << dAxi << "\n";
        out << "(" << specs[si].ax << "," << specs[si].ay << "," << specs[si].az
            << ") " << sThetas[ti] << " " << sChis[ci] << " S GENERAL " << dGen << "\n";
        maxS_axi = std::max(maxS_axi, dAxi);
        CHECK_MSG(dAxi < dGen, "S theta=" << sThetas[ti] << " chi=" << sChis[ci]
                 << ": axi=" << dAxi << " general=" << dGen);
      }
    }

    CHECK_MSG(maxP_axi <= specs[si].limitP,
             "spec " << si << ": max P draws " << maxP_axi
             << " exceeds limit " << specs[si].limitP);
    CHECK_MSG(maxS_axi <= specs[si].limitS,
             "spec " << si << ": max S draws " << maxS_axi
             << " exceeds limit " << specs[si].limitS);
  }

  out.close();

}


TEST(axi_mfp_curves) {

  UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  Scatterer scat(sp, Scatterer::SAMPLER_AXISYMMETRIC);
  const S2::TesselSphere dirs6(S2::TESS_ICO, 6);

  std::ofstream out("build/aniso_mfp_curves.txt");
  out << "# theta mfp_P mfp_SV mfp_SH\n";

  Real mfpP0 = 0, mfpP90 = 0;

  for (int theta = 0; theta <= 180; theta += 5) {
    const S2::ThetaPhi dirTP(theta*Geometry::DtoR, 0.0);
    const R3::XYZ dXYZ = dirTP.XYZ();
    const R3::XYZ svPol = dirTP.ThetaHat().XYZ();
    const R3::XYZ shPol = dirTP.PhiHat().XYZ();

    const Real rateP  = scat.IntegratedRate(RAY_P, dXYZ, dXYZ, dirs6);
    const Real rateSV = scat.IntegratedRate(RAY_S, dXYZ, svPol, dirs6);
    const Real rateSH = scat.IntegratedRate(RAY_S, dXYZ, shPol, dirs6);
    const Real mfpP = 1.0/rateP, mfpSV = 1.0/rateSV, mfpSH = 1.0/rateSH;

    CHECK(std::isfinite(mfpP)  && mfpP  > 0.0);
    CHECK(std::isfinite(mfpSV) && mfpSV > 0.0);
    CHECK(std::isfinite(mfpSH) && mfpSH > 0.0);

    out << theta << " " << mfpP << " " << mfpSV << " " << mfpSH << "\n";

    if (theta == 0)  mfpP0  = mfpP;
    if (theta == 90) mfpP90 = mfpP;
  }
  out.close();

  CHECK_MSG(mfpP90 < mfpP0, "P(90)=" << mfpP90 << " should be < P(0)=" << mfpP0);

}


TEST(axi_vertical_incidence) {

  UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  Scatterer scat(sp, Scatterer::SAMPLER_AXISYMMETRIC);

  const S2::ThetaPhi dir(0.0, 0.0);
  const R3::XYZ origin(0,0,0);

  for (int i = 0; i < 1000; ++i) {
    Scatterer::AnisoAttempt att =
      scat.AttemptAnisotropicScatter(RAY_P, origin, dir, 0.37);
    CHECK(std::isfinite(att.InPol));
    if (att.IsReal) {
      CHECK(std::isfinite(att.Rel.GetDirection().Theta()));
      CHECK(std::isfinite(att.Rel.GetDirection().Phi()));
    }
  }
  for (int i = 0; i < 1000; ++i) {
    Scatterer::AnisoAttempt att =
      scat.AttemptAnisotropicScatter(RAY_S, origin, dir, 0.37);
    CHECK(std::isfinite(att.InPol));
    if (att.IsReal) {
      CHECK(std::isfinite(att.Rel.GetDirection().Theta()));
      CHECK(std::isfinite(att.Rel.GetDirection().Phi()));
      CHECK(std::isfinite(att.Rel.GetPolarization()));
    }
  }

}


TEST(axi_curved_frame) {

  ECS.SetMapping(EarthCoords::RAE_CURVED);

  const R3::XYZ loc(950.0, 0.0, 0.0);
  const R3::XYZ up = ECS.GetUp(loc);

  const R3::XYZ dXYZ = Dir(70, 20);
  const S2::ThetaPhi dir(dXYZ.Theta(), dXYZ.Phi());

  const Real rot = Scatterer::CanonicalPPolarization(dir, up);
  const R3::OrthoAxes AA(dir.Theta(), dir.Phi(), rot);

  const R3::XYZ d = dir.XYZ();
  const R3::XYZ target = d.ScaledBy(up.Dot(d)) + up.ScaledBy(-1.0);
  const R3::XYZ targetUnit = target.Unit();
  const R3::XYZ s1 = AA.S1();

  CHECK(s1.Cross(targetUnit).Mag() < 1e-12);
  CHECK(s1.Dot(targetUnit) > 0.0);

  ECS.SetMapping(EarthCoords::ENU_ORTHO);

}


TEST(axi_low_degree) {

  const S2::TesselSphere & toas = UseTOADegree(2);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);

  std::vector<ProbDist> pdists;
  BuildIsoPDists(sp, toas, pdists);
  AxiProposal axi(sp, toas, pdists);
  CHECK_MSG(axi.NumCells() == 320, "expected 320 cells, got " << axi.NumCells());

  std::vector<Real> thetaDegsP;
  for (Real t = 0; t <= 180; t += 15) {thetaDegsP.push_back(t);}
  std::vector<std::pair<Real,Real> > thetaChiS;
  for (Real t = 0; t <= 90; t += 30) {
    for (Real c = 0; c <= 90; c += 30) {
      thetaChiS.push_back(std::make_pair(t,c));
    }
  }
  CheckAxiBoundsValid(sp, toas, thetaDegsP, thetaChiS, 1);

  Scatterer scat(sp, Scatterer::SAMPLER_AXISYMMETRIC);
  const S2::ThetaPhi dir(37.0*Geometry::DtoR, 11.0*Geometry::DtoR);
  const R3::XYZ origin(0,0,0);
  int reals = 0;
  for (int i = 0; i < 1000 && reals < 1000; ) {
    Scatterer::AnisoAttempt att =
      scat.AttemptAnisotropicScatter(RAY_P, origin, dir, 0.0);
    ++i;
    if (att.IsReal) {++reals;}
  }
  CHECK_MSG(reals > 0, "no real events observed at low degree");

}


TEST(axi_eps_zero) {

  UseTOADegree(5);
  ScatterParams sp = LopNorCrust(2, 2, 0.5, 0.5, 0.0);  // eps=0
  bool threw = false;
  Scatterer * scat = 0;
  try {
    scat = new Scatterer(sp, Scatterer::SAMPLER_AXISYMMETRIC);
  } catch (...) {
    threw = true;
  }
  CHECK_MSG(!threw, "construction with eps=0 must not throw");
  const Real mfpP = scat->GetMajorantMeanFreePath(RAY_P);
  const Real mfpS = scat->GetMajorantMeanFreePath(RAY_S);
  CHECK(std::isinf(mfpP) && mfpP > 0.0);
  CHECK(std::isinf(mfpS) && mfpS > 0.0);

  // I3: no NaN in the reported directional MFPs either (see
  // general_eps_zero's comment -- same AnisoRatio() fix, exercised here
  // via IntegrateRateAndDipole() on the axisymmetric sampler's path).
  const std::vector<Real> & mfps = scat->GetReportMFPs();
  CHECK_MSG(mfps.size() == 8, "expected 8 AXI directional MFPs, got "
           << mfps.size());
  for (std::size_t i = 0; i < mfps.size(); ++i) {
    CHECK_MSG(!std::isnan(mfps[i]), "GetReportMFPs()[" << i << "] is NaN");
    CHECK_MSG(mfps[i] > 0.0, "GetReportMFPs()[" << i << "] is not positive");
  }
  delete scat;

}


TEST(axi_extreme_ratio) {

  const S2::TesselSphere & toas = UseTOADegree(4);
  ScatterParams sp = LopNorCrust(10, 10, 0.5);

  std::vector<ProbDist> pdists;
  BuildIsoPDists(sp, toas, pdists);
  AxiProposal axi(sp, toas, pdists);  // must not throw

  std::vector<Real> thetaDegsP;
  for (Real t = 0; t <= 180; t += 20) {thetaDegsP.push_back(t);}
  std::vector<std::pair<Real,Real> > thetaChiS;
  for (Real t = 0; t <= 90; t += 30) {
    for (Real c = 0; c <= 90; c += 30) {
      thetaChiS.push_back(std::make_pair(t,c));
    }
  }
  CheckAxiBoundsValid(sp, toas, thetaDegsP, thetaChiS, 3);

  Scatterer scat(sp, Scatterer::SAMPLER_AXISYMMETRIC);
  const S2::ThetaPhi dirP(60.0*Geometry::DtoR, 0.0);
  const Real meanP = MeanDrawsPerRealEvent(scat, RAY_P, dirP, dirP.XYZ(), 500);
  CHECK_MSG(std::isfinite(meanP) && meanP < 50.0,
           "mean P draws/event=" << meanP);

  const S2::ThetaPhi dirS(60.0*Geometry::DtoR, 0.0);
  const R3::OrthoAxes rotAA(dirS.Theta(), dirS.Phi(), 45.0*Geometry::DtoR);
  const Real meanS = MeanDrawsPerRealEvent(scat, RAY_S, dirS, rotAA.S1(), 500);
  CHECK_MSG(std::isfinite(meanS) && meanS < 50.0,
           "mean S draws/event=" << meanS);

}


TEST(isolated_axi_overridemfp) {

  UseTOADegree(5);
  Scatterer::OverrideMFP(100, 50);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  Scatterer axiScat(sp, Scatterer::SAMPLER_AXISYMMETRIC);
  Scatterer genScat(sp, Scatterer::SAMPLER_GENERAL);

  const int n = 20000;
  Real total = 0;
  for (int i = 0; i < n; ++i) {total += axiScat.GetRandomPathLength(RAY_P);}
  const Real mean = total/n;
  CHECK(NearRel(mean, 100.0, 0.03));

  const R3::XYZ origin(0,0,0);
  const S2::ThetaPhi dirP(60.0*Geometry::DtoR, 30.0*Geometry::DtoR);
  int divisor = 1;
  for (int i = 0; i < 5-2; ++i) divisor *= 4;

  std::map<int,int> countsA, countsG;
  for (int i = 0; i < 20000; ++i) {
    Scatterer::AnisoAttempt a =
      axiScat.AttemptAnisotropicScatter(RAY_P, origin, dirP, 0.0);
    CHECK(a.IsReal);
    countsA[a.Conv*100000 + (int)(a.Toa/divisor)]++;

    Scatterer::AnisoAttempt g =
      genScat.AttemptAnisotropicScatter(RAY_P, origin, dirP, 0.0);
    CHECK(g.IsReal);
    countsG[g.Conv*100000 + (int)(g.Toa/divisor)]++;
  }

  std::map<int,int> keys = countsA;
  for (std::map<int,int>::const_iterator it = countsG.begin();
       it != countsG.end(); ++it) {keys[it->first] = 0;}
  Real chi2 = 0; int dof = 0;
  for (std::map<int,int>::const_iterator it = keys.begin();
       it != keys.end(); ++it) {
    const Real a = countsA.count(it->first) ? countsA[it->first] : 0;
    const Real g = countsG.count(it->first) ? countsG[it->first] : 0;
    if (a+g < 20.0) continue;
    const Real num = a-g;
    chi2 += num*num/(a+g);
    ++dof;
  }
  CHECK_MSG(dof > 0, "no bins had combined count >= 20");
  const Real chi2_per_dof = chi2/dof;
  CHECK_MSG(chi2_per_dof < 1.5,
           "chi2/dof too large: " << chi2_per_dof << " (dof=" << dof << ")");

}


// ::: Fix round 1, controller ruling R11 :::
//
//   --overridemfp + eps==0 on an anisotropic Scatterer must throw
//   std::logic_error at construction (both samplers), rather than the
//   axisymmetric sampler's accept/reject step computing a silent 0/0
//   NaN that can't trip the bound-violation check and hangs
//   --overridemfp's "loop until accepted" forever (every AxiProposal
//   bin weight, and every general-sampler mPDists entry, is identically
//   0 when eps==0 -- see ScatterParams::PSATO()'s eps^2 prefactor).
//   Each needs its own fresh process (Scatterer::OverrideMFP() is
//   class-static and unresettable within one run), hence isolated_*.
//

TEST(isolated_axi_overridemfp_eps_zero_throws) {

  UseTOADegree(5);
  Scatterer::OverrideMFP(100, 50);
  ScatterParams sp = LopNorCrust(2, 2, 0.5, 0.5, 0.0);  // eps=0

  bool threw = false;
  std::string what;
  try {
    Scatterer scat(sp, Scatterer::SAMPLER_AXISYMMETRIC);
  } catch (const std::logic_error & e) {
    threw = true;
    what = e.what();
  }
  CHECK_MSG(threw, "SAMPLER_AXISYMMETRIC construction under "
           "--overridemfp + eps=0 should throw std::logic_error");
  CHECK_MSG(what.find("eps") != std::string::npos,
           "exception message should mention eps: " << what);

}


TEST(isolated_general_overridemfp_eps_zero_throws) {

  UseTOADegree(5);
  Scatterer::OverrideMFP(100, 50);
  ScatterParams sp = LopNorCrust(2, 2, 0.5, 0.5, 0.0);  // eps=0

  bool threw = false;
  std::string what;
  try {
    Scatterer scat(sp, Scatterer::SAMPLER_GENERAL);
  } catch (const std::logic_error & e) {
    threw = true;
    what = e.what();
  }
  CHECK_MSG(threw, "SAMPLER_GENERAL construction under "
           "--overridemfp + eps=0 should throw std::logic_error");
  CHECK_MSG(what.find("eps") != std::string::npos,
           "exception message should mention eps: " << what);

}


TEST(isolated_isotropic_overridemfp_eps_zero_unaffected) {

  // Controller ruling R11 only targets anisotropic scatterers:
  // isotropic behavior (today's exact path) must stay untouched, even
  // under --overridemfp with eps=0.
  UseTOADegree(5);
  Scatterer::OverrideMFP(100, 50);
  ScatterParams sp = LopNorCrust(0.5, 0.5, 0.5, 0.5, 0.0);  // isotropic, eps=0

  bool threw = false;
  Scatterer * scat = 0;
  try {
    scat = new Scatterer(sp, Scatterer::SAMPLER_AUTO);
  } catch (...) {
    threw = true;
  }
  CHECK_MSG(!threw, "isotropic construction under --overridemfp + eps=0 "
           "must not throw");
  CHECK(!scat->IsAnisotropic());
  CHECK(NearRel(scat->GetMajorantMeanFreePath(RAY_P), 100.0, 1e-12));
  delete scat;

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 2 (PERFORMANCE PLAN) TESTS: Bound-kernel restructure   ****
// ****  and geometry cache                                          ****
// ****                                                              ****
//
// See docs/superpowers/specs/2026-10-17-axisymmetric-sampler-
// performance-design.md, decision 1, and
// .superpowers/sdd/2026-10-17-axisymmetric-sampler-performance/
// task-2-brief.md.
//

// ::: ReferenceBound3db253a(...) :::
//
//   Ruling R-P1's frozen oracle: a free-standing, literal transcription
//   of AxiProposal::ComputeBound()'s body exactly as it existed at
//   commit 3db253a (the commit this branch -- and the performance plan
//   -- starts from). Confirmed identical to aniso_axi.cpp at this
//   branch's HEAD (a8ee04b) by `git diff 3db253a HEAD -- aniso_axi.cpp`
//   (empty) before this task began, so freezing "3db253a's ComputeBound"
//   and "today's ComputeBound" are the same thing as of this writing.
//   This is the oracle test_axi_kernel_matches_reference checks the
//   restructured PrepareBoundRep()/BoundRatio() kernel against: every
//   member-variable read ('this->mAh' etc.) is replaced by a plain
//   parameter of the same name, nothing else is changed. Do not "fix"
//   or simplify this function if it starts looking stale -- it must
//   stay a literal transcription of the frozen commit, bugs and all.
//
static Real ReferenceBound3db253a(const R3::XYZ & t, Real r_i,
                                  Real k_in, Real k_out,
                                  const R3::XYZ & u_c, Real delta,
                                  Real ah, Real az, Real aMin,
                                  Real M, Real kappa) {

  const R3::XYZ zhat(0.0, 0.0, 1.0);
  const R3::XYZ m_rep = t.ScaledBy(k_out) + zhat.ScaledBy(-k_in);
  const Real mr  = m_rep.Mag();
  const Real rho = k_out * 2.0 * sin(r_i/2.0);

  Real alpha_c = 0.0;
  if (mr > 0.0) {
    Real cosang = m_rep.Unit().Dot(u_c);
    if (cosang > 1.0)  cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;
    alpha_c = acos(cosang);
  }
  const Real alpha_lo = std::max((Real)0.0, alpha_c - delta);
  const Real alpha_hi = std::min((Real)Geometry::Pi, alpha_c + delta);

  const Real s_lo = mr*cos(alpha_hi) - rho;
  const Real s_hi = mr*cos(alpha_lo) + rho;

  const Real s2_hi = std::max(s_lo*s_lo, s_hi*s_hi);
  const Real s2_lo = (s_lo <= 0.0 && 0.0 <= s_hi) ?
                      (Real)0.0 : std::min(s_lo*s_lo, s_hi*s_hi);

  const Real mlo   = std::max(mr-rho, (Real)0.0);
  const Real m2_lo = mlo*mlo;
  const Real mhi   = mr+rho;
  const Real m2_hi = mhi*mhi;

  Real Q_lo;
  if (ah >= az) {
    Q_lo = std::max(ah*ah*m2_lo - (ah*ah-az*az)*std::min(s2_hi,m2_hi),
                    az*az*m2_lo);
  } else {
    Q_lo = ah*ah*m2_lo + (az*az-ah*ah)*s2_lo;
  }

  const Real ratioBase = (1.0 + aMin*aMin*m2_hi) / (1.0 + Q_lo);
  const Real bound = M * pow(ratioBase, kappa+1.5);

  return std::min(M, bound);

}


TEST(axi_kernel_matches_reference) {

  // Representative directions: every TESS_ICO degree-3 centre (1,280 of
  // them), plus zhat on its own -- t==zhat with k_in==k_out makes
  // m_rep (== k_out*t - k_in*zhat) exactly the zero vector, exercising
  // PrepareBoundRep()'s "MHat = zero vector if Mr == 0" special case
  // (today's ComputeBound() never calls Unit() on m_rep in that case
  // either; see the 'mr > 0.0' guard above).
  std::vector<R3::XYZ> ts;
  {
    S2::TesselSphere tess(S2::TESS_ICO, 3);
    for (Index i = 0; i < (Index) tess.size(); ++i) {ts.push_back(tess[i].XYZ());}
  }
  ts.push_back(R3::XYZ(0.0, 0.0, 1.0));

  const Real r_is[3] = {0.0, 0.003, 0.05};

  struct KPair {Real k_in, k_out;};
  const KPair kpairs[4] = { {1.0,1.0}, {1.0,1.73}, {1.73,1.0}, {1.73,1.73} };

  // Incidence-bin caps, built with exactly the u_c/delta expressions
  // the AxiProposal constructor uses for P bins {0,22,44} (bin corners,
  // mid-range, and the 90-degree edge) and S bins {(0,0),(22,4),(44,8)}
  // (same theta corners, spread across chi).
  struct CapCase {R3::XYZ u_c; Real delta;};
  std::vector<CapCase> caps;
  const int pBins[3] = {0, 22, 44};
  for (int i = 0; i < 3; ++i) {
    const int b = pBins[i];
    const Real theta_lo = (2.0*b)     * Geometry::DtoR;
    const Real theta_hi = (2.0*b+2.0) * Geometry::DtoR;
    const Real theta_c  = 0.5*(theta_lo+theta_hi);
    const Real delta    = 0.5*(theta_hi-theta_lo);
    CapCase cc;
    cc.u_c = R3::XYZ(-sin(theta_c), 0.0, cos(theta_c));
    cc.delta = delta;
    caps.push_back(cc);
  }
  const int sBC[3][2] = {{0,0}, {22,4}, {44,8}};
  for (int i = 0; i < 3; ++i) {
    const int b = sBC[i][0], c = sBC[i][1];
    const Real theta_lo = (2.0*b)     * Geometry::DtoR;
    const Real theta_hi = (2.0*b+2.0) * Geometry::DtoR;
    const Real theta_c  = 0.5*(theta_lo+theta_hi);
    const Real chi_lo = (10.0*c)      * Geometry::DtoR;
    const Real chi_hi = (10.0*c+10.0) * Geometry::DtoR;
    const Real chi_c  = 0.5*(chi_lo+chi_hi);
    const Real delta = 0.5*(theta_hi-theta_lo)
                      + sin(theta_hi)*0.5*(chi_hi-chi_lo);
    CapCase cc;
    cc.u_c = R3::XYZ(-sin(theta_c)*cos(chi_c),
                     sin(theta_c)*sin(chi_c),
                     cos(theta_c));
    cc.delta = delta;
    caps.push_back(cc);
  }

  struct Spec {Real ah, az, kappa;};
  const Spec specs[4] = {
    {2.0, 0.5, 0.3}, {2.0, 0.5, 0.5}, {0.5, 2.0, 0.5}, {10.0, 0.5, 0.5}
  };

  long long count = 0;
  for (std::size_t ti = 0; ti < ts.size(); ++ti) {
    for (int ri = 0; ri < 3; ++ri) {
      for (int ki = 0; ki < 4; ++ki) {
        for (std::size_t ci = 0; ci < caps.size(); ++ci) {
          for (int si = 0; si < 4; ++si) {

            const Real aMin = std::min(specs[si].ah, specs[si].az);
            const Real M = specs[si].ah*specs[si].ah*specs[si].az
                         / (aMin*aMin*aMin);

            const AxiProposal::BoundRep rep =
              AxiProposal::PrepareBoundRep(ts[ti], r_is[ri],
                                           kpairs[ki].k_in, kpairs[ki].k_out,
                                           aMin);
            const AxiProposal::BinCap cap =
              AxiProposal::MakeBinCap(caps[ci].u_c, caps[ci].delta);
            const Real ratio = AxiProposal::BoundRatio(rep, cap,
                                                        specs[si].ah,
                                                        specs[si].az);
            const Real got = std::min(M, M*pow(ratio, specs[si].kappa+1.5));

            const Real want = ReferenceBound3db253a(
              ts[ti], r_is[ri], kpairs[ki].k_in, kpairs[ki].k_out,
              caps[ci].u_c, caps[ci].delta,
              specs[si].ah, specs[si].az, aMin, M, specs[si].kappa);

            CHECK_MSG(std::fabs(got-want) <= 1e-12*std::fabs(want),
                     "mismatch at t#" << ti << " r_i=" << r_is[ri]
                     << " k_in=" << kpairs[ki].k_in
                     << " k_out=" << kpairs[ki].k_out
                     << " cap#" << ci << " spec#" << si
                     << ": got=" << got << " want=" << want);
            ++count;
          }
        }
      }
    }
  }

  // Edge cases for the trig-free rewrite (design decision 3): caps at
  // and across the c=+1/c=-1 boundaries, where the closed-form cos(alpha
  // -+ delta) identities replace acos/cos and a catastrophic-cancellation
  // regression (near c = +-1) would first show up. Checked against the
  // same acos/cos oracle above, which is untouched by the rewrite.
  //
  //   - u_c == +MHat (c=+1 exactly) and u_c == -MHat (c=-1 exactly);
  //   - u_c rotated 0.5*delta off +MHat (alpha_c == 0.5*delta, just
  //     inside the "within delta of 0" branch) and off -MHat (alpha_c ==
  //     pi-0.5*delta, just inside "within delta of pi");
  //   - the same four cap directions again, against the Mr==0
  //     representative (t=zhat, k_in=k_out=1.73), which must return the
  //     same value for every one of them (PrepareBoundRep()/BoundRatio()
  //     force c=1 whenever Mr==0, regardless of u_c).
  //
  // t=(0.6,0,0.8) (unit), k_in=1, k_out=1.73 gives a representative
  // M-vector that is neither zero nor aligned with any axis, so +-MHat
  // and the rotated directions below are non-trivial unit vectors.
  {
    const R3::XYZ t_edge(0.6, 0.0, 0.8);
    const Real k_in_edge = 1.0, k_out_edge = 1.73;
    const R3::XYZ zhat_edge(0.0, 0.0, 1.0);
    const R3::XYZ mhat_edge =
      AxiProposal::PrepareBoundRep(t_edge, 0.0, k_in_edge, k_out_edge, 1.0).MHat;
    const R3::XYZ perp_edge = mhat_edge.GetInPlaneUnitPerpendicular(zhat_edge);

    std::vector<CapCase> edgeCaps;
    for (std::size_t ci = 0; ci < caps.size(); ++ci) {
      const Real delta = caps[ci].delta;
      const Real half  = 0.5*delta;
      CapCase cc;

      cc.u_c = mhat_edge;                       // c == +1
      cc.delta = delta;
      edgeCaps.push_back(cc);

      cc.u_c = mhat_edge.ScaledBy(-1.0);        // c == -1
      cc.delta = delta;
      edgeCaps.push_back(cc);

      cc.u_c = mhat_edge.ScaledBy(cos(half))
             + perp_edge.ScaledBy(sin(half));   // alpha_c == 0.5*delta
      cc.delta = delta;
      edgeCaps.push_back(cc);

      cc.u_c = mhat_edge.ScaledBy(-cos(half))
             + perp_edge.ScaledBy(sin(half));   // alpha_c == pi-0.5*delta
      cc.delta = delta;
      edgeCaps.push_back(cc);
    }

    for (std::size_t ei = 0; ei < edgeCaps.size(); ++ei) {
      const AxiProposal::BinCap cap =
        AxiProposal::MakeBinCap(edgeCaps[ei].u_c, edgeCaps[ei].delta);
      for (int ri = 0; ri < 3; ++ri) {
        for (int si = 0; si < 4; ++si) {

          const Real aMin = std::min(specs[si].ah, specs[si].az);
          const Real M = specs[si].ah*specs[si].ah*specs[si].az
                       / (aMin*aMin*aMin);

          // Mr > 0 (the t_edge/k_in_edge/k_out_edge representative):
          const AxiProposal::BoundRep rep =
            AxiProposal::PrepareBoundRep(t_edge, r_is[ri],
                                         k_in_edge, k_out_edge, aMin);
          const Real ratio =
            AxiProposal::BoundRatio(rep, cap, specs[si].ah, specs[si].az);
          const Real got = std::min(M, M*pow(ratio, specs[si].kappa+1.5));
          const Real want = ReferenceBound3db253a(
            t_edge, r_is[ri], k_in_edge, k_out_edge,
            edgeCaps[ei].u_c, edgeCaps[ei].delta,
            specs[si].ah, specs[si].az, aMin, M, specs[si].kappa);
          CHECK_MSG(std::fabs(got-want) <= 1e-12*std::fabs(want),
                   "edge(Mr>0) mismatch at edge#" << ei << " r_i=" << r_is[ri]
                   << " spec#" << si << ": got=" << got << " want=" << want);
          ++count;

          // Mr == 0 (t=zhat, k_in=k_out=1.73), same cap:
          const AxiProposal::BoundRep rep0 =
            AxiProposal::PrepareBoundRep(zhat_edge, r_is[ri], 1.73, 1.73, aMin);
          const Real ratio0 =
            AxiProposal::BoundRatio(rep0, cap, specs[si].ah, specs[si].az);
          const Real got0 = std::min(M, M*pow(ratio0, specs[si].kappa+1.5));
          const Real want0 = ReferenceBound3db253a(
            zhat_edge, r_is[ri], 1.73, 1.73,
            edgeCaps[ei].u_c, edgeCaps[ei].delta,
            specs[si].ah, specs[si].az, aMin, M, specs[si].kappa);
          CHECK_MSG(std::fabs(got0-want0) <= 1e-12*std::fabs(want0),
                   "edge(Mr==0) mismatch at edge#" << ei << " r_i=" << r_is[ri]
                   << " spec#" << si << ": got=" << got0 << " want=" << want0);
          ++count;
        }
      }
    }
  }

  CHECK(count > 0);

}


TEST(axi_geometry_cache_degrees) {

  // Builds AxiProposal at degree 5, then degree 2 (a different cache
  // key, interleaved to make sure it can't corrupt the degree-5 entry),
  // then degree 5 again. UseTOADegree() caches one TesselSphere per
  // degree for the program's life, so the two degree-5 builds below are
  // handed the exact same take-off-angle set both times -- this exists
  // to catch a geometry cache that corrupts or wrongly shares state
  // across degrees, not to prove determinism (a pure, uncached rebuild
  // would also be bit-identical to itself).
  ScatterParams sp = LopNorCrust(2, 2, 0.5);

  const S2::TesselSphere & toas5a = UseTOADegree(5);
  std::vector<ProbDist> pdists5a;
  BuildIsoPDists(sp, toas5a, pdists5a);
  AxiProposal axi5a(sp, toas5a, pdists5a);

  const S2::TesselSphere & toas2 = UseTOADegree(2);
  std::vector<ProbDist> pdists2;
  BuildIsoPDists(sp, toas2, pdists2);
  AxiProposal axi2(sp, toas2, pdists2);
  (void) axi2;

  const S2::TesselSphere & toas5b = UseTOADegree(5);
  std::vector<ProbDist> pdists5b;
  BuildIsoPDists(sp, toas5b, pdists5b);
  AxiProposal axi5b(sp, toas5b, pdists5b);

  CHECK(axi5a.MajorantRate(RAY_P) == axi5b.MajorantRate(RAY_P));
  CHECK(axi5a.MajorantRate(RAY_S) == axi5b.MajorantRate(RAY_S));
  CHECK(axi5a.NumCells() == axi5b.NumCells());

  for (int b = 0; b < AxiProposal::kThetaBins; ++b) {
    const Real theta = (2.0*b + 1.0) * Geometry::DtoR;  // bin centre
    const R3::XYZ u_c(-sin(theta), 0.0, cos(theta));
    for (int conv = 0; conv < 2; ++conv) {  // GPP, GPS
      for (Index cell = 0; cell < axi5a.NumCells(); ++cell) {
        const Real va = axi5a.AcceptanceBound(RAY_P, u_c, conv, cell);
        const Real vb = axi5b.AcceptanceBound(RAY_P, u_c, conv, cell);
        CHECK_MSG(va == vb,
                 "bin=" << b << " conv=" << conv << " cell=" << cell
                 << ": first-build=" << va << " second-build=" << vb);
      }
    }
  }

}


TEST(axi_geometry_cache_mismatch) {

  // Exercises the fingerprint-MISMATCH branch of GetOrBuildAxiGeometry()
  // (aniso_axi.cpp): a take-off-angle set that shares degree 5's cache
  // key (toa.size()) with the standard set, but isn't it. Per the brief
  // ("On a mismatch the geometry is rebuilt without caching"), this must
  // rebuild geometry for that one construction *without* disturbing the
  // entry already cached under that size -- checked here by sandwiching
  // the mismatched build between two standard-set builds and requiring
  // the first and third to agree exactly (if the mismatch path instead
  // overwrote the cached entry, or read/wrote it concurrently with a
  // corrupted size, the third build could come back wrong) -- and by
  // requiring the mismatched build itself not to throw.
  //
  // Perturbation choice: each direction's phi is shifted by
  // (index mod 7)*1e-4 radians. A rigid rotation about z (this review's
  // first suggestion) was tried and rejected: a standalone check showed
  // TESS_ICO's vertex directions sum to (numerically) zero by
  // icosahedral symmetry, and rotation is linear, so a rotated copy's
  // component sum is *also* (numerically) zero -- a 0.01 rad z-rotation
  // of this exact degree-5 set changed FingerprintOf()'s sum by only
  // ~5e-13, which is sub-floating-point-noise, not a dependable
  // mismatch. Rotation does *not* break mirror matching, for what it's
  // worth: GetOrBuildAxiGeometry() builds the mirror-cell permutations
  // from a fresh canonical TesselSphere, never from 'toa' itself. The
  // index-dependent phi shift used here is not a rigid rotation, so it
  // doesn't inherit that zero-sum cancellation; a standalone check
  // confirmed it changes the fingerprint by about 2e-3, many orders of
  // magnitude above floating-point noise.
  ScatterParams sp = LopNorCrust(2, 2, 0.5);

  const S2::TesselSphere & toas5 = UseTOADegree(5);
  std::vector<ProbDist> pdists5a;
  BuildIsoPDists(sp, toas5, pdists5a);
  AxiProposal axi5a(sp, toas5, pdists5a);

  S2::S2Set perturbed;
  for (Index i = 0; i < (Index) toas5.size(); ++i) {
    const Real shift = 1.0e-4 * (Real) (i % 7);
    perturbed.push_back(S2::ThetaPhi(toas5[i].Theta(), toas5[i].Phi() + shift));
  }
  CHECK(perturbed.size() == toas5.size());

  std::vector<ProbDist> pdistsPerturbed;
  BuildIsoPDists(sp, perturbed, pdistsPerturbed);
  bool threwOnPerturbed = false;
  try {
    AxiProposal axiPerturbed(sp, perturbed, pdistsPerturbed);
    (void) axiPerturbed;
  } catch (...) {
    threwOnPerturbed = true;
  }
  CHECK_MSG(!threwOnPerturbed,
           "AxiProposal on a same-size, perturbed take-off-angle set "
           "should not throw");

  std::vector<ProbDist> pdists5b;
  BuildIsoPDists(sp, toas5, pdists5b);
  AxiProposal axi5b(sp, toas5, pdists5b);

  CHECK(axi5a.MajorantRate(RAY_P) == axi5b.MajorantRate(RAY_P));
  CHECK(axi5a.MajorantRate(RAY_S) == axi5b.MajorantRate(RAY_S));
  CHECK(axi5a.NumCells() == axi5b.NumCells());

  for (int b = 0; b < AxiProposal::kThetaBins; ++b) {
    const Real theta = (2.0*b + 1.0) * Geometry::DtoR;  // bin centre
    const R3::XYZ u_c(-sin(theta), 0.0, cos(theta));
    for (int conv = 0; conv < 2; ++conv) {  // GPP, GPS
      for (Index cell = 0; cell < axi5a.NumCells(); ++cell) {
        const Real va = axi5a.AcceptanceBound(RAY_P, u_c, conv, cell);
        const Real vb = axi5b.AcceptanceBound(RAY_P, u_c, conv, cell);
        CHECK_MSG(va == vb,
                 "bin=" << b << " conv=" << conv << " cell=" << cell
                 << ": first-build=" << va << " second-build=" << vb);
      }
    }
  }

}


TEST(axi_geometry_cache_mismatch_rotation) {

  // Controller ruling R2: FingerprintOf() (aniso_axi.cpp) now uses an
  // index-weighted sum specifically so a rotated take-off-angle set is
  // recognized as different from the un-rotated one -- a plain
  // (unweighted) component sum could not do this reliably, since
  // TESS_ICO's vertex directions sum to (numerically) zero by
  // icosahedral symmetry and rotation preserves that exactly. This test
  // sandwiches a 0.01 rad-about-z rotation of the standard degree-5 set
  // between two standard-set builds, as test_axi_geometry_cache_mismatch
  // does for its own (non-rotation) perturbation, but additionally
  // confirms the rotated build actually took the mismatch path via an
  // observable effect: its AcceptanceBound table must differ from the
  // standard build's *somewhere*, which can only happen if the mismatch
  // path really rebuilt geometry from the rotated directions rather than
  // (bug) silently reusing the standard set's cached geometry.
  ScatterParams sp = LopNorCrust(2, 2, 0.5);

  const S2::TesselSphere & toas5 = UseTOADegree(5);
  std::vector<ProbDist> pdists5a;
  BuildIsoPDists(sp, toas5, pdists5a);
  AxiProposal axi5a(sp, toas5, pdists5a);

  S2::S2Set rotated;
  const Real phi0 = 0.01;  // small fixed rotation about z, radians
  for (Index i = 0; i < (Index) toas5.size(); ++i) {
    rotated.push_back(S2::ThetaPhi(toas5[i].Theta(), toas5[i].Phi() + phi0));
  }
  CHECK(rotated.size() == toas5.size());

  std::vector<ProbDist> pdistsRotated;
  BuildIsoPDists(sp, rotated, pdistsRotated);
  bool threwOnRotated = false;
  AxiProposal * axiRotated = 0;
  try {
    axiRotated = new AxiProposal(sp, rotated, pdistsRotated);
  } catch (...) {
    threwOnRotated = true;
  }
  CHECK_MSG(!threwOnRotated,
           "AxiProposal on a same-size, rotated take-off-angle set "
           "should not throw");

  std::vector<ProbDist> pdists5b;
  BuildIsoPDists(sp, toas5, pdists5b);
  AxiProposal axi5b(sp, toas5, pdists5b);

  // (1) Observable effect proving the mismatch path really rebuilt
  // geometry from the rotated set.
  bool foundDifference = false;
  Real maxAbsDiff = 0.0;
  for (int b = 0; b < AxiProposal::kThetaBins; ++b) {
    const Real theta = (2.0*b + 1.0) * Geometry::DtoR;
    const R3::XYZ u_c(-sin(theta), 0.0, cos(theta));
    for (int conv = 0; conv < 2; ++conv) {
      for (Index cell = 0; cell < axi5a.NumCells(); ++cell) {
        const Real vStd = axi5a.AcceptanceBound(RAY_P, u_c, conv, cell);
        const Real vRot = axiRotated->AcceptanceBound(RAY_P, u_c, conv, cell);
        if (vRot != vStd) {foundDifference = true;}
        maxAbsDiff = std::max(maxAbsDiff, std::fabs(vRot - vStd));
      }
    }
  }
  CHECK_MSG(foundDifference,
           "rotated build's AcceptanceBound table is bit-identical to "
           "the standard build's at every P bin x conv x cell checked -- "
           "the mismatch path may not be rebuilding geometry from the "
           "rotated set (max |diff| observed = " << maxAbsDiff << ")");
  delete axiRotated;

  // (2) The standard set's cache entry survived the rotated build
  // untouched: a fresh standard-set build still matches the first
  // exactly.
  CHECK(axi5a.MajorantRate(RAY_P) == axi5b.MajorantRate(RAY_P));
  CHECK(axi5a.MajorantRate(RAY_S) == axi5b.MajorantRate(RAY_S));
  CHECK(axi5a.NumCells() == axi5b.NumCells());

  for (int b = 0; b < AxiProposal::kThetaBins; ++b) {
    const Real theta = (2.0*b + 1.0) * Geometry::DtoR;
    const R3::XYZ u_c(-sin(theta), 0.0, cos(theta));
    for (int conv = 0; conv < 2; ++conv) {
      for (Index cell = 0; cell < axi5a.NumCells(); ++cell) {
        const Real va = axi5a.AcceptanceBound(RAY_P, u_c, conv, cell);
        const Real vb = axi5b.AcceptanceBound(RAY_P, u_c, conv, cell);
        CHECK_MSG(va == vb,
                 "bin=" << b << " conv=" << conv << " cell=" << cell
                 << ": first-build=" << va << " second-build=" << vb);
      }
    }
  }

}


TEST(axi_geometry_cache_survives_build_exception) {

  // Ledger M3: an exception during BuildAxiGeometry() (aniso_axi.cpp)
  // must not leave a partial entry in GetOrBuildAxiGeometry()'s cache
  // under a fingerprint that a later, correct build at the same
  // take-off-angle count would then wrongly reuse. Neither of
  // BuildAxiGeometry()'s own throws is reachable through AxiProposal's
  // public constructor with real take-off-angle data (the non-
  // icosahedral guard in AxiProposal's own constructor runs *before*
  // the geometry cache is ever touched, by design, so it can't
  // demonstrate cache poisoning either way; the mirror-cell "no match"
  // throw depends only on a canonical, always-exactly-symmetric
  // TesselSphere built from 'coarseDeg' alone, never from the caller's
  // actual directions, so no take-off-angle set that survives the
  // non-icosahedral guard can trigger it) -- see
  // AxiProposal_TestForceNextGeometryBuildToThrow()'s doc comment
  // (aniso_axi.hpp) for why this test uses that one-shot hook instead.
  //
  // Take-off-angle count: a fresh TESS_ICO degree-3 set (1,280
  // directions). No other test in this file builds an AxiProposal at
  // this size (UseTOADegree() is never called with 3 elsewhere; the one
  // other degree-3 TesselSphere in this file, in
  // axi_kernel_matches_reference, feeds PrepareBoundRep()/BoundRatio()
  // directly, never an AxiProposal), so GetOrBuildAxiGeometry() has
  // never built this size before -- the forced throw below is a
  // genuine first-time cache miss, not a fingerprint mismatch against
  // an existing entry.
  S2::S2Set toa;
  {
    const S2::TesselSphere freshToa(S2::TESS_ICO, 3);
    for (Index i = 0; i < (Index) freshToa.size(); ++i) {
      toa.push_back(S2::ThetaPhi(freshToa[i].Theta(), freshToa[i].Phi()));
    }
  }

  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  std::vector<ProbDist> pdists;
  BuildIsoPDists(sp, toa, pdists);

  AxiProposal_TestForceNextGeometryBuildToThrow();
  bool threw = false;
  try {
    AxiProposal axiFail(sp, toa, pdists);
    (void) axiFail;
  } catch (const std::logic_error &) {
    threw = true;
  }
  CHECK_MSG(threw, "forced geometry build failure should propagate out "
           "of AxiProposal's constructor");

  // The real check: a second, unforced build at the exact same
  // take-off-angle set must succeed and produce a fully working
  // AxiProposal. If the cache had been poisoned by the failed build
  // above (the pre-fix bug), this would either throw (MxPerm/MyPerm
  // sized for a different build, or empty) or silently read back
  // zeroed/undersized data.
  bool threwSecond = false;
  AxiProposal * axiOk = 0;
  try {
    axiOk = new AxiProposal(sp, toa, pdists);
  } catch (...) {
    threwSecond = true;
  }
  CHECK_MSG(!threwSecond, "a correct build at the same take-off-angle "
           "count, after a forced failure, must not throw");

  if (axiOk != 0) {
    CHECK(axiOk->NumCells() > 0);
    // Exercise a bound query at every cell -- a partial or zeroed cache
    // entry is the likeliest way a surviving regression would still
    // misbehave silently rather than throw or size-mismatch outright.
    const R3::XYZ u_c(0.0, 0.0, 1.0);
    for (Index cell = 0; cell < axiOk->NumCells(); ++cell) {
      const Real v = axiOk->AcceptanceBound(RAY_P, u_c, 0, cell);
      CHECK_MSG(v >= 0.0, "cell " << cell << ": negative bound " << v);
    }
    delete axiOk;
  }

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 3 TESTS: Parallel construction, bit-identical tables   ****
// ****                                                              ****
//

// ::: struct ThrowingSpawn ::: (ledger I1 test support)
//
//   A ParallelFor() thread-spawn strategy -- parallel_for_detail::
//   DefaultSpawn's counterpart (parallel.hpp) -- that throws
//   std::system_error on its 'failAtCall'-th invocation (0-based)
//   instead of constructing a std::thread, so
//   test_parallel_for_survives_spawn_failure (below) can force,
//   deterministically, the std::thread-construction failure ledger I1
//   hardens ParallelFor against, without needing to actually exhaust
//   the host's real OS thread limit (the only way that failure arises
//   in production -- see ParallelFor's class-level doc comment,
//   parallel.hpp). Passed to the test-only ParallelForSpawn() entry
//   point, never to the public, four-argument ParallelFor().
struct ThrowingSpawn {
  unsigned failAtCall;     // 0-based spawn attempt to fail.
  mutable unsigned calls;  // Spawning runs single-threaded (the calling
                           // thread's own loop in ParallelForSpawn), so
                           // this needs no synchronization.
  explicit ThrowingSpawn(unsigned k) : failAtCall(k), calls(0) {}
  template <class Worker>
  std::thread operator()(Worker & w) const {
    if (calls++ == failAtCall) {
      throw std::system_error(
        std::make_error_code(std::errc::resource_unavailable_try_again),
        "ParallelFor test hook: forced thread spawn failure");
    }
    return std::thread(w);
  }
};

TEST(parallel_for_covers_and_rethrows) {

  // Non-throwing case: for every thread count, a body that increments a
  // dedicated atomic counter per index must visit every index in
  // [0,1000) exactly once -- checked directly (not just inferred from a
  // total), so a double-dispatch of some index racing with a skipped
  // other index (which could otherwise cancel out in a simple sum)
  // can't hide.
  const unsigned threadCounts[] = {0, 1, 3, 64};
  for (std::size_t ti = 0; ti < sizeof(threadCounts)/sizeof(threadCounts[0]); ++ti) {
    const unsigned threads = threadCounts[ti];

    std::vector<std::atomic<int> > hits(1000);
    for (std::size_t i = 0; i < hits.size(); ++i) {hits[i].store(0);}

    ParallelFor((Index) 0, (Index) 1000, threads, [&](Index i) {
      hits[i].fetch_add(1);
    });

    for (Index i = 0; i < 1000; ++i) {
      CHECK_MSG(hits[i].load() == 1,
               "threads=" << threads << " index=" << i
               << " hit-count=" << hits[i].load());
    }
  }

  // Throwing case (threads 3 and 64, per the brief): a body that throws
  // std::logic_error("boom") at i==500 must make ParallelFor rethrow
  // exactly that exception, and no index outside [0,1000) is ever
  // visited.
  const unsigned throwThreadCounts[] = {3, 64};
  for (std::size_t ti = 0; ti < sizeof(throwThreadCounts)/sizeof(throwThreadCounts[0]); ++ti) {
    const unsigned threads = throwThreadCounts[ti];

    std::atomic<Index> maxSeen(0);
    bool threw = false;
    std::string what;
    try {
      ParallelFor((Index) 0, (Index) 1000, threads, [&](Index i) {
        Index prev = maxSeen.load();
        while (i > prev && !maxSeen.compare_exchange_weak(prev, i)) {}
        if (i == 500) {throw std::logic_error("boom");}
      });
    } catch (const std::logic_error & e) {
      threw = true;
      what = e.what();
    } catch (...) {
      threw = true;
      what = "<unexpected exception type>";
    }

    CHECK_MSG(threw, "threads=" << threads << ": expected ParallelFor to throw");
    CHECK_MSG(what == "boom",
             "threads=" << threads << ": what()=\"" << what << "\"");
    CHECK_MSG(maxSeen.load() <= 999,
             "threads=" << threads << ": visited index " << maxSeen.load()
             << " > 999");
  }

  // Broader (n, threads) coverage (brief's Step 1 extension): n values
  // spanning smaller-than-, comparable-to-, and much-larger-than- the
  // chunked claiming grain ('max(1, n/(threads*8))', parallel.hpp), and
  // not always a multiple of it (100003 is prime). For each combination:
  // (1) a non-throwing run visits every index in [0,n) exactly once
  // (same per-index check as above, not just a total); (2) a run that
  // throws std::logic_error("boom") at i==n/2 is rethrown as exactly
  // that exception, and no index is ever visited more than once --
  // regardless of how many higher indices a chunk straddling n/2 skips,
  // or how many chunks past it other threads had already claimed before
  // the failure was recorded.
  const Index extN[] = {1, 7, 1000, 100003};
  const unsigned extThreads[] = {2, 16, 64};
  for (std::size_t ni = 0; ni < sizeof(extN)/sizeof(extN[0]); ++ni) {
    const Index n = extN[ni];
    for (std::size_t ti = 0; ti < sizeof(extThreads)/sizeof(extThreads[0]); ++ti) {
      const unsigned threads = extThreads[ti];

      std::vector<std::atomic<int> > hits(n);
      for (Index i = 0; i < n; ++i) {hits[i].store(0);}
      ParallelFor((Index) 0, n, threads, [&](Index i) {
        hits[i].fetch_add(1);
      });
      for (Index i = 0; i < n; ++i) {
        CHECK_MSG(hits[i].load() == 1,
                 "n=" << n << " threads=" << threads << " index=" << i
                 << " hit-count=" << hits[i].load());
      }

      const Index mid = n / 2;
      std::vector<std::atomic<int> > hits2(n);
      for (Index i = 0; i < n; ++i) {hits2[i].store(0);}
      bool threw2 = false;
      std::string what2;
      try {
        ParallelFor((Index) 0, n, threads, [&](Index i) {
          hits2[i].fetch_add(1);
          if (i == mid) {throw std::logic_error("boom");}
        });
      } catch (const std::logic_error & e) {
        threw2 = true;
        what2 = e.what();
      } catch (...) {
        threw2 = true;
        what2 = "<unexpected exception type>";
      }
      CHECK_MSG(threw2, "n=" << n << " threads=" << threads
               << ": expected ParallelFor to throw");
      CHECK_MSG(what2 == "boom",
               "n=" << n << " threads=" << threads
               << ": what()=\"" << what2 << "\"");
      for (Index i = 0; i < n; ++i) {
        CHECK_MSG(hits2[i].load() <= 1,
                 "n=" << n << " threads=" << threads << " index=" << i
                 << " visited " << hits2[i].load() << " times");
      }
    }
  }

}


// ::: CheckAxiTablesThreadInvariant(sp) :::
//
//   Builds AxiProposal at degree 5 twice -- threads=1 and threads=7 --
//   from the same ScatterParams/take-off-angle set/pdists, and requires
//   MajorantRate(P/S), BinRate() at every P and S bin centre, and
//   AcceptanceBound() for every (bin, conv, cell) to match exactly: the
//   thread count used to build the tables must never change a single
//   bit of them (brief Step 1).
//
static void CheckAxiTablesThreadInvariant(const ScatterParams & sp) {

  const S2::TesselSphere & toas = UseTOADegree(5);
  std::vector<ProbDist> pdists;
  BuildIsoPDists(sp, toas, pdists);

  AxiProposal axi1(sp, toas, pdists, 1);
  AxiProposal axi7(sp, toas, pdists, 7);

  CHECK(axi1.MajorantRate(RAY_P) == axi7.MajorantRate(RAY_P));
  CHECK(axi1.MajorantRate(RAY_S) == axi7.MajorantRate(RAY_S));
  CHECK(axi1.NumCells() == axi7.NumCells());

  for (int b = 0; b < AxiProposal::kThetaBins; ++b) {
    const Real theta = (2.0*b + 1.0) * Geometry::DtoR;  // bin centre
    const R3::XYZ u_c(-sin(theta), 0.0, cos(theta));

    const Real rate1 = axi1.BinRate(RAY_P, u_c);
    const Real rate7 = axi7.BinRate(RAY_P, u_c);
    CHECK_MSG(rate1 == rate7,
             "P BinRate bin=" << b << ": threads=1 " << rate1
             << " threads=7 " << rate7);

    for (int conv = 0; conv < 2; ++conv) {  // GPP, GPS
      for (Index cell = 0; cell < axi1.NumCells(); ++cell) {
        const Real va = axi1.AcceptanceBound(RAY_P, u_c, conv, cell);
        const Real vb = axi7.AcceptanceBound(RAY_P, u_c, conv, cell);
        CHECK_MSG(va == vb,
                 "P bin=" << b << " conv=" << conv << " cell=" << cell
                 << ": threads=1 " << va << " threads=7 " << vb);
      }
    }
  }

  for (int b = 0; b < AxiProposal::kThetaBins; ++b) {
    for (int c = 0; c < AxiProposal::kChiBins; ++c) {
      const Real theta = (2.0*b + 1.0)  * Geometry::DtoR;  // bin centre
      const Real chi   = (10.0*c + 5.0) * Geometry::DtoR;  // bin centre
      const R3::XYZ u_c(-sin(theta)*cos(chi), sin(theta)*sin(chi), cos(theta));

      const Real rate1 = axi1.BinRate(RAY_S, u_c);
      const Real rate7 = axi7.BinRate(RAY_S, u_c);
      CHECK_MSG(rate1 == rate7,
               "S BinRate bin=" << b << "," << c << ": threads=1 " << rate1
               << " threads=7 " << rate7);

      for (int conv = 2; conv < 4; ++conv) {  // GSP, GSS
        for (Index cell = 0; cell < axi1.NumCells(); ++cell) {
          const Real va = axi1.AcceptanceBound(RAY_S, u_c, conv, cell);
          const Real vb = axi7.AcceptanceBound(RAY_S, u_c, conv, cell);
          CHECK_MSG(va == vb,
                   "S bin=" << b << "," << c << " conv=" << conv
                   << " cell=" << cell << ": threads=1 " << va
                   << " threads=7 " << vb);
        }
      }
    }
  }

}


TEST(parallel_for_survives_spawn_failure) {

  // Ledger I1: if std::thread construction throws partway through
  // ParallelFor's spawn loop, ParallelFor must not let that reach the
  // caller while it still holds already-started, un-joined threads --
  // std::thread's destructor calls std::terminate() in that case,
  // turning a resource-limit error (the production cause: an
  // unreasonably large --threads value exceeding the OS thread-creation
  // limit, since main.cpp accepts --threads up to UINT_MAX with no
  // upper-bound check) into an abort with no diagnostic. This test
  // forces that failure deterministically via the test-only
  // ParallelForSpawn() entry point and ThrowingSpawn (above), instead of
  // actually exhausting the host's real OS thread limit. If
  // ParallelFor's fix regressed, this test would not report a CHECK
  // failure -- std::terminate() would abort the whole test binary before
  // reaching any CHECK_MSG below -- so merely completing this test is
  // itself part of what it checks.
  //
  // Per-index coverage, every index visited exactly once, at several
  // (threads, fail-at-call) combinations -- including failing on the
  // very first spawn attempt (failAt=0: zero real threads ever start,
  // the calling thread alone must drain everything) and failing after
  // a few threads have already started (failAt=2: those threads keep
  // running concurrently with the calling thread's own drain):
  const unsigned threadCounts[] = {4, 16, 64};
  const unsigned failAtCalls[]  = {0, 2};
  for (std::size_t ti = 0; ti < sizeof(threadCounts)/sizeof(threadCounts[0]); ++ti) {
    for (std::size_t fi = 0; fi < sizeof(failAtCalls)/sizeof(failAtCalls[0]); ++fi) {
      const unsigned threads = threadCounts[ti];
      const unsigned failAt  = failAtCalls[fi];
      if (failAt >= threads) {continue;}  // Nothing to skip past.

      const Index n = 1000;
      std::vector<std::atomic<int> > hits(n);
      for (Index i = 0; i < n; ++i) {hits[i].store(0);}

      ParallelForSpawn((Index) 0, n, threads, [&](Index i) {
        hits[i].fetch_add(1);
      }, ThrowingSpawn(failAt));

      for (Index i = 0; i < n; ++i) {
        CHECK_MSG(hits[i].load() == 1,
                 "threads=" << threads << " failAt=" << failAt
                 << " index=" << i << " hit-count=" << hits[i].load());
      }
    }
  }

  // The forced spawn failure must not itself surface as an exception to
  // the caller (parallel.hpp's class-level doc comment: the call
  // completes "with fewer effective workers", not an error) -- distinct
  // from a genuine body(i) exception, which must still propagate exactly
  // as parallel_for_covers_and_rethrows (above) already checks with the
  // default spawn strategy; ParallelForSpawn() shares that exact same
  // body/exception-handling code, parametrized only on how threads get
  // spawned, so this exercises the one combination that test couldn't
  // reach (a body(i) exception *and* a spawn failure in the same call).
  {
    bool threw = false;
    try {
      ParallelForSpawn((Index) 0, (Index) 1000, 8u,
                       [&](Index) {}, ThrowingSpawn(1));
    } catch (...) {
      threw = true;
    }
    CHECK_MSG(!threw, "a forced thread-spawn failure must not itself "
             "propagate out of ParallelFor");
  }
  {
    bool threw = false;
    std::string what;
    try {
      ParallelForSpawn((Index) 0, (Index) 1000, 8u, [&](Index i) {
        if (i == 500) {throw std::logic_error("boom");}
      }, ThrowingSpawn(1));
    } catch (const std::logic_error & e) {
      threw = true;
      what = e.what();
    } catch (...) {
      threw = true;
      what = "<unexpected exception type>";
    }
    CHECK_MSG(threw, "a genuine body(i) exception must still propagate "
             "even when a thread-spawn failure also occurred");
    CHECK_MSG(what == "boom", "what()=\"" << what << "\"");
  }

}


TEST(axi_tables_thread_invariant) {
  CheckAxiTablesThreadInvariant(LopNorCrust(2, 2, 0.5));
  CheckAxiTablesThreadInvariant(LopNorCrust(0.5, 0.5, 2));
}


TEST(report_stats_thread_invariant) {

  UseTOADegree(5);  // Degree for each Scatterer's own (general-sampler
                    // and/or axisymmetric) probability tables;
                    // ComputeReportStats() always integrates over the
                    // fixed degree-6 ReportDirs() regardless.

  ScatterParams spAxi = LopNorCrust(2, 2, 0.5);  // axisymmetric: 8-state
                                                 // GetReportMFPs() branch
  ScatterParams spTri = LopNorCrust(2, 1, 0.5);  // triaxial: 9-state
                                                 // GetReportMFPs() branch

  Scatterer::SetBuildThreads(1);
  Scatterer axiScat1(spAxi);
  Scatterer triScat1(spTri);
  const std::vector<Real> mfpsAxi1 = axiScat1.GetReportMFPs();
  const std::vector<Real> mfpsTri1 = triScat1.GetReportMFPs();

  Scatterer::SetBuildThreads(7);
  Scatterer axiScat7(spAxi);
  Scatterer triScat7(spTri);
  const std::vector<Real> mfpsAxi7 = axiScat7.GetReportMFPs();
  const std::vector<Real> mfpsTri7 = triScat7.GetReportMFPs();

  Scatterer::SetBuildThreads(1);  // restore

  CHECK(mfpsAxi1.size() == 8);
  CHECK(mfpsAxi7.size() == 8);
  CHECK(mfpsTri1.size() == 9);
  CHECK(mfpsTri7.size() == 9);

  for (std::size_t i = 0; i < mfpsAxi1.size(); ++i) {
    CHECK_MSG(mfpsAxi1[i] == mfpsAxi7[i],
             "axisymmetric GetReportMFPs()[" << i << "]: threads=1 "
             << mfpsAxi1[i] << " threads=7 " << mfpsAxi7[i]);
  }
  for (std::size_t i = 0; i < mfpsTri1.size(); ++i) {
    CHECK_MSG(mfpsTri1[i] == mfpsTri7[i],
             "triaxial GetReportMFPs()[" << i << "]: threads=1 "
             << mfpsTri1[i] << " threads=7 " << mfpsTri7[i]);
  }

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 4b TESTS: Parallel geometry, bit-identical            ****
// ****                                                              ****
//

// ::: CheckAxiGeometryThreadInvariant(degree) :::
//
//   Proves BuildAxiGeometry() (aniso_axi.cpp) -- specifically Task 4b's
//   two newly-parallelized loops, the representative-radius loop and
//   the mirror-cell match -- is bit-identical for threads=1 vs.
//   threads=7, by forcing two independent, *uncached* builds of the
//   same take-off-angle content.
//
//   Why "uncached" matters: GetOrBuildAxiGeometry()'s cache is keyed
//   only by toa.size(), so simply building twice with the *same*
//   S2::S2Set (as CheckAxiTablesThreadInvariant, above, does) would, on
//   the second call, just hand back the geometry already cached from
//   the first -- whatever thread count that first call happened to use
//   -- regardless of the 'threads' argument the second call passes. That
//   proves the (already-parallel, Task 3) bin tables are thread-
//   invariant, but says nothing about BuildAxiGeometry() itself, which
//   never even reruns on a cache hit.
//
//   The fix (brief Step 1): build once with the *standard* degree-D
//   take-off-angle set, fixing (or merely confirming, if some earlier
//   test already did) that size's cache entry to the standard
//   fingerprint; then build twice more with a same-size but
//   differently-fingerprinted set -- the 0.01 rad-about-z rotation from
//   test_axi_geometry_cache_mismatch_rotation, confirmed there to take
//   the fingerprint-mismatch path observably -- once per thread count.
//   GetOrBuildAxiGeometry()'s mismatch path *always* builds fresh into a
//   scratch object and never writes the cache (by design, so a
//   same-size-but-different set can never evict or overwrite the real
//   entry), so every call with the rotated set is an independent,
//   from-scratch BuildAxiGeometry() call no matter how many times it's
//   been called before -- the threads=1 and threads=7 rotated builds
//   below are two genuinely separate computations of the same geometry,
//   not one computation compared against itself.
//
//   Checks (brief Step 1) MajorantRate(P/S) and AcceptanceBound() at
//   every P bin x conv x cell, both at the P-bin-centre directions
//   themselves (x<=0, so Fold()'s Mx is always false -- this exercises
//   PrepareBoundRep()'s per-representative data, including the
//   representative-radius loop's output, since every cell's bound
//   depends on it) and, additionally, at their x-mirrored counterparts
//   (x>=0, Mx always true for these bin centres), which route every
//   AcceptanceBound() call through MapCell() and so exercise the
//   mirror-match loop's MxPerm output too -- the brief's P-bin-centre
//   check alone cannot do this, since Fold() never sets Mx for a
//   direction already at x<=0.
//
static void CheckAxiGeometryThreadInvariant(int degree) {

  const S2::TesselSphere & toas = UseTOADegree(degree);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);

  // Fixes (or confirms) this take-off-angle size's cache entry to the
  // standard fingerprint, so every build below with the rotated set is
  // guaranteed to take the mismatch (uncached, from-scratch) path.
  std::vector<ProbDist> pdistsStd;
  BuildIsoPDists(sp, toas, pdistsStd);
  AxiProposal axiWarm(sp, toas, pdistsStd, 1);
  (void) axiWarm;

  S2::S2Set rotated;
  const Real phi0 = 0.01;  // Same rotation as
                          // test_axi_geometry_cache_mismatch_rotation;
                          // confirmed there to change the fingerprint
                          // (observable AcceptanceBound difference from
                          // the standard set).
  for (Index i = 0; i < (Index) toas.size(); ++i) {
    rotated.push_back(S2::ThetaPhi(toas[i].Theta(), toas[i].Phi() + phi0));
  }
  CHECK(rotated.size() == toas.size());

  std::vector<ProbDist> pdistsRot1;
  BuildIsoPDists(sp, rotated, pdistsRot1);
  AxiProposal axi1(sp, rotated, pdistsRot1, 1);

  std::vector<ProbDist> pdistsRot7;
  BuildIsoPDists(sp, rotated, pdistsRot7);
  AxiProposal axi7(sp, rotated, pdistsRot7, 7);

  CHECK_MSG(axi1.MajorantRate(RAY_P) == axi7.MajorantRate(RAY_P),
           "degree=" << degree << " MajorantRate(P): threads=1 "
           << axi1.MajorantRate(RAY_P) << " threads=7 "
           << axi7.MajorantRate(RAY_P));
  CHECK_MSG(axi1.MajorantRate(RAY_S) == axi7.MajorantRate(RAY_S),
           "degree=" << degree << " MajorantRate(S): threads=1 "
           << axi1.MajorantRate(RAY_S) << " threads=7 "
           << axi7.MajorantRate(RAY_S));
  CHECK(axi1.NumCells() == axi7.NumCells());

  for (int b = 0; b < AxiProposal::kThetaBins; ++b) {
    const Real theta = (2.0*b + 1.0) * Geometry::DtoR;  // bin centre
    const R3::XYZ u_c      (-sin(theta), 0.0, cos(theta));  // Mx=false
    const R3::XYZ u_cMirror( sin(theta), 0.0, cos(theta));  // Mx=true

    for (int conv = 0; conv < 2; ++conv) {  // GPP, GPS
      for (Index cell = 0; cell < axi1.NumCells(); ++cell) {

        const Real va = axi1.AcceptanceBound(RAY_P, u_c, conv, cell);
        const Real vb = axi7.AcceptanceBound(RAY_P, u_c, conv, cell);
        CHECK_MSG(va == vb,
                 "degree=" << degree << " bin=" << b << " conv=" << conv
                 << " cell=" << cell << ": threads=1 " << va
                 << " threads=7 " << vb);

        const Real vam = axi1.AcceptanceBound(RAY_P, u_cMirror, conv, cell);
        const Real vbm = axi7.AcceptanceBound(RAY_P, u_cMirror, conv, cell);
        CHECK_MSG(vam == vbm,
                 "degree=" << degree << " bin=" << b << " conv=" << conv
                 << " cell=" << cell << " (x-mirrored): threads=1 " << vam
                 << " threads=7 " << vbm);
      }
    }
  }

}

TEST(axi_geometry_thread_invariant) {
  CheckAxiGeometryThreadInvariant(5);  // repEqualsToa: exercises the
                                       // mirror match only (brief Step 1).
  CheckAxiGeometryThreadInvariant(7);  // exercises the representative-
                                       // radius acos loop too (nonzero
                                       // repBlock, same as
                                       // compare_aniso_reference.sh's
                                       // choice of degree 7).
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 4c TESTS: Parallel isotropic table construction,       ****
// ****                 bit-identical, with cached PSATO constants   ****
// ****                                                              ****
//

// ::: PhononObs / DrawPhononSequence(s, intype, seed, count) :::
//
//   Draws 'count' Phonons from s.GetRandomScatteredRelativePhonon(intype)
//   and records each one's (theta, phi, polarization, full raytype) for
//   later comparison. The draws run under a freshly-seeded RandomContext
//   (random.hpp), installed and torn down by a small RAII guard, so:
//     (a) the same 'seed' always reproduces the same draw sequence,
//         regardless of what randomness (if any) ran before this call --
//         letting a caller replay an *identical* sequence of
//         RandomUniform()/RandomInteger() draws against two differently-
//         built Scatterers and attribute any difference in the results
//         to the tables, not to the draws; and
//     (b) the thread_local active-context pointer is always restored to
//         null before returning -- including if a draw throws -- so a
//         failure here can't leave a later test running against a
//         dangling pointer to this function's now-destroyed local
//         RandomContext.
//
struct PhononObs {
  Real theta, phi, pol;
  raytype rt;
};

static std::vector<PhononObs> DrawPhononSequence(Scatterer & s,
                                                  raytype intype,
                                                  std::uint64_t seed,
                                                  int count) {
  struct ScopedRandomContext {
    RandomContext ctx;
    explicit ScopedRandomContext(std::uint64_t sd) : ctx(sd) {
      SetActiveRandomContext(&ctx);
    }
    ~ScopedRandomContext() {SetActiveRandomContext(0);}
  } guard(seed);

  std::vector<PhononObs> obs;
  obs.reserve(count);
  for (int i = 0; i < count; ++i) {
    const Phonon p = s.GetRandomScatteredRelativePhonon(intype);
    PhononObs o;
    o.theta = p.GetDirection().Theta();
    o.phi   = p.GetDirection().Phi();
    o.pol   = p.GetPolarization();
    o.rt    = p.GetFullRaytype();
    obs.push_back(o);
  }
  return obs;
}

// TEST(iso_tables_thread_invariant):
//
//   Guard test (brief Step 1): construction-time parallelism must never
//   change a table, a report value, or any output of a subsequently-
//   constructed Scatterer, whatever Scatterer::GetBuildThreads() was set
//   to during its build. Builds an isotropic and an axisymmetric
//   Scatterer at threads=1, then again at threads=7, and requires
//   GetMajorantMeanFreePath(P/S) (both specs), GetReportMFPs() (the
//   axisymmetric spec only -- empty, and so trivially equal, for an
//   isotropic one), and the first 2000 results of a fixed-seed
//   RandomUniform()/RandomInteger() draw sequence through
//   GetRandomScatteredRelativePhonon(RAY_P) and (RAY_S) on the isotropic
//   Scatterer (the only sampling path that method supports; see its
//   guard against anisotropic ScatterParams) to come out exactly equal.
//
//   This exercises both halves of Task 4c: a wrong ScatterParams::
//   PSATO() numerator cache would perturb every GSATO() value feeding
//   PopulateProbDists() (and, transitively, ComputeReportStats() and the
//   axisymmetric sampler's AxiProposal, built from the same mPDists); a
//   data race or ordering bug in PopulateProbDists()'s new ParallelFor
//   would perturb mPDists/m_spol directly. Either would change which
//   take-off-angle index, conversion, and/or polarization get drawn.
//
//   Like task-4b-report.md's analogous test, this is RED-free: it can
//   only fail if the new parallel implementation itself has a bug, not
//   because today's (pre-Task-4c) construction already varied with
//   thread count -- recorded as such per the brief's own instruction for
//   this situation.
//
TEST(iso_tables_thread_invariant) {

  UseTOADegree(5);

  Scatterer::SetBuildThreads(1);
  Scatterer isoScat1(LopNorCrust(0.5, 0.5, 0.5));   // isotropic
  Scatterer axiScat1(LopNorCrust(2,   2,   0.5));   // axisymmetric

  const Real mfpIsoP1 = isoScat1.GetMajorantMeanFreePath(RAY_P);
  const Real mfpIsoS1 = isoScat1.GetMajorantMeanFreePath(RAY_S);
  const Real mfpAxiP1 = axiScat1.GetMajorantMeanFreePath(RAY_P);
  const Real mfpAxiS1 = axiScat1.GetMajorantMeanFreePath(RAY_S);
  const std::vector<Real> reportMfps1 = axiScat1.GetReportMFPs();

  const std::uint64_t seedP = 1357111317ULL;
  const std::uint64_t seedS = 2468121416ULL;
  const int kDraws = 2000;

  const std::vector<PhononObs> drawsP1 =
    DrawPhononSequence(isoScat1, RAY_P, seedP, kDraws);
  const std::vector<PhononObs> drawsS1 =
    DrawPhononSequence(isoScat1, RAY_S, seedS, kDraws);

  Scatterer::SetBuildThreads(7);
  Scatterer isoScat7(LopNorCrust(0.5, 0.5, 0.5));
  Scatterer axiScat7(LopNorCrust(2,   2,   0.5));

  const Real mfpIsoP7 = isoScat7.GetMajorantMeanFreePath(RAY_P);
  const Real mfpIsoS7 = isoScat7.GetMajorantMeanFreePath(RAY_S);
  const Real mfpAxiP7 = axiScat7.GetMajorantMeanFreePath(RAY_P);
  const Real mfpAxiS7 = axiScat7.GetMajorantMeanFreePath(RAY_S);
  const std::vector<Real> reportMfps7 = axiScat7.GetReportMFPs();

  const std::vector<PhononObs> drawsP7 =
    DrawPhononSequence(isoScat7, RAY_P, seedP, kDraws);
  const std::vector<PhononObs> drawsS7 =
    DrawPhononSequence(isoScat7, RAY_S, seedS, kDraws);

  Scatterer::SetBuildThreads(1);  // restore

  CHECK_MSG(mfpIsoP1 == mfpIsoP7, "isotropic MFP(P): threads=1 " << mfpIsoP1
           << " threads=7 " << mfpIsoP7);
  CHECK_MSG(mfpIsoS1 == mfpIsoS7, "isotropic MFP(S): threads=1 " << mfpIsoS1
           << " threads=7 " << mfpIsoS7);
  CHECK_MSG(mfpAxiP1 == mfpAxiP7, "axisymmetric MFP(P): threads=1 "
           << mfpAxiP1 << " threads=7 " << mfpAxiP7);
  CHECK_MSG(mfpAxiS1 == mfpAxiS7, "axisymmetric MFP(S): threads=1 "
           << mfpAxiS1 << " threads=7 " << mfpAxiS7);

  CHECK(reportMfps1.size() == reportMfps7.size());
  for (std::size_t i = 0; i < reportMfps1.size(); ++i) {
    CHECK_MSG(reportMfps1[i] == reportMfps7[i],
             "GetReportMFPs()[" << i << "]: threads=1 " << reportMfps1[i]
             << " threads=7 " << reportMfps7[i]);
  }

  CHECK(drawsP1.size() == drawsP7.size());
  for (std::size_t i = 0; i < drawsP1.size(); ++i) {
    CHECK_MSG(drawsP1[i].theta == drawsP7[i].theta &&
             drawsP1[i].phi   == drawsP7[i].phi   &&
             drawsP1[i].pol   == drawsP7[i].pol   &&
             drawsP1[i].rt    == drawsP7[i].rt,
             "RAY_P draw " << i << " differs: threads=1 (theta="
             << drawsP1[i].theta << " phi=" << drawsP1[i].phi << " pol="
             << drawsP1[i].pol << " rt=" << drawsP1[i].rt
             << ") threads=7 (theta=" << drawsP7[i].theta << " phi="
             << drawsP7[i].phi << " pol=" << drawsP7[i].pol << " rt="
             << drawsP7[i].rt << ")");
  }

  CHECK(drawsS1.size() == drawsS7.size());
  for (std::size_t i = 0; i < drawsS1.size(); ++i) {
    CHECK_MSG(drawsS1[i].theta == drawsS7[i].theta &&
             drawsS1[i].phi   == drawsS7[i].phi   &&
             drawsS1[i].pol   == drawsS7[i].pol   &&
             drawsS1[i].rt    == drawsS7[i].rt,
             "RAY_S draw " << i << " differs: threads=1 (theta="
             << drawsS1[i].theta << " phi=" << drawsS1[i].phi << " pol="
             << drawsS1[i].pol << " rt=" << drawsS1[i].rt
             << ") threads=7 (theta=" << drawsS7[i].theta << " phi="
             << drawsS7[i].phi << " pol=" << drawsS7[i].pol << " rt="
             << drawsS7[i].rt << ")");
  }

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 5 (PERFORMANCE PLAN) TESTS: Faster axisymmetric        ****
// ****                 sampling path                                ****
// ****                                                              ****
//
// See docs/superpowers/specs/2026-10-17-axisymmetric-sampler-
// performance-design.md, decision 4, and .superpowers/sdd/
// 2026-10-17-axisymmetric-sampler-performance/task-5-brief.md.
//

// ::: ReferenceLocate3db253a(intype, u) :::
//
//   Ruling R-P1's frozen oracle for AxiProposal::Locate(): a free-
//   standing, literal transcription of AxiProposal::Fold()'s body
//   exactly as it existed before this task -- confirmed identical to
//   3db253a's Fold() (the commit this branch, and the performance plan,
//   starts from) via `git show 3db253a:aniso_axi.cpp` before this task
//   began, the same freezing rationale as ReferenceBound3db253a above.
//   Fold() was private and returned (private) FoldResult; this
//   free-standing copy takes 'intype' and 'u' as plain parameters and
//   returns the now-public AxiProposal::BinRef instead, so it can be
//   compared directly against Locate()'s return value -- nothing else
//   is changed. Do not "fix" or simplify this function if it starts
//   looking stale -- it must stay a literal transcription of the frozen
//   commit, bugs and all.
//
static AxiProposal::BinRef ReferenceLocate3db253a(raytype intype, R3::XYZ u) {

  if (u.z() < 0.0) {u = u.Negative();}

  const bool mx = (u.x() > 0.0);
  const Real ux = mx ? -u.x() : u.x();

  bool my = false;
  Real uy = 0.0;
  if (intype == RAY_S) {
    my = (u.y() < 0.0);
    uy = my ? -u.y() : u.y();
  }

  Real uz = u.z();
  if (uz > 1.0)  uz = 1.0;
  if (uz < -1.0) uz = -1.0;
  const Real theta = acos(uz);
  const Real sinTheta = sin(theta);

  Real chi = 0.0;
  if (intype == RAY_S && sinTheta > 1e-12) {
    chi = atan2(uy, -ux);
  }

  int thetaBin = (int) floor(theta*Geometry::RtoD/2.0);
  if (thetaBin < 0) thetaBin = 0;
  if (thetaBin > AxiProposal::kThetaBins-1) thetaBin = AxiProposal::kThetaBins-1;

  int chiBin = 0;
  if (intype == RAY_S) {
    chiBin = (int) floor(chi*Geometry::RtoD/10.0);
    if (chiBin < 0) chiBin = 0;
    if (chiBin > AxiProposal::kChiBins-1) chiBin = AxiProposal::kChiBins-1;
  }

  AxiProposal::BinRef r;
  r.ThetaBin = thetaBin;
  r.ChiBin   = chiBin;
  r.Mx       = mx;
  r.My       = my;
  return r;

}


TEST(axi_locate_matches_reference) {

  // A minimal AxiProposal just to call Locate() on -- its geometry and
  // tables don't affect Locate(), which depends only on kThetaBins/
  // kChiBins (fixed constants) and 'u_rel'/'intype'.
  const S2::TesselSphere & toas = UseTOADegree(2);
  ScatterParams sp = LopNorCrust(2, 2, 0.5);
  std::vector<ProbDist> pdists;
  BuildIsoPDists(sp, toas, pdists);
  AxiProposal axi(sp, toas, pdists);

  std::vector<R3::XYZ> dirs;

  // 200,000 Fibonacci-sphere directions.
  const int N = 200000;
  const Real goldenAngle = Geometry::Pi * (3.0 - sqrt(5.0));
  for (int i = 0; i < N; ++i) {
    const Real zz  = 1.0 - 2.0*((Real) i + 0.5)/(Real) N;
    const Real rr  = sqrt(std::max((Real)0.0, 1.0 - zz*zz));
    const Real phi = goldenAngle * (Real) i;
    dirs.push_back(R3::XYZ(rr*cos(phi), rr*sin(phi), zz));
  }

  // Every theta edge (2b deg), at chi==45deg, +-1e-9 rad:
  for (int b = 0; b <= AxiProposal::kThetaBins; ++b) {
    const Real chi0 = 45.0*Geometry::DtoR;
    for (int s = -1; s <= 1; s += 2) {
      const Real theta = 2.0*b*Geometry::DtoR + s*1e-9;
      dirs.push_back(R3::XYZ(-sin(theta)*cos(chi0),
                             sin(theta)*sin(chi0),
                             cos(theta)));
    }
  }

  // Every chi edge (10c deg), at theta==45deg, +-1e-9 rad:
  for (int c = 0; c <= AxiProposal::kChiBins; ++c) {
    const Real theta0 = 45.0*Geometry::DtoR;
    for (int s = -1; s <= 1; s += 2) {
      const Real chi = 10.0*c*Geometry::DtoR + s*1e-9;
      dirs.push_back(R3::XYZ(-sin(theta0)*cos(chi),
                             sin(theta0)*sin(chi),
                             cos(theta0)));
    }
  }

  // Axis-aligned directions.
  dirs.push_back(R3::XYZ(0,0,1));  dirs.push_back(R3::XYZ(0,0,-1));
  dirs.push_back(R3::XYZ(1,0,0));  dirs.push_back(R3::XYZ(-1,0,0));
  dirs.push_back(R3::XYZ(0,1,0));  dirs.push_back(R3::XYZ(0,-1,0));

  const raytype intypes[2] = {RAY_P, RAY_S};
  for (int ii = 0; ii < 2; ++ii) {
    const raytype intype = intypes[ii];
    for (std::size_t i = 0; i < dirs.size(); ++i) {
      const R3::XYZ & u = dirs[i];
      const AxiProposal::BinRef got  = axi.Locate(intype, u);
      const AxiProposal::BinRef want = ReferenceLocate3db253a(intype, u);

      CHECK_MSG(got.Mx == want.Mx,
               "Mx mismatch: intype=" << ii << " i=" << i);
      CHECK_MSG(got.My == want.My,
               "My mismatch: intype=" << ii << " i=" << i);

      if (got.ThetaBin != want.ThetaBin) {
        // Reconstruct theta via acos (fine in a test) to confirm we're
        // within the 1e-12 rad of a bin edge the brief allows either
        // neighbor at.
        R3::XYZ v = u; if (v.z() < 0.0) v = v.Negative();
        Real vz = v.z(); if (vz > 1.0) vz = 1.0; if (vz < -1.0) vz = -1.0;
        const Real theta = acos(vz);
        const Real edgeStep = 2.0*Geometry::DtoR;
        const Real nearestEdge = edgeStep * floor(theta/edgeStep + 0.5);
        const int dTheta = got.ThetaBin - want.ThetaBin;
        CHECK_MSG((dTheta == 1 || dTheta == -1) &&
                 std::fabs(theta-nearestEdge) <= 1e-12,
                 "ThetaBin mismatch far from an edge: intype=" << ii
                 << " i=" << i << " got=" << got.ThetaBin
                 << " want=" << want.ThetaBin << " theta=" << theta);
      }

      if (intype == RAY_S && got.ChiBin != want.ChiBin) {
        R3::XYZ v = u; if (v.z() < 0.0) v = v.Negative();
        const bool vmx = (v.x() > 0.0);
        const Real vux = vmx ? -v.x() : v.x();
        const bool vmy = (v.y() < 0.0);
        const Real vuy = vmy ? -v.y() : v.y();
        Real vz = v.z(); if (vz > 1.0) vz = 1.0; if (vz < -1.0) vz = -1.0;
        const Real theta = acos(vz);
        const Real sinTheta = sin(theta);

        if (sinTheta <= 1e-12) {
          // Within about 1.5e-8 rad of a pole, cos(theta) itself rounds
          // to exactly +-1 in double precision (1-cos(theta) underflows
          // relative to 1's ULP), so ReferenceLocate3db253a's acos/sin
          // round-trip always recovers sinTheta==0 and so always
          // defaults to chi==0 (its own 'sinTheta > 1e-12' guard) --
          // regardless of the true azimuth. That default is an artifact
          // of the lossy trig round-trip, not a meaningful fold: this
          // direction is at/within the chi guard's own degeneracy
          // radius of a pole either way (both agree ThetaBin==0, above).
          // Locate() instead reads the azimuth from u.x()/u.y() directly
          // (no round-trip), so it can resolve a real ChiBin here
          // without contradicting the reference's own notion of when
          // chi is undefined. Any got.ChiBin is accepted in this case.
        } else {
          const Real chi = atan2(vuy, -vux);
          const Real edgeStep = 10.0*Geometry::DtoR;
          const Real nearestEdge = edgeStep * floor(chi/edgeStep + 0.5);
          const int dChi = got.ChiBin - want.ChiBin;
          CHECK_MSG((dChi == 1 || dChi == -1) &&
                   std::fabs(chi-nearestEdge) <= 1e-12,
                   "ChiBin mismatch far from an edge: intype=" << ii
                   << " i=" << i << " got=" << got.ChiBin
                   << " want=" << want.ChiBin << " chi=" << chi);
        }
      }
    }
  }

}


TEST(axi_ratio_matches_fabric) {

  ScatterParams sp = LopNorCrust(2, 2, 0.5);  // axisymmetric

  struct MapCase {bool curved; R3::XYZ loc;};
  const MapCase cases[2] = {
    {false, R3::XYZ(0,0,0)},
    {true,  R3::XYZ(950.0,0.0,0.0)},
  };

  for (int ci = 0; ci < 2; ++ci) {
    if (cases[ci].curved) {ECS.SetMapping(EarthCoords::RAE_CURVED);}
    const R3::XYZ loc = cases[ci].loc;
    const R3::XYZ up = ECS.GetUp(loc);

    std::mt19937 rng(4040u + (unsigned) ci);
    std::uniform_real_distribution<Real> uTheta(0.0, Geometry::Pi);
    std::uniform_real_distribution<Real> uPhi(0.0, 2.0*Geometry::Pi);
    std::uniform_real_distribution<Real> uPol(0.0, 2.0*Geometry::Pi);
    std::uniform_int_distribution<int> uConv(0, 3);

    const Real kP = sp.GetKP(), kS = sp.GetKS();

    for (int i = 0; i < 1000; ++i) {
      const S2::ThetaPhi dir(uTheta(rng), uPhi(rng));
      const Real pol = uPol(rng);
      const S2::ThetaPhi t(uTheta(rng), uPhi(rng));
      const int conv = uConv(rng);

      const Real k_in  = (conv==0||conv==1) ? kP : kS;
      const Real k_out = (conv==0||conv==2) ? kP : kS;

      const R3::OrthoAxes AA(dir.Theta(), dir.Phi(), pol);
      const R3::XYZ m = AA.Express(t.XYZ()).ScaledBy(k_out)
                      + dir.XYZ().ScaledBy(-k_in);

      const Real mu = m.Dot(up);
      const Real lhs = sp.AxiRatio(m.MagSquared(), mu*mu);
      const Real rhs = sp.AnisoRatio(Scatterer::ToFabric(loc, m));

      CHECK_MSG(NearRel(lhs, rhs, 1e-12),
               "mapping=" << ci << " i=" << i << " conv=" << conv
               << ": lhs=" << lhs << " rhs=" << rhs);
    }
  }

  ECS.SetMapping(EarthCoords::ENU_ORTHO);  // restore for later tests

}


TEST(axi_p_urel_matches_canonical) {

  struct MapCase {bool curved; R3::XYZ loc;};
  const MapCase cases[2] = {
    {false, R3::XYZ(0,0,0)},
    {true,  R3::XYZ(950.0,0.0,0.0)},
  };

  for (int ci = 0; ci < 2; ++ci) {
    if (cases[ci].curved) {ECS.SetMapping(EarthCoords::RAE_CURVED);}
    const R3::XYZ loc = cases[ci].loc;
    const R3::XYZ up = ECS.GetUp(loc);

    std::vector<S2::ThetaPhi> dirs;
    std::mt19937 rng(5050u + (unsigned) ci);
    std::uniform_real_distribution<Real> uz(-1.0, 1.0);
    std::uniform_real_distribution<Real> uphi(0.0, 2.0*Geometry::Pi);
    for (int i = 0; i < 1000; ++i) {
      const Real z = uz(rng);
      dirs.push_back(S2::ThetaPhi(acos(z), uphi(rng)));
    }
    dirs.push_back(S2::ThetaPhi(up.Theta(), up.Phi()));           // dir == +up
    const R3::XYZ negUp = up.ScaledBy(-1.0);
    dirs.push_back(S2::ThetaPhi(negUp.Theta(), negUp.Phi()));     // dir == -up

    for (std::size_t i = 0; i < dirs.size(); ++i) {
      const S2::ThetaPhi & dir = dirs[i];
      Real cz = up.Dot(dir.XYZ());
      if (cz > 1.0)  cz = 1.0;
      if (cz < -1.0) cz = -1.0;
      const R3::XYZ expected(-sqrt(std::max((Real)0.0, 1.0-cz*cz)), 0.0, cz);

      const Real rot = Scatterer::CanonicalPPolarization(dir, up);
      const R3::OrthoAxes AA(dir.Theta(), dir.Phi(), rot);
      const R3::XYZ actual(up.Dot(AA.S1()), up.Dot(AA.S2()), up.Dot(AA.S3()));

      CHECK_MSG(std::fabs(actual.x()-expected.x()) < 1e-12 &&
               std::fabs(actual.y()-expected.y()) < 1e-12 &&
               std::fabs(actual.z()-expected.z()) < 1e-12,
               "mapping=" << ci << " i=" << i << ": expected=("
               << expected.x() << "," << expected.y() << "," << expected.z()
               << ") actual=(" << actual.x() << "," << actual.y() << ","
               << actual.z() << ")");
    }
  }

  ECS.SetMapping(EarthCoords::ENU_ORTHO);  // restore for later tests

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 5b TESTS: Vectorized bound kernel, bit-identical        ****
// ****                                                              ****
//

// ::: CheckAxiRatioRunMatchesScalar<AhGeAz>(sp) :::
//
//   Brief's Step 1 test: builds the representative geometry a degree-7
//   AxiProposal construction would use internally (AxiProposal::
//   kRepDegree==6 caps the representative count below degree 7's own
//   327,680 fine directions, so repBlock==4 and every representative's
//   radius is nonzero -- the same reason
//   CheckAxiGeometryThreadInvariant(7), above, picks degree 7).
//   GetOrBuildAxiGeometry()/BuildAxiGeometry() (aniso_axi.cpp) that would
//   normally compute this are file-local (anonymous namespace) and
//   keyed to a process-wide cache this test shouldn't disturb, so the
//   geometry is reproduced here directly from AxiProposal's own public
//   degree constants instead -- a fresh degree-kRepDegree TesselSphere
//   for the centres, each one's radius the largest angle to any of its
//   own block of degree-7 fine directions, exactly BuildAxiGeometry()'s
//   formula.
//
//   For every one of AxiProposal's own 45 P and 405 S incidence bins
//   (same u_c/delta formulas as the constructor and
//   CheckAxiTablesThreadInvariant, above) and both of that bin type's
//   (k_in,k_out) pairs (P: (kP,kP),(kP,kS); S: (kS,kP),(kS,kS)),
//   RatioRun<AhGeAz>()'s output must equal BoundRatio()'s, exactly, for
//   every representative -- not just within a tolerance, since the only
//   difference between them is meant to be scalar vs. (GCC-)vectorized
//   code generation for the *same* arithmetic (see aniso_axi.cpp's
//   RatioRun() doc comment): any mismatch would mean the vectorized loop
//   rounds some representative/bin combination differently, which the
//   task's "never accept a mismatch" rule treats as a bug to fix, not a
//   tolerance to loosen.
//
//   Each RatioRun<AhGeAz>() call below covers the full [0,numReps)
//   range in one call, unlike BuildBinTable()'s own per-cell calls
//   (aniso_axi.cpp, repsPerCell at a time, to keep its scratch buffer
//   small) -- equivalent per representative either way, since
//   RatioRun()'s loop carries no state across iterations (no reduction,
//   no element depends on any other), so index i's result is the same
//   whichever [lo,hi) range containing it was passed.
//
template <bool AhGeAz>
static void CheckAxiRatioRunMatchesScalar(const ScatterParams & sp) {

  const S2::TesselSphere & toa7 = UseTOADegree(7);

  S2::TesselSphere repTess(S2::TESS_ICO, AxiProposal::kRepDegree);
  const Index numReps  = (Index) repTess.size();
  const Index repBlock = (Index) toa7.size() / numReps;
  CHECK(repBlock * numReps == (Index) toa7.size());
  CHECK(repBlock > 1);  // degree 7 > kRepDegree: radii are nonzero.

  std::vector<R3::XYZ> centers(numReps);
  std::vector<Real> radius(numReps, 0.0);
  for (Index i = 0; i < numReps; ++i) {
    centers[i] = repTess[i].XYZ();
    const Index base = i * repBlock;
    Real maxAngle = 0.0;
    for (Index k = base; k < base + repBlock; ++k) {
      Real c = centers[i].Dot(toa7[k].XYZ());
      if (c > 1.0)  c = 1.0;
      if (c < -1.0) c = -1.0;
      const Real ang = acos(c);
      if (ang > maxAngle) {maxAngle = ang;}
    }
    radius[i] = maxAngle;
  }

  const Real aMin = sp.GetA();
  const Real ah   = sp.GetAx();
  const Real az   = sp.GetAz();
  CHECK(AhGeAz == (ah >= az));

  std::vector<AxiProposal::BoundRep> repsPP(numReps), repsPS(numReps),
                                     repsSP(numReps), repsSS(numReps);
  for (Index i = 0; i < numReps; ++i) {
    repsPP[i] = AxiProposal::PrepareBoundRep(
      centers[i], radius[i], sp.GetKP(), sp.GetKP(), aMin);
    repsPS[i] = AxiProposal::PrepareBoundRep(
      centers[i], radius[i], sp.GetKP(), sp.GetKS(), aMin);
    repsSP[i] = AxiProposal::PrepareBoundRep(
      centers[i], radius[i], sp.GetKS(), sp.GetKP(), aMin);
    repsSS[i] = AxiProposal::PrepareBoundRep(
      centers[i], radius[i], sp.GetKS(), sp.GetKS(), aMin);
  }

  AxiProposal::RepSoA soaPP, soaPS, soaSP, soaSS;
  soaPP.Assign(repsPP); soaPS.Assign(repsPS);
  soaSP.Assign(repsSP); soaSS.Assign(repsSS);

  std::vector<Real> ratio0(numReps), ratio1(numReps);
  long long checked = 0;

  // P bins: theta only (u_c/delta exactly as the AxiProposal
  // constructor's own P-bin loop), 'repsPP'/'repsPS'.
  for (int b = 0; b < AxiProposal::kThetaBins; ++b) {
    const Real theta_lo = (2.0*b)     * Geometry::DtoR;
    const Real theta_hi = (2.0*b+2.0) * Geometry::DtoR;
    const Real theta_c  = 0.5*(theta_lo+theta_hi);
    const Real delta    = 0.5*(theta_hi-theta_lo);
    const R3::XYZ u_c(-sin(theta_c), 0.0, cos(theta_c));
    const AxiProposal::BinCap cap = AxiProposal::MakeBinCap(u_c, delta);

    AxiProposal::RatioRun<AhGeAz>(soaPP, 0, numReps, cap, ah, az, &ratio0[0]);
    AxiProposal::RatioRun<AhGeAz>(soaPS, 0, numReps, cap, ah, az, &ratio1[0]);

    for (Index i = 0; i < numReps; ++i) {
      const Real wantPP = AxiProposal::BoundRatio(repsPP[i], cap, ah, az);
      CHECK_MSG(ratio0[i] == wantPP,
               "P bin=" << b << " PP rep=" << i << ": RatioRun="
               << ratio0[i] << " BoundRatio=" << wantPP);
      const Real wantPS = AxiProposal::BoundRatio(repsPS[i], cap, ah, az);
      CHECK_MSG(ratio1[i] == wantPS,
               "P bin=" << b << " PS rep=" << i << ": RatioRun="
               << ratio1[i] << " BoundRatio=" << wantPS);
      checked += 2;
    }
  }

  // S bins: theta and chi (same formulas as the constructor's S-bin
  // loop), 'repsSP'/'repsSS'.
  for (int b = 0; b < AxiProposal::kThetaBins; ++b) {
    for (int c = 0; c < AxiProposal::kChiBins; ++c) {
      const Real theta_lo = (2.0*b)     * Geometry::DtoR;
      const Real theta_hi = (2.0*b+2.0) * Geometry::DtoR;
      const Real theta_c  = 0.5*(theta_lo+theta_hi);
      const Real chi_lo = (10.0*c)      * Geometry::DtoR;
      const Real chi_hi = (10.0*c+10.0) * Geometry::DtoR;
      const Real chi_c  = 0.5*(chi_lo+chi_hi);
      const Real delta = 0.5*(theta_hi-theta_lo)
                        + sin(theta_hi)*0.5*(chi_hi-chi_lo);
      const R3::XYZ u_c(-sin(theta_c)*cos(chi_c),
                        sin(theta_c)*sin(chi_c),
                        cos(theta_c));
      const AxiProposal::BinCap cap = AxiProposal::MakeBinCap(u_c, delta);

      AxiProposal::RatioRun<AhGeAz>(soaSP, 0, numReps, cap, ah, az, &ratio0[0]);
      AxiProposal::RatioRun<AhGeAz>(soaSS, 0, numReps, cap, ah, az, &ratio1[0]);

      for (Index i = 0; i < numReps; ++i) {
        const Real wantSP = AxiProposal::BoundRatio(repsSP[i], cap, ah, az);
        CHECK_MSG(ratio0[i] == wantSP,
                 "S bin=" << b << "," << c << " SP rep=" << i
                 << ": RatioRun=" << ratio0[i] << " BoundRatio=" << wantSP);
        const Real wantSS = AxiProposal::BoundRatio(repsSS[i], cap, ah, az);
        CHECK_MSG(ratio1[i] == wantSS,
                 "S bin=" << b << "," << c << " SS rep=" << i
                 << ": RatioRun=" << ratio1[i] << " BoundRatio=" << wantSS);
        checked += 2;
      }
    }
  }

  CHECK(checked == 2LL * (AxiProposal::kThetaBins
                         + AxiProposal::kThetaBins*AxiProposal::kChiBins)
                  * numReps);

}

TEST(axi_ratio_run_matches_scalar) {
  CheckAxiRatioRunMatchesScalar<true>(LopNorCrust(2, 2, 0.5));     // ah >= az
  CheckAxiRatioRunMatchesScalar<false>(LopNorCrust(0.5, 0.5, 2));  // ah <  az
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 6 (PERFORMANCE PLAN) TESTS: Near-axisymmetric          ****
// ****  validation                                                  ****
// ****                                                              ****
//
// Validates the theory behind decisions 4-5 (the faster axisymmetric
// sampler and its same-medium performance gate) by checking that a
// *near*-axisymmetric medium (ay = ax*(1+delta), delta small) behaves
// consistently with the axisymmetric (ay == ax) case: the exact rate
// varies smoothly (linearly, to leading order) in delta, the
// axisymmetric and general samplers -- via AUTO, so the routing itself
// is exercised too -- agree statistically on a medium so close to
// axisymmetric that they are effectively sampling the same physics, the
// two GetReportMFPs() branches (AXI's 8 states, ANI's 9) agree on the
// states that coincide by the axisymmetric spec's symmetry, and the
// per-real-event cost (tentative collisions and table draws) of both
// samplers is recorded for a range of incidence states. See
// docs/superpowers/specs/2026-10-17-axisymmetric-sampler-performance-
// design.md, decision 6, and .superpowers/sdd/2026-10-17-axisymmetric-
// sampler-performance/task-6-brief.md.
//

TEST(near_axi_rate_matches) {

  const S2::TesselSphere & toas = UseTOADegree(5);

  struct SpecCase { Real ax, ay, az; };
  const SpecCase specs[] = {
    {2, 2, 0.5}, {5, 5, 0.5}, {0.5, 0.5, 2},
  };
  const Real kappas[] = {0.3, 0.5};
  const Real deltas[] = {1e-6, 1e-4, 1e-2};

  const R3::XYZ E(1,0,0), N(0,1,0), U(0,0,1), D(Dir(60,30));
  struct Case { raytype in; R3::XYZ dir; R3::XYZ pol; };
  const Case cases[] = {
    {RAY_P, E, E}, {RAY_P, N, N}, {RAY_P, U, U}, {RAY_P, D, D},
    {RAY_S, E, N}, {RAY_S, E, U},
    {RAY_S, N, E}, {RAY_S, N, U},
    {RAY_S, U, E}, {RAY_S, U, N},
  };
  const std::size_t nCases  = sizeof(cases)/sizeof(cases[0]);
  const std::size_t nDeltas = sizeof(deltas)/sizeof(deltas[0]);

  for (std::size_t si = 0; si < sizeof(specs)/sizeof(specs[0]); ++si) {
    for (std::size_t ki = 0; ki < sizeof(kappas)/sizeof(kappas[0]); ++ki) {

      ScatterParams spBase = LopNorCrust(specs[si].ax, specs[si].ay,
                                         specs[si].az, kappas[ki]);
      Scatterer scatBase(spBase, Scatterer::SAMPLER_GENERAL);
                            // SAMPLER_GENERAL forced: IntegratedRate()
                            // is a pure function of mParams and 'toas'
                            // (IntegrateRateAndDipole() -> GSATO()/
                            // AnisoRatio()), independent of which
                            // sampler the Scatterer would use to draw
                            // events, so forcing GENERAL (rather than
                            // letting AUTO pick AXISYMMETRIC for these
                            // ax==ay specs) just avoids paying for an
                            // unused AxiProposal build.

      std::vector<Real> R0(nCases);
      for (std::size_t ci = 0; ci < nCases; ++ci) {
        R0[ci] = scatBase.IntegratedRate(cases[ci].in, cases[ci].dir,
                                         cases[ci].pol, toas);
      }

      std::vector<std::vector<Real> > dR(nCases, std::vector<Real>(nDeltas));
      for (std::size_t di = 0; di < nDeltas; ++di) {
        ScatterParams spNear = LopNorCrust(specs[si].ax,
                                           specs[si].ax * (1.0 + deltas[di]),
                                           specs[si].az, kappas[ki]);
        Scatterer scatNear(spNear, Scatterer::SAMPLER_GENERAL);

        for (std::size_t ci = 0; ci < nCases; ++ci) {
          const Real Rn = scatNear.IntegratedRate(cases[ci].in, cases[ci].dir,
                                                  cases[ci].pol, toas);
          dR[ci][di] = Rn - R0[ci];
          const Real bound = 5.0 * deltas[di] * R0[ci];
          CHECK_MSG(std::fabs(dR[ci][di]) <= bound,
                   "spec " << si << " kappa " << kappas[ki] << " case " << ci
                   << " delta " << deltas[di] << ": |dR|="
                   << std::fabs(dR[ci][di]) << " bound=" << bound);
        }
      }

      // Where the smallest-delta signal is resolvable above rounding
      // noise, confirm it scales close to linearly: a 100x increase in
      // delta (1e-6 -> 1e-4) should give roughly a 100x increase in
      // |dR|, loosely bounded to [50,200] per the brief.
      for (std::size_t ci = 0; ci < nCases; ++ci) {
        if (std::fabs(dR[ci][0]) > 1e-13 * R0[ci]) {
          const Real ratio = dR[ci][1] / dR[ci][0];
          CHECK_MSG(ratio >= 50.0 && ratio <= 200.0,
                   "spec " << si << " kappa " << kappas[ki] << " case " << ci
                   << ": dR(1e-4)/dR(1e-6)=" << ratio
                   << " expected in [50,200]");
        }
      }

    }
  }

}


TEST(near_axi_sampler_matches) {

  UseTOADegree(5);
  const Real delta = 1e-5;
  ScatterParams spAxi = LopNorCrust(2, 2, 0.5);
  ScatterParams spGen = LopNorCrust(2, 2.0 * (1.0 + delta), 0.5);
  Scatterer axiScat(spAxi, Scatterer::SAMPLER_AUTO);
  Scatterer genScat(spGen, Scatterer::SAMPLER_AUTO);

  CHECK_MSG(axiScat.GetReportMFPs().size() == 8,
           "AUTO should route ax==ay to the 8-state AXI branch, got "
           << axiScat.GetReportMFPs().size());
  CHECK_MSG(genScat.GetReportMFPs().size() == 9,
           "AUTO should route ax!=ay to the 9-state ANI branch, got "
           << genScat.GetReportMFPs().size());

  const R3::XYZ E(1,0,0), N(0,1,0), U(0,0,1), D(Dir(60,30));
  struct Case { raytype in; R3::XYZ dir; R3::XYZ pol; };
  const Case cases[] = {
    {RAY_P, E, E}, {RAY_P, N, N}, {RAY_P, U, U}, {RAY_P, D, D},
    {RAY_S, E, N}, {RAY_S, E, U},
    {RAY_S, N, E}, {RAY_S, N, U},
    {RAY_S, U, E}, {RAY_S, U, N},
  };
  const int n = 20000;

  for (std::size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i) {
    const S2::ThetaPhi dirTP(cases[i].dir.Theta(), cases[i].dir.Phi());
    const Real meanA = MeanDistanceToRealEvent(axiScat, cases[i].in, dirTP,
                                               cases[i].pol, n);
    const Real meanG = MeanDistanceToRealEvent(genScat, cases[i].in, dirTP,
                                               cases[i].pol, n);
    const Real sigA = meanA / std::sqrt((Real) n);
    const Real sigG = meanG / std::sqrt((Real) n);
    const Real tol = 4.0 * std::sqrt(sigA*sigA + sigG*sigG);
    CHECK_MSG(std::fabs(meanA - meanG) < tol,
             "case " << i << ": meanA=" << meanA << " meanG=" << meanG
             << " tol=" << tol);
  }

  // (Conv, degree-2 cell) two-sample chi-square -- the same unequal-N
  // statistic axi_matches_general() uses above -- for P at Dir(60,30)
  // and S along E with polarization U:
  int divisor = 1;
  for (int i = 0; i < 5-2; ++i) divisor *= 4;  // 4^(5-2)=64

  const S2::ThetaPhi dirP(60.0*Geometry::DtoR, 30.0*Geometry::DtoR);
  const S2::ThetaPhi dirS(Geometry::Pi90, 0.0);  // E
  const R3::XYZ polU(0, 0, 1);                   // U

  struct DistCase { raytype in; S2::ThetaPhi dir; R3::XYZ pol; };
  const DistCase distCases[] = {
    {RAY_P, dirP, dirP.XYZ()},
    {RAY_S, dirS, polU},
  };
  const int targetReal = 200000;

  std::map<std::pair<int,int>, int> countsA, countsG;

  for (std::size_t i = 0; i < sizeof(distCases)/sizeof(distCases[0]); ++i) {
    const Real rot = Scatterer::RotationForS1(distCases[i].dir,
                                              distCases[i].pol);
    const R3::XYZ origin(0,0,0);
    int collectedA = 0, collectedG = 0;
    while (collectedA < targetReal || collectedG < targetReal) {
      if (collectedA < targetReal) {
        Scatterer::AnisoAttempt a = axiScat.AttemptAnisotropicScatter(
          distCases[i].in, origin, distCases[i].dir, rot);
        if (a.IsReal) {
          countsA[std::make_pair((int)i*8+a.Conv, (int)(a.Toa/divisor))]++;
          ++collectedA;
        }
      }
      if (collectedG < targetReal) {
        Scatterer::AnisoAttempt g = genScat.AttemptAnisotropicScatter(
          distCases[i].in, origin, distCases[i].dir, rot);
        if (g.IsReal) {
          countsG[std::make_pair((int)i*8+g.Conv, (int)(g.Toa/divisor))]++;
          ++collectedG;
        }
      }
    }
  }

  Real Na = 0, Ng = 0;
  for (std::map<std::pair<int,int>,int>::const_iterator it = countsA.begin();
       it != countsA.end(); ++it) {Na += it->second;}
  for (std::map<std::pair<int,int>,int>::const_iterator it = countsG.begin();
       it != countsG.end(); ++it) {Ng += it->second;}

  Real chi2 = 0; int dof = 0;
  std::map<std::pair<int,int>,int> keys = countsA;
  for (std::map<std::pair<int,int>,int>::const_iterator it = countsG.begin();
       it != countsG.end(); ++it) {keys[it->first] = 0;}
  for (std::map<std::pair<int,int>,int>::const_iterator it = keys.begin();
       it != keys.end(); ++it) {
    const Real a = countsA.count(it->first) ? countsA[it->first] : 0;
    const Real g = countsG.count(it->first) ? countsG[it->first] : 0;
    if (a+g < 20.0) continue;
    const Real num = a*Ng - g*Na;
    chi2 += num*num/(Na*Ng*(a+g));
    ++dof;
  }
  CHECK_MSG(dof > 0, "no bins had combined count >= 20");
  const Real chi2_per_dof = chi2/dof;
  CHECK_MSG(chi2_per_dof < 1.5,
           "chi2/dof too large: " << chi2_per_dof << " (dof=" << dof << ")");

}


TEST(near_axi_report_matches) {

  UseTOADegree(5);  // GetReportMFPs() always integrates over its own
                    // fixed degree-6 ReportDirs() regardless (see
                    // report_directional_mfps, above).
  const Real delta = 1e-5;
  ScatterParams spAxi  = LopNorCrust(2, 2, 0.5);
  ScatterParams spNear = LopNorCrust(2, 2.0 * (1.0 + delta), 0.5);
  Scatterer axiScat(spAxi,  Scatterer::SAMPLER_AUTO);
  Scatterer genScat(spNear, Scatterer::SAMPLER_AUTO);

  const std::vector<Real> & axi = axiScat.GetReportMFPs();
  const std::vector<Real> & ani = genScat.GetReportMFPs();
  CHECK_MSG(axi.size() == 8, "expected 8 AXI directional MFPs, got "
           << axi.size());
  CHECK_MSG(ani.size() == 9, "expected 9 ANI directional MFPs, got "
           << ani.size());

  // AXI order: P0 P45 P90 S0 S45SV S45SH S90SV S90SH
  // ANI order: P_E P_N P_U S_E/N S_E/U S_N/E S_N/U S_U/E S_U/N
  //
  // The axisymmetric spec's symmetry axis is z (up), so AXI's
  // "incidence angle from the axis" states map onto ANI's "incidence
  // along a named fabric axis" states as follows: propagation along the
  // axis (AXI theta=0) is ANI's U axis; propagation perpendicular to the
  // axis (AXI theta=90) is ANI's E and N axes (indistinguishable by
  // symmetry when ax==ay); and for S90 (propagation perpendicular to the
  // axis), "SV" (polarization containing the axis) and "SH"
  // (polarization perpendicular to the axis) map onto whichever of
  // ANI's two propagation/polarization axis pairs puts the polarization
  // along/across U respectively.
  const Real P0    = axi[0], P90   = axi[2], S0 = axi[3],
             S90SV = axi[6], S90SH = axi[7];
  const Real P_E   = ani[0], P_N   = ani[1], P_U   = ani[2],
             S_E_N = ani[3], S_E_U = ani[4], S_N_E = ani[5],
             S_N_U = ani[6], S_U_E = ani[7], S_U_N = ani[8];

  CHECK_MSG(NearRel(P0, P_U, 1e-4), "P0=" << P0 << " P_U=" << P_U);
  CHECK_MSG(NearRel(P90, P_E, 1e-4), "P90=" << P90 << " P_E=" << P_E);
  CHECK_MSG(NearRel(P90, P_N, 1e-4), "P90=" << P90 << " P_N=" << P_N);
  CHECK_MSG(NearRel(S0, S_U_E, 1e-4), "S0=" << S0 << " S_U_E=" << S_U_E);
  CHECK_MSG(NearRel(S0, S_U_N, 1e-4), "S0=" << S0 << " S_U_N=" << S_U_N);
  CHECK_MSG(NearRel(S90SV, S_E_U, 1e-4),
           "S90SV=" << S90SV << " S_E_U=" << S_E_U);
  CHECK_MSG(NearRel(S90SV, S_N_U, 1e-4),
           "S90SV=" << S90SV << " S_N_U=" << S_N_U);
  CHECK_MSG(NearRel(S90SH, S_E_N, 1e-4),
           "S90SH=" << S90SH << " S_E_N=" << S_E_N);
  CHECK_MSG(NearRel(S90SH, S_N_E, 1e-4),
           "S90SH=" << S90SH << " S_N_E=" << S_N_E);

}


// ::: MeasureCollisionCost(s, in, dir, pol, n, variance, meanDraws) :::
//
//   Chases n real (accepted) events exactly as MeanDistanceToRealEvent()
//   does -- one GetRandomPathLength() call per tentative collision,
//   immediately followed by AttemptAnisotropicScatter() -- but records
//   the per-event tentative-collision count instead of distance,
//   returning its mean and, via 'variance', the sample variance needed
//   to compare two samplers' means against their combined standard
//   error. 'meanDraws' receives the mean of AnisoAttempt.Draws summed
//   over the same tentative collisions.
//
//   A "tentative collision" is one GetRandomPathLength() draw (one
//   majorant-distance proposal; see AttemptAnisotropicScatter()'s
//   header comment -- without --overridemfp, each proposal/accept trial
//   corresponds to exactly one such draw). A "draw" is
//   AnisoAttempt.Draws: for the general sampler this is always 1 (every
//   tentative collision reaches the (conversion, take-off-angle) table
//   draw), but for the axisymmetric sampler it is 0 on a tentative
//   collision thinned by the cheap top-level bin-rate test -- see
//   AttemptAxisymmetricScatter() -- and 1 otherwise. So "draws" can run
//   well below "tentative collisions" for the axisymmetric sampler: that
//   gap is the performance effect test_near_axi_cost's cost table
//   documents, while the two stay equal for the general sampler.
//
static Real MeasureCollisionCost(Scatterer & s, raytype in,
                                 const S2::ThetaPhi & dir,
                                 const R3::XYZ & pol, int n,
                                 Real & variance, Real & meanDraws) {
  const Real rot = Scatterer::RotationForS1(dir, pol);
  const R3::XYZ origin(0, 0, 0);
  std::vector<Real> tentative(n);
  long long totalDraws = 0;

  for (int i = 0; i < n; ++i) {
    int tcount = 0;
    Scatterer::AnisoAttempt att;
    do {
      s.GetRandomPathLength(in);
      att = s.AttemptAnisotropicScatter(in, origin, dir, rot);
      ++tcount;
      totalDraws += att.Draws;
    } while (!att.IsReal);
    tentative[i] = (Real) tcount;
  }

  Real sum = 0;
  for (int i = 0; i < n; ++i) sum += tentative[i];
  const Real mean = sum / n;

  Real sqsum = 0;
  for (int i = 0; i < n; ++i) {
    const Real d = tentative[i] - mean;
    sqsum += d * d;
  }
  variance = sqsum / (n - 1);
  meanDraws = (Real) totalDraws / n;
  return mean;
}


TEST(near_axi_cost) {

  UseTOADegree(7);
  Scatterer::SetBuildThreads(8);  // Degree 7 builds several general-
                                  // sampler Scatterers below; table
                                  // construction is bit-identical for
                                  // any thread count (design decision
                                  // 2), so this only shortens the test.

  const Real delta = 1e-5;

  struct SpecCase { Real ax, ay, az; };
  const SpecCase specs[] = {
    {1, 1, 0.5}, {2, 2, 0.5}, {5, 5, 0.5},
  };

  struct StateCase { raytype in; Real theta; Real chi; };
  // P at theta 0,30,60,89; S at theta 30,89 with chi 0,45,90 -- the same
  // incidence grid as test_axi_draws_per_event, above.
  const StateCase states[] = {
    {RAY_P, 0, 0}, {RAY_P, 30, 0}, {RAY_P, 60, 0}, {RAY_P, 89, 0},
    {RAY_S, 30, 0}, {RAY_S, 30, 45}, {RAY_S, 30, 90},
    {RAY_S, 89, 0}, {RAY_S, 89, 45}, {RAY_S, 89, 90},
  };

  const int n = 5000;

  std::ofstream out("build/aniso_near_axi_cost.txt");
  out << "# spec state raytype theta chi sampler mean_tentative mean_draws\n";

  for (std::size_t si = 0; si < sizeof(specs)/sizeof(specs[0]); ++si) {

    ScatterParams spAxi  = LopNorCrust(specs[si].ax, specs[si].ay,
                                       specs[si].az);
    ScatterParams spNear = LopNorCrust(specs[si].ax,
                                       specs[si].ay * (1.0 + delta),
                                       specs[si].az);
    Scatterer axiScat(spAxi,  Scatterer::SAMPLER_AXISYMMETRIC);
    Scatterer genScat(spNear, Scatterer::SAMPLER_GENERAL);

    // Deterministic majorant-MFP check (no RNG involved): the
    // axisymmetric sampler's bin-tailored majorant must not be looser
    // (smaller MFP) than the near-axisymmetric general sampler's, to
    // within the near medium's own 1e-5 asymmetry.
    const Real axiMfpP = axiScat.GetMajorantMeanFreePath(RAY_P);
    const Real genMfpP = genScat.GetMajorantMeanFreePath(RAY_P);
    CHECK_MSG(axiMfpP >= genMfpP * (1.0 - 1e-5),
             "spec " << si << ": axi P majorant MFP " << axiMfpP
             << " < near-general " << genMfpP << " * (1-1e-5)");
    const Real axiMfpS = axiScat.GetMajorantMeanFreePath(RAY_S);
    const Real genMfpS = genScat.GetMajorantMeanFreePath(RAY_S);
    CHECK_MSG(axiMfpS >= genMfpS * (1.0 - 1e-5),
             "spec " << si << ": axi S majorant MFP " << axiMfpS
             << " < near-general " << genMfpS << " * (1-1e-5)");

    for (std::size_t ti = 0; ti < sizeof(states)/sizeof(states[0]); ++ti) {

      const S2::ThetaPhi dir(states[ti].theta * Geometry::DtoR, 0.0);
      R3::XYZ pol;
      if (states[ti].in == RAY_P) {
        pol = dir.XYZ();
      } else {
        const R3::OrthoAxes rotAA(dir.Theta(), dir.Phi(),
                                  states[ti].chi * Geometry::DtoR);
        pol = rotAA.S1();
      }

      Real varA = 0, drawsA = 0, varG = 0, drawsG = 0;
      const Real meanA = MeasureCollisionCost(axiScat, states[ti].in, dir,
                                              pol, n, varA, drawsA);
      const Real meanG = MeasureCollisionCost(genScat, states[ti].in, dir,
                                              pol, n, varG, drawsG);

      const char * rtName = (states[ti].in == RAY_P) ? "P" : "S";
      out << "(" << specs[si].ax << "," << specs[si].ay << "," << specs[si].az
          << ") " << ti << " " << rtName << " " << states[ti].theta << " "
          << states[ti].chi << " AXI " << meanA << " " << drawsA << "\n";
      out << "(" << specs[si].ax << "," << specs[si].ay << "," << specs[si].az
          << ") " << ti << " " << rtName << " " << states[ti].theta << " "
          << states[ti].chi << " GENERAL " << meanG << " " << drawsG << "\n";

      const Real sigmaDiff = std::sqrt(varA/n + varG/n);
      CHECK_MSG(meanA <= meanG + 4.0*sigmaDiff,
               "spec " << si << " state " << ti << " (" << rtName
               << " theta=" << states[ti].theta << " chi=" << states[ti].chi
               << "): axi mean tentative " << meanA << " exceeds general "
               << meanG << " + 4*sigma=" << (4.0*sigmaDiff));
    }
  }

  out.close();
  Scatterer::SetBuildThreads(1);  // restore

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TASK 6b: Canonical P frame in the general sampler           ****
// ****                                                              ****
//
// Investigation (2026-10-17-axisymmetric-sampler-performance, ruling
// R10): for P incidence the equal-weight take-off-angle grid is not
// rotation-invariant, so the integrated rate depends on the (unphysical)
// P frame angle about the ray. The axisymmetric sampler always uses the
// canonical frame (CanonicalPPolarization()); before this task the
// general sampler used whatever polarization-rotation angle the
// incoming phonon happened to carry, so already-scattered P phonons hit
// the two samplers with different effective rates. See
// docs/toa-quadrature-bias.md for the mechanism and example magnitudes.
//

// ::: TwoSampleChiSquare(A, G, dof) :::
//
//   The unequal-N two-sample chi-square statistic test_axi_matches_
//   general() and test_near_axi_sampler_matches() use inline, factored
//   out here for test_p_frame_consistency(): chi2 = sum_i (a_i*Ng -
//   g_i*Na)^2 / (Na*Ng*(a_i+g_i)), summed over bins whose combined count
//   a_i+g_i is at least 20. Returns the total chi2 and sets 'dof' to the
//   number of bins included.
//
static Real TwoSampleChiSquare(const std::map<std::pair<int,int>,int> & A,
                               const std::map<std::pair<int,int>,int> & G,
                               int & dof) {
  Real Na = 0, Ng = 0;
  for (std::map<std::pair<int,int>,int>::const_iterator it = A.begin();
       it != A.end(); ++it) {Na += it->second;}
  for (std::map<std::pair<int,int>,int>::const_iterator it = G.begin();
       it != G.end(); ++it) {Ng += it->second;}

  Real chi2 = 0; dof = 0;
  std::map<std::pair<int,int>,int> keys = A;
  for (std::map<std::pair<int,int>,int>::const_iterator it = G.begin();
       it != G.end(); ++it) {keys[it->first] = 0;}
  for (std::map<std::pair<int,int>,int>::const_iterator it = keys.begin();
       it != keys.end(); ++it) {
    const Real a = A.count(it->first) ? A.at(it->first) : 0;
    const Real g = G.count(it->first) ? G.at(it->first) : 0;
    if (a+g < 20.0) continue;
    const Real num = a*Ng - g*Na;
    chi2 += num*num/(Na*Ng*(a+g));
    ++dof;
  }
  return chi2;
}


TEST(p_frame_consistency) {

  UseTOADegree(7);  // Matches the investigation's diagnostic (debug-
                    // scratch.patch, scratch_lopnor_r10_abs).

  // Lop Nor r=10-like crust and mantle (task-6b-brief.md, Step 1):
  struct Layer {const char * name; Real vp, vs, a, az, kappa;};
  const Layer layers[] = {
    {"crust",  6.40, 3.63, 5.0, 0.5, 0.3},
    {"mantle", 8.07, 4.63, 7.0, 0.7, 0.5},
  };

  const Real thetasDeg[] = {30.0, 60.0, 89.0};
  const Real phiDeg = 20.0;        // Arbitrary, non-axis-aligned, as in
                                   // the investigation's diagnostic.
  const Real rot = 0.7;            // Non-canonical carried P
                                   // polarization-rotation angle (rad);
                                   // the investigation's diagnostic
                                   // tested P only at this angle.
  const int nDist = 20000;
  const int targetReal = 100000;

  const R3::XYZ origin(0, 0, 0);
  S2::TesselSphere deg2(S2::TESS_ICO, 2);  // For binning absolute
                              // outgoing directions, independent of the
                              // degree-7 array UseTOADegree(7) pointed
                              // PhononSource at above -- built directly
                              // rather than through UseTOADegree() so it
                              // never becomes the global TOA array (that
                              // would desync it from the degree-7
                              // mPDists/mpAxi tables the scatterers
                              // below are built against).

  for (std::size_t li = 0; li < sizeof(layers)/sizeof(layers[0]); ++li) {
    const Layer & L = layers[li];

    ScatterParams spAxi = LopNorCrust(L.a, L.a, L.az, L.kappa, 0.01,
                                      L.vp, L.vs);                // ax==ay
    ScatterParams spGen = LopNorCrust(L.a, L.a*(1.0+1.0e-5), L.az, L.kappa,
                                      0.01, L.vp, L.vs);          // ax!=ay
    Scatterer axiScat(spAxi, Scatterer::SAMPLER_AUTO);
    Scatterer genScat(spGen, Scatterer::SAMPLER_AUTO);
    CHECK_MSG(axiScat.GetReportMFPs().size() == 8,
             L.name << ": expected AUTO to route ax==ay to the 8-state "
             "AXI branch, got " << axiScat.GetReportMFPs().size());
    CHECK_MSG(genScat.GetReportMFPs().size() == 9,
             L.name << ": expected AUTO to route ax!=ay to the 9-state "
             "GENERAL branch, got " << genScat.GetReportMFPs().size());

    for (std::size_t ti = 0; ti < sizeof(thetasDeg)/sizeof(thetasDeg[0]);
         ++ti) {
      const S2::ThetaPhi dirTP(thetasDeg[ti]*Geometry::DtoR,
                               phiDeg*Geometry::DtoR);

      // A 'pol' vector whose RotationForS1(dirTP, pol) is 'rot': the
      // same construction test_axi_draws_per_event() and others use to
      // hand a chosen rotation angle to the R3::XYZ-taking helpers.
      const R3::OrthoAxes polAA(dirTP.Theta(), dirTP.Phi(), rot);
      const R3::XYZ pol = polAA.S1();

      // --- mean distance to a real event, axi vs general ---
      const Real meanA = MeanDistanceToRealEvent(axiScat, RAY_P, dirTP,
                                                 pol, nDist);
      const Real meanG = MeanDistanceToRealEvent(genScat, RAY_P, dirTP,
                                                 pol, nDist);
      const Real sigA = meanA / std::sqrt((Real) nDist);
      const Real sigG = meanG / std::sqrt((Real) nDist);
      const Real tol = 4.0 * std::sqrt(sigA*sigA + sigG*sigG);
      CHECK_MSG(std::fabs(meanA - meanG) < tol,
               L.name << " theta=" << thetasDeg[ti] << ": meanA=" << meanA
               << " meanG=" << meanG << " tol=" << tol << " z="
               << (meanA-meanG)/std::sqrt(sigA*sigA+sigG*sigG));

      // --- joint (Conv, absolute outgoing direction cell) two-sample
      //     chi-square, applying InPol/Rel through Phonon::Transform()
      //     exactly as Phonon::Propagate() does on a real event ---
      std::map<std::pair<int,int>, int> countsA, countsG;
      int collectedA = 0, collectedG = 0;
      while (collectedA < targetReal || collectedG < targetReal) {
        if (collectedA < targetReal) {
          const Scatterer::AnisoAttempt a =
            axiScat.AttemptAnisotropicScatter(RAY_P, origin, dirTP, rot);
          if (a.IsReal) {
            Phonon ph(dirTP, RAY_P);
            ph.SetPolarization(a.InPol);
            ph.Transform(a.Rel);
            const R3::XYZ e = ph.GetDirection().XYZ();
            const int cell =
              deg2.GetBestIndexFromPoint(S2::ThetaPhi(e.Theta(), e.Phi()));
            countsA[std::make_pair(a.Conv, cell)]++;
            ++collectedA;
          }
        }
        if (collectedG < targetReal) {
          const Scatterer::AnisoAttempt g =
            genScat.AttemptAnisotropicScatter(RAY_P, origin, dirTP, rot);
          if (g.IsReal) {
            Phonon ph(dirTP, RAY_P);
            ph.SetPolarization(g.InPol);
            ph.Transform(g.Rel);
            const R3::XYZ e = ph.GetDirection().XYZ();
            const int cell =
              deg2.GetBestIndexFromPoint(S2::ThetaPhi(e.Theta(), e.Phi()));
            countsG[std::make_pair(g.Conv, cell)]++;
            ++collectedG;
          }
        }
      }

      int dof = 0;
      const Real chi2 = TwoSampleChiSquare(countsA, countsG, dof);
      CHECK_MSG(dof > 0, L.name << " theta=" << thetasDeg[ti]
               << ": no bins had combined count >= 20");
      const Real chi2_per_dof = chi2 / dof;
      CHECK_MSG(chi2_per_dof < 1.5,
               L.name << " theta=" << thetasDeg[ti]
               << ": chi2/dof too large: " << chi2_per_dof
               << " (dof=" << dof << ")");
    }
  }

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TOA QUADRATURE TESTS                                        ****
// ****                                                              ****
//
// See docs/superpowers/plans/2026-10-25-toa-area-quadrature.md and
// docs/superpowers/specs/2026-10-25-toa-area-quadrature-design.md, and
// .superpowers/sdd/2026-10-25-toa-area-quadrature/task-1-brief.md (Task
// 1: per-direction solid angles in S2::TesselSphere).
//

TEST(tessel_solid_angles) {

  // S2::TesselSphere's icosahedral quadrisection spreads take-off
  // directions over the sphere, but not perfectly uniformly: each
  // direction is really the centerpoint of a small spherical triangle,
  // and those triangles' solid angles (areas) vary somewhat from leaf
  // to leaf. This test builds, at every degree from 0 through 8, both a
  // plain TesselSphere and one with per-direction solid angles turned
  // on, and checks:
  //
  //   o  the plain one (withSolidAngles defaulted to false) is
  //      completely unaffected: HasSolidAngles() is false and
  //      SolidAngles() is empty;
  //   o  the solid-angle-reporting one has exactly the same directions,
  //      in the same order, as the plain one -- ie, turning the
  //      feature on never perturbs the direction set itself;
  //   o  every reported solid angle is positive, and the full set sums
  //      to 4*pi steradians (the area of the whole unit sphere) to
  //      within tight relative tolerance;
  //   o  at degree 0 (20 congruent icosahedral faces, no subdivision
  //      yet), every leaf has the same solid angle to near machine
  //      precision;
  //   o  at degrees 6 through 8, the leaf-to-leaf spread of solid
  //      angles (relative to the mean) matches the fixed statistics a
  //      throwaway prototype measured at degree 6 (task-1-brief.md) --
  //      empirical confirmation that the leaf formula in
  //      Triangle::populate_lists() (geom_s2.cpp) is being applied
  //      correctly and consistently across the whole recursive
  //      subdivision, not just at the one degree it was measured at.
  //
  // Per-degree summary statistics are also logged to
  // build/toa_quadrature_solid_angles.txt for human inspection.

  std::ofstream out("build/toa_quadrature_solid_angles.txt");
  out << "# degree N sum_rel_err min_over_mean max_over_mean rms\n";

  const Real fourPi = 4.0 * Geometry::Pi;

  for (int d = 0; d <= 8; ++d) {

    S2::TesselSphere plain(S2::TESS_ICO, d);      // withSolidAngles
    S2::TesselSphere sa(S2::TESS_ICO, d, true);   // defaults to false

    CHECK_MSG(!plain.HasSolidAngles() && plain.SolidAngles().empty(),
             "degree " << d << ": plain TesselSphere unexpectedly reports "
             "solid angles");
    CHECK_MSG(sa.HasSolidAngles()
             && sa.size() == plain.size()
             && sa.SolidAngles().size() == sa.size(),
             "degree " << d << ": sa.HasSolidAngles()=" << sa.HasSolidAngles()
             << " sa.size()=" << sa.size() << " plain.size()=" << plain.size()
             << " sa.SolidAngles().size()=" << sa.SolidAngles().size());

    Real sum = 0.0;
    for (std::size_t k = 0; k < sa.size(); ++k) {
      CHECK_MSG(sa[k].Theta() == plain[k].Theta()
               && sa[k].Phi() == plain[k].Phi(),
               "degree " << d << " k=" << k << ": direction differs between "
               "plain and solid-angle TesselSpheres");
      CHECK_MSG(sa.SolidAngles()[k] > 0.0,
               "degree " << d << " k=" << k << ": Omega="
               << sa.SolidAngles()[k]);
      sum += sa.SolidAngles()[k];
    }

    const Real mean = sum / (Real) sa.size();
    const Real sumRelErr = std::fabs(sum - fourPi) / fourPi;
    CHECK_MSG(sumRelErr <= 1e-12,
             "degree " << d << ": sum=" << sum << " expected " << fourPi
             << " relerr=" << sumRelErr);

    Real minOmega = sa.SolidAngles()[0];
    Real maxOmega = sa.SolidAngles()[0];
    Real sumSqDev = 0.0;
    for (std::size_t k = 0; k < sa.SolidAngles().size(); ++k) {
      const Real omega = sa.SolidAngles()[k];
      minOmega = std::min(minOmega, omega);
      maxOmega = std::max(maxOmega, omega);
      const Real dev = omega/mean - 1.0;
      sumSqDev += dev*dev;
    }
    const Real minOverMean = minOmega / mean;
    const Real maxOverMean = maxOmega / mean;
    const Real rms = std::sqrt(sumSqDev / (Real) sa.SolidAngles().size());

    out << d << " " << sa.size() << " " << sumRelErr << " "
        << minOverMean << " " << maxOverMean << " " << rms << "\n";

    if (d == 0) {
      // Twenty congruent icosahedral faces, undivided: every leaf
      // triangle is the same shape and size, so the raw min and max
      // (not just their ratios to the mean) should agree to near
      // machine precision.
      CHECK_MSG(maxOmega/minOmega - 1.0 <= 1e-14,
               "degree 0: minOmega=" << minOmega << " maxOmega=" << maxOmega);
    }
    if (d >= 6) {
      // Recursive quadrisection has run deep enough that the
      // leaf-shape distribution has settled into the fixed pattern the
      // throwaway prototype measured at degree 6 (task-1-brief.md):
      // min/mean=0.9277266, max/mean=1.2065026, rms=0.0864447.
      CHECK_MSG(minOverMean >= 0.9276 && minOverMean <= 0.9278,
               "degree " << d << ": min/mean=" << minOverMean);
      CHECK_MSG(maxOverMean >= 1.2064 && maxOverMean <= 1.2067,
               "degree " << d << ": max/mean=" << maxOverMean);
      CHECK_MSG(rms >= 0.0864 && rms <= 0.0865,
               "degree " << d << ": rms=" << rms);
    }
  }

  out.close();

}


TEST(tessel_solid_angles_mirror) {

  // The icosahedron's 12 vertices (init_icosahedron(), geom_s2.cpp) sit
  // symmetrically under reflection through any of the three coordinate
  // planes, so the whole take-off-angle set -- and, if the leaf formula
  // is being applied correctly, each direction's solid angle -- should
  // be invariant under those reflections too: a solid angle depends
  // only on the pairwise angles between its three vertex vectors (via
  // dot products) and the unsigned volume of the parallelepiped they
  // span (via the absolute value in the formula), and both of those are
  // preserved by any orthogonal map, reflections included. This checks
  // that at degree 3, for every direction, the images mirrored through
  // the y-z plane and through the x-z plane each land exactly on some
  // other direction in the set, with exactly the same solid angle.

  S2::TesselSphere sa(S2::TESS_ICO, 3, true);

  for (std::size_t k = 0; k < sa.size(); ++k) {
    const R3::XYZ p = sa[k].XYZ();
    const R3::XYZ images[2] = {
      R3::XYZ(-p.x(),  p.y(), p.z()),   // Mirror through the y-z plane
      R3::XYZ( p.x(), -p.y(), p.z()),   // Mirror through the x-z plane
    };
    for (int im = 0; im < 2; ++im) {
      const S2::ThetaPhi imageTP(images[im].Theta(), images[im].Phi());
      const int j = sa.GetBestIndexFromPoint(imageTP);
      CHECK_MSG(sa[j].XYZ().DistFrom(images[im]) < 1e-9,
               "k=" << k << " image=" << im << ": nearest direction j=" << j
               << " is " << sa[j].XYZ().DistFrom(images[im])
               << " away from the mirrored point");
      CHECK_MSG(std::fabs(sa.SolidAngles()[j] - sa.SolidAngles()[k])
               <= 1e-12 * sa.SolidAngles()[k],
               "k=" << k << " j=" << j << " image=" << im
               << ": Omega[k]=" << sa.SolidAngles()[k]
               << " Omega[j]=" << sa.SolidAngles()[j]);
    }
  }

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  TEST RUNNER                                                 ****
// ****                                                              ****
//

int main(int argc, char ** argv) {

  std::vector<TestCase> & tests = Registry();

  if (argc >= 2 && std::string(argv[1]) == "--list") {
    for (std::size_t i = 0; i < tests.size(); ++i) {
      std::cout << tests[i].name << "\n";
    }
    return 0;
  }

  const std::string filter = (argc >= 2) ? argv[1] : "";
  bool any_failed = false;

  for (std::size_t i = 0; i < tests.size(); ++i) {
    const std::string & name = tests[i].name;
    const bool isolated = (name.compare(0, 9, "isolated_") == 0);

    bool should_run;
    if (filter.empty()) {
      should_run = !isolated;                      // Default run: skip
                                                    // isolated_* tests.
    } else if (isolated) {
      should_run = (name == filter);                // Isolated tests only
                                                    // run on exact match.
    } else {
      should_run = (name.find(filter) != std::string::npos);
    }
    if (!should_run) continue;

    try {
      tests[i].func();
      std::cout << "PASS " << name << "\n";
    } catch (const TestFailure & f) {
      std::cout << "FAIL " << name << ": " << f.message << "\n";
      any_failed = true;
    } catch (const std::exception & e) {
      std::cout << "FAIL " << name << ": unexpected exception: "
                << e.what() << "\n";
      any_failed = true;
    } catch (...) {
      std::cout << "FAIL " << name << ": unknown exception\n";
      any_failed = true;
    }
  }

  return any_failed ? 1 : 0;

}
