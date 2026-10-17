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
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "elastic.hpp"
#include "geom.hpp"
#include "geom_s2.hpp"
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

// ::: LopNorCrust(ax,ay,az,kappa,eps) :::
//
//   A ScatterParams reference fixture at Lop Nor crustal velocities
//   (Vp 6.2, Vs 3.6 km/s) and 2 Hz, with caller-chosen heterogeneity
//   scale lengths. nu is fixed at 0.8, matching the Lop Nor crust layer
//   (see user_LopNorCylMoho_inc.cpp).
//
static ScatterParams LopNorCrust(Real ax, Real ay, Real az,
                                  Real kappa=0.5, Real eps=0.01) {
  ScatterParams::SetFrequencyHertz(2.0);
  return ScatterParams(Elastic::VpVs(6.2,3.6),
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

  int realCountP = 0;
  for (int i = 0; i < 2000; ++i) {
    Scatterer::AnisoAttempt att =
      scat.AttemptAnisotropicScatter(RAY_P, origin, dirP, pol);
    if (att.IsReal) {
      ++realCountP;
      CHECK(att.Rel.GetDirection().Theta() < 1.e-6);  // nudge_if_singular() keeps
                                                      // theta slightly off zero
      CHECK(att.Rel.GetRaytype() == RAY_P);
      CHECK(att.InPol == pol);
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
  // rate the simulation actually realizes under its own quadrature.
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
