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
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "elastic.hpp"
#include "geom.hpp"
#include "geom_s2.hpp"
#include "scatparams.hpp"
#include "sources.hpp"
#include "phonons.hpp"
#include "scatterers.hpp"


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
