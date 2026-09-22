#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "scatparams.hpp"

static R3::XYZ RotateY(const R3::XYZ & v, Real angle) {
  const Real c = std::cos(angle);
  const Real s = std::sin(angle);
  return R3::XYZ(c*v.x() + s*v.z(), v.y(), -s*v.x() + c*v.z());
}

static R3::XYZ RotateX(const R3::XYZ & v, Real angle) {
  const Real c = std::cos(angle);
  const Real s = std::sin(angle);
  return R3::XYZ(v.x(), c*v.y() - s*v.z(), s*v.y() + c*v.z());
}

static void AssertClose(Real lhs, Real rhs) {
  assert(std::abs(lhs - rhs) < 1.0e-11);
}

int main() {
  Elastic::HSneak spectrum(0.8, 0.01, 1.0, 0.5);
  ScatterParams isotropic(spectrum, 1.0, 1.7, 2.0, 2.0);
  ScatterParams anisotropic(spectrum, 1.0, 1.7, 2.0, 1.0);

  assert(!isotropic.IsAnisotropic());
  assert(anisotropic.IsAnisotropic());
  assert(isotropic.GetHorizontalCorrelationLength() == 2.0);
  assert(isotropic.GetVerticalCorrelationLength() == 2.0);

  const R3::XYZ vertical_q(0.0, 0.0, 1.0);
  const R3::XYZ horizontal_q(1.0, 0.0, 0.0);
  const Real iso_vertical = isotropic.PowerSpectralDensity(vertical_q);
  const Real iso_horizontal = isotropic.PowerSpectralDensity(horizontal_q);
  assert(std::abs(iso_vertical - iso_horizontal) < 1.0e-12);

  const Real aniso_vertical = anisotropic.PowerSpectralDensity(vertical_q);
  const Real aniso_horizontal = anisotropic.PowerSpectralDensity(horizontal_q);
  assert(std::abs(aniso_vertical - aniso_horizontal) > 1.0e-8);

  const R3::XYZ incoming(0.6, 0.0, 0.8);
  const R3::XYZ vertical(0.0, 0.0, 1.0);
  const S2::S2Point toa = S2::ThetaPhi(0.9, 1.1);
  Real gpp, gps, gsp, gss, spol;
  anisotropic.GSATO(incoming, vertical, toa, gpp, gps, gsp, gss, spol);

  const Real rotation = 0.7;
  const R3::XYZ rotated_incoming =
      RotateX(RotateY(incoming, rotation), 0.4);
  const R3::XYZ rotated_vertical =
      RotateX(RotateY(vertical, rotation), 0.4);
  Real rotated_gpp, rotated_gps, rotated_gsp, rotated_gss, rotated_spol;
  anisotropic.GSATO(rotated_incoming, rotated_vertical, toa,
                    rotated_gpp, rotated_gps, rotated_gsp, rotated_gss,
                    rotated_spol);
  AssertClose(gpp, rotated_gpp);
  AssertClose(gps, rotated_gps);
  AssertClose(gsp, rotated_gsp);
  AssertClose(gss, rotated_gss);
  AssertClose(spol, rotated_spol);

  bool rejected = false;
  try {
    ScatterParams::SetGlobalCorrelationLengths(1.0, 0.0);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  assert(rejected);
  ScatterParams::SetGlobalCorrelationLengths(0.0, 0.0);

  std::cout << "anisotropic scattering tests passed\n";
  return 0;
}
