#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "scatparams.hpp"

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
