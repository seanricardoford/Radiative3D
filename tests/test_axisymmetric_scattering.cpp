#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "axisymmetric_scattering.hpp"

struct ReferenceTotals {
  Real gpp;
  Real gps;
  Real gsp;
  Real gss;

  Real inverse_p() const { return gpp + gps; }
  Real inverse_s() const { return gsp + gss; }
};

static ReferenceTotals IntegrateReference(const ScatterParams & params,
                                          Real mu_in,
                                          int npsi,
                                          int nzeta) {
  const Real dmu = 2.0 / npsi;
  const Real dzeta = 2.0 * Geometry::Pi / nzeta;
  ReferenceTotals result = {0.0, 0.0, 0.0, 0.0};
  const R3::XYZ incoming(
      std::sqrt(std::max(static_cast<Real>(0.0), 1.0 - mu_in*mu_in)),
      0.0, mu_in);
  const R3::XYZ vertical(0.0, 0.0, 1.0);

  for (int imu = 0; imu < npsi; ++imu) {
    const Real cos_psi = -1.0 + (imu + 0.5) * dmu;
    const Real psi = std::acos(cos_psi);
    for (int izeta = 0; izeta < nzeta; ++izeta) {
      const Real zeta = (izeta + 0.5) * dzeta;
      Real gpp, gps, gsp, gss, spol;
      params.GSATO(incoming, vertical, S2::ThetaPhi(psi, zeta),
                   gpp, gps, gsp, gss, spol);
      const Real weight = 1.0 / (npsi * nzeta);
      result.gpp += gpp * weight;
      result.gps += gps * weight;
      result.gsp += gsp * weight;
      result.gss += gss * weight;
    }
  }
  return result;
}

static void AssertRelativeClose(Real actual, Real expected, Real tolerance) {
  assert(std::abs(actual - expected) /
         std::max(static_cast<Real>(1.0e-30), std::abs(expected))
         < tolerance);
}

int main() {
  Elastic::HSneak spectrum(0.8, 0.01, 1.0, 0.5);
  ScatterParams params(spectrum, 1.0, 1.7, 1.25, 0.625);
  AxisymmetricScatteringKernel kernel(params);

  assert(kernel.IncomingBinCount() >= 32);
  assert(kernel.DirectionNodeCount() >= 512);
  assert(kernel.CDFsAreValid(1.0e-12));

  const Real mus[] = {-0.8, -0.2, 0.0, 0.4, 0.9};
  for (Real mu : mus) {
    const ReferenceTotals reference =
        IntegrateReference(params, mu, 96, 192);
    AssertRelativeClose(kernel.GetInverseMeanFreePath(RAY_P, mu),
                        reference.inverse_p(), 5.0e-3);
    AssertRelativeClose(kernel.GetInverseMeanFreePath(RAY_S, mu),
                        reference.inverse_s(), 5.0e-3);
    AssertRelativeClose(kernel.GetConversionWeight(
                            AxisymmetricScatteringKernel::GPP, mu),
                        reference.gpp, 5.0e-3);
    AssertRelativeClose(kernel.GetConversionWeight(
                            AxisymmetricScatteringKernel::GPS, mu),
                        reference.gps, 5.0e-3);
    AssertRelativeClose(kernel.GetConversionWeight(
                            AxisymmetricScatteringKernel::GSP, mu),
                        reference.gsp, 5.0e-3);
    AssertRelativeClose(kernel.GetConversionWeight(
                            AxisymmetricScatteringKernel::GSS, mu),
                        reference.gss, 5.0e-3);
  }

  std::cout << "axisymmetric scattering tests passed\n";
  return 0;
}
