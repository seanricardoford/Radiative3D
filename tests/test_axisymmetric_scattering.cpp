#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "axisymmetric_scattering.hpp"
#include "phonons.hpp"
#include "scatterers.hpp"

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

static R3::XYZ RotateZ(const R3::XYZ & v, Real angle) {
  const Real c = std::cos(angle);
  const Real s = std::sin(angle);
  return R3::XYZ(c*v.x() - s*v.y(), s*v.x() + c*v.y(), v.z());
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

  const Real invariant_mu = 0.13;
  const R3::XYZ canonical_incoming(
      std::sqrt(1.0 - invariant_mu*invariant_mu), 0.0, invariant_mu);
  const R3::XYZ rotated_incoming = RotateZ(canonical_incoming, 0.37);
  assert(std::abs(kernel.GetInverseMeanFreePath(
      RAY_P, invariant_mu)
      - kernel.GetInverseMeanFreePath(RAY_P, rotated_incoming.z()))
      < 1.0e-12);
  for (int conversion = 0;
       conversion < AxisymmetricScatteringKernel::NUM_CONVERSIONS;
       ++conversion) {
    const AxisymmetricScatteringKernel::conversion_e kind =
        static_cast<AxisymmetricScatteringKernel::conversion_e>(conversion);
    assert(std::isfinite(kernel.GetConversionWeight(kind, invariant_mu)));
  }
  const Real position = (invariant_mu + 1.0) * 0.5 * 64.0;
  const Real lower_mu = -1.0 + 2.0 * std::floor(position) / 64.0;
  const Real upper_mu = -1.0 + 2.0 * std::ceil(position) / 64.0;
  const Real interpolated = kernel.GetInverseMeanFreePath(
      RAY_P, invariant_mu);
  const Real lower_value = kernel.GetInverseMeanFreePath(RAY_P, lower_mu);
  const Real upper_value = kernel.GetInverseMeanFreePath(RAY_P, upper_mu);
  assert(interpolated >= std::min(lower_value, upper_value));
  assert(interpolated <= std::max(lower_value, upper_value));

  S2::TesselSphere toa(S2::TESS_ICO, 3);
  PhononSource::Set_TOA_Array(&toa);
  Scatterer scatterer(params);
  const R3::XYZ incoming(0.6, 0.2, std::sqrt(0.6));
  const R3::XYZ vertical(0.3, -0.4, std::sqrt(0.75));
  RandomEngine rng_a(0x1234);
  RandomEngine rng_b(0x1234);
  Phonon rel_a = scatterer.GetRandomScatteredRelativePhonon(
      RAY_S, incoming, vertical, 0.0, rng_a);
  Phonon rel_b = scatterer.GetRandomScatteredRelativePhonon(
      RAY_S, incoming, vertical, Geometry::Pi90, rng_b);
  const S2::ThetaPhi incoming_direction(incoming.Theta(), incoming.Phi());
  Phonon out_a(incoming_direction, RAY_S);
  Phonon out_b(incoming_direction, RAY_S);
  out_a.SetPolarization(0.0);
  out_b.SetPolarization(Geometry::Pi90);
  out_a.Transform(rel_a);
  out_b.Transform(rel_b);
  assert(std::abs(rel_a.GetDirection().Phi()
                 - rel_b.GetDirection().Phi()) > 1.0e-6);
  assert(out_a.GetDirection().XYZ().DistFrom(
             out_b.GetDirection().XYZ()) < 1.0e-11);
  assert(std::abs(out_a.DirectionOfMotion().Dot(
             out_a.GetDirection().XYZ())) < 1.0e-11);
  assert(std::abs(out_b.DirectionOfMotion().Dot(
             out_b.GetDirection().XYZ())) < 1.0e-11);

  Scatterer::OverrideMFP(12.0, 34.0);
  Scatterer override_scatterer(params);
  RandomEngine path_rng(0x5678);
  Real path_sum_p = 0.0;
  Real path_sum_s = 0.0;
  const int path_count = 4096;
  for (int i = 0; i < path_count; ++i) {
    path_sum_p += override_scatterer.GetRandomPathLength(RAY_P, path_rng);
    path_sum_s += override_scatterer.GetRandomPathLength(RAY_S, path_rng);
  }
  assert(path_sum_p / path_count > 10.5);
  assert(path_sum_p / path_count < 13.5);
  assert(path_sum_s / path_count > 30.0);
  assert(path_sum_s / path_count < 38.0);

  Scatterer::SetNoDeflect();
  Scatterer no_deflect_scatterer(params);
  RandomEngine no_deflect_rng(0x9abc);
  Phonon no_deflect = no_deflect_scatterer.GetRandomScatteredRelativePhonon(
      RAY_S, incoming, vertical, Geometry::Pi90, no_deflect_rng);
  assert(no_deflect.GetDirection().Theta() <= 1.0e-6);
  assert(no_deflect.GetRaytype() == RAY_S);

  std::cout << "axisymmetric scattering tests passed\n";
  return 0;
}
