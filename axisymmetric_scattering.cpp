// axisymmetric_scattering.cpp
//
#include <algorithm>
#include <cmath>
#include <limits>

#include "axisymmetric_scattering.hpp"

namespace {

const std::size_t kIncomingBins = 65;
const int kPsiNodes = 48;
const int kZetaNodes = 64;

void BuildGaussLegendre(int count, std::vector<Real> & nodes,
                        std::vector<Real> & weights) {
  nodes.resize(count);
  weights.resize(count);
  const Real tolerance = 1.0e-15;
  const int half = (count + 1) / 2;

  for (int i = 0; i < half; ++i) {
    Real root = std::cos(Geometry::Pi * (i + 0.75) / (count + 0.5));
    Real previous = 0.0;
    while (std::abs(root - previous) > tolerance) {
      Real p0 = 1.0;
      Real p1 = root;
      for (int degree = 2; degree <= count; ++degree) {
        const Real p2 = ((2.0 * degree - 1.0) * root * p1
                         - (degree - 1.0) * p0) / degree;
        p0 = p1;
        p1 = p2;
      }
      const Real derivative = count * (root * p1 - p0)
                            / (root * root - 1.0);
      previous = root;
      root -= p1 / derivative;
    }

    Real p0 = 1.0;
    Real p1 = root;
    for (int degree = 2; degree <= count; ++degree) {
      const Real p2 = ((2.0 * degree - 1.0) * root * p1
                       - (degree - 1.0) * p0) / degree;
      p0 = p1;
      p1 = p2;
    }
    const Real derivative = count * (root * p1 - p0)
                          / (root * root - 1.0);
    const Real weight = 2.0 / ((1.0 - root * root)
                               * derivative * derivative);
    nodes[i] = -root;
    nodes[count - 1 - i] = root;
    weights[i] = weight;
    weights[count - 1 - i] = weight;
  }
}

int InputIndex(raytype intype) {
  return intype == RAY_P ? 0 : 1;
}

bool IsAllowed(raytype intype,
               AxisymmetricScatteringKernel::conversion_e conversion) {
  if (intype == RAY_P) {
    return conversion == AxisymmetricScatteringKernel::GPP
        || conversion == AxisymmetricScatteringKernel::GPS;
  }
  return conversion == AxisymmetricScatteringKernel::GSP
      || conversion == AxisymmetricScatteringKernel::GSS;
}

} // namespace

AxisymmetricScatteringKernel::AxisymmetricScatteringKernel(
    const ScatterParams & params) {
  std::vector<Real> psi_nodes;
  std::vector<Real> psi_weights;
  BuildGaussLegendre(kPsiNodes, psi_nodes, psi_weights);

  mDirections.reserve(kPsiNodes * kZetaNodes);
  mWeights.reserve(kPsiNodes * kZetaNodes);
  mPolarizations.reserve(kPsiNodes * kZetaNodes);
  const R3::XYZ canonical_incoming(0.0, 0.0, 1.0);
  for (int ipsi = 0; ipsi < kPsiNodes; ++ipsi) {
    const Real psi = std::acos(psi_nodes[ipsi]);
    for (int izeta = 0; izeta < kZetaNodes; ++izeta) {
      const Real zeta = 2.0 * Geometry::Pi * izeta / kZetaNodes;
      Real gpp, gps, gsp, gss, spol;
      params.GSATO(canonical_incoming, R3::XYZ(0,0,1),
                   S2::ThetaPhi(psi, zeta),
                   gpp, gps, gsp, gss, spol);
      // The values stored in this shared node table are re-evaluated below
      // for each incoming bin.  Keep the node and S/S polarization here.
      (void)gpp;
      (void)gps;
      (void)gsp;
      (void)gss;
      mDirections.push_back(S2::ThetaPhi(psi, zeta));
      mWeights.push_back(0.5 * psi_weights[ipsi] / kZetaNodes);
      mPolarizations.push_back(spol);
    }
  }

  mBins.resize(kIncomingBins);
  for (std::size_t ibin = 0; ibin < kIncomingBins; ++ibin) {
    Bin & bin = mBins[ibin];
    bin.mu = -1.0 + 2.0 * ibin / (kIncomingBins - 1.0);
    for (int intype = 0; intype < RAY_NBT; ++intype) {
      bin.inverse_mfp[intype] = 0.0;
      bin.dipole[intype] = 0.0;
    }
    for (int conversion = 0; conversion < NUM_CONVERSIONS; ++conversion) {
      bin.totals[conversion] = 0.0;
      bin.cdf[conversion].resize(mDirections.size(), 0.0);
    }

    const Real transverse = std::sqrt(
        std::max(static_cast<Real>(0.0), 1.0 - bin.mu * bin.mu));
    const R3::XYZ incoming(transverse, 0.0, bin.mu);
    const R3::XYZ vertical(0.0, 0.0, 1.0);
    Real moments[NUM_CONVERSIONS] = {0.0, 0.0, 0.0, 0.0};

    for (std::size_t node = 0; node < mDirections.size(); ++node) {
      Real gpp, gps, gsp, gss, spol;
      params.GSATO(incoming, vertical, mDirections[node],
                   gpp, gps, gsp, gss, spol);
      const Real values[NUM_CONVERSIONS] = {gpp, gps, gsp, gss};
      for (int conversion = 0; conversion < NUM_CONVERSIONS; ++conversion) {
        const Real weighted = values[conversion] * mWeights[node];
        bin.totals[conversion] += weighted;
        moments[conversion] += weighted * mDirections[node].z();
        if (node == 0) {
          bin.cdf[conversion][node] = weighted;
        } else {
          bin.cdf[conversion][node] =
              bin.cdf[conversion][node - 1] + weighted;
        }
      }
    }

    bin.inverse_mfp[0] = bin.totals[GPP] + bin.totals[GPS];
    bin.inverse_mfp[1] = bin.totals[GSP] + bin.totals[GSS];
    for (int intype = 0; intype < RAY_NBT; ++intype) {
      const int first = intype == 0 ? GPP : GSP;
      const int second = intype == 0 ? GPS : GSS;
      const Real total = bin.totals[first] + bin.totals[second];
      bin.dipole[intype] = total > 0.0
          ? (moments[first] + moments[second]) / total : 0.0;
    }
  }
}

Real AxisymmetricScatteringKernel::Interpolate(Real lower, Real upper,
                                                Real fraction) {
  return lower + fraction * (upper - lower);
}

void AxisymmetricScatteringKernel::Bracket(Real mu_in,
                                            std::size_t & lower,
                                            std::size_t & upper,
                                            Real & fraction) const {
  const Real mu = std::max(static_cast<Real>(-1.0),
                           std::min(static_cast<Real>(1.0), mu_in));
  const Real position = (mu + 1.0) * 0.5 * (mBins.size() - 1.0);
  lower = static_cast<std::size_t>(std::floor(position));
  upper = lower + 1 < mBins.size() ? lower + 1 : lower;
  fraction = upper == lower ? 0.0 : position - lower;
}

Real AxisymmetricScatteringKernel::GetInverseMeanFreePath(
    raytype intype, Real mu_in) const {
  std::size_t lower, upper;
  Real fraction;
  Bracket(mu_in, lower, upper, fraction);
  const int index = InputIndex(intype);
  return Interpolate(mBins[lower].inverse_mfp[index],
                     mBins[upper].inverse_mfp[index], fraction);
}

Real AxisymmetricScatteringKernel::GetConversionWeight(
    conversion_e conversion, Real mu_in) const {
  std::size_t lower, upper;
  Real fraction;
  Bracket(mu_in, lower, upper, fraction);
  return Interpolate(mBins[lower].totals[conversion],
                     mBins[upper].totals[conversion], fraction);
}

Real AxisymmetricScatteringKernel::GetAverageInverseMeanFreePath(
    raytype intype) const {
  Real result = 0.0;
  const int index = InputIndex(intype);
  for (std::size_t i = 0; i < mBins.size(); ++i) {
    const Real weight = (i == 0 || i + 1 == mBins.size()) ? 0.5 : 1.0;
    result += weight * mBins[i].inverse_mfp[index];
  }
  return result / (mBins.size() - 1.0);
}

Real AxisymmetricScatteringKernel::GetAverageConversionWeight(
    conversion_e conversion) const {
  Real result = 0.0;
  for (std::size_t i = 0; i < mBins.size(); ++i) {
    const Real weight = (i == 0 || i + 1 == mBins.size()) ? 0.5 : 1.0;
    result += weight * mBins[i].totals[conversion];
  }
  return result / (mBins.size() - 1.0);
}

Real AxisymmetricScatteringKernel::GetAverageDipole(raytype intype) const {
  Real result = 0.0;
  const int index = InputIndex(intype);
  for (std::size_t i = 0; i < mBins.size(); ++i) {
    const Real weight = (i == 0 || i + 1 == mBins.size()) ? 0.5 : 1.0;
    result += weight * mBins[i].dipole[index];
  }
  return result / (mBins.size() - 1.0);
}

AxisymmetricScatteringKernel::Sample
AxisymmetricScatteringKernel::GetRandomSample(raytype intype, Real mu_in,
                                              RandomEngine & rng) const {
  std::size_t lower, upper;
  Real fraction;
  Bracket(mu_in, lower, upper, fraction);

  Real conversion_weights[NUM_CONVERSIONS] = {0.0, 0.0, 0.0, 0.0};
  Real total = 0.0;
  for (int conversion = 0; conversion < NUM_CONVERSIONS; ++conversion) {
    if (IsAllowed(intype,
                  static_cast<conversion_e>(conversion))) {
      conversion_weights[conversion] =
          Interpolate(mBins[lower].totals[conversion],
                      mBins[upper].totals[conversion], fraction);
      total += conversion_weights[conversion];
    }
  }

  Sample fallback = {0.0, 0.0, 0.0, intype == RAY_P ? RAY_P : RAY_S};
  if (!(total > 0.0) || !std::isfinite(total)) {
    return fallback;
  }

  const Real conversion_pick = rng.Uniform01() * total;
  Real cumulative = 0.0;
  conversion_e selected = GSS;
  for (int conversion = 0; conversion < NUM_CONVERSIONS; ++conversion) {
    cumulative += conversion_weights[conversion];
    if (conversion_pick <= cumulative) {
      selected = static_cast<conversion_e>(conversion);
      break;
    }
  }

  const Real lower_weight = (1.0 - fraction) * mBins[lower].totals[selected];
  const Real upper_weight = fraction * mBins[upper].totals[selected];
  const Real bin_pick = rng.Uniform01() * (lower_weight + upper_weight);
  const Bin & bin = bin_pick <= lower_weight
      ? mBins[lower] : mBins[upper];
  const std::vector<Real> & cdf = bin.cdf[selected];
  if (cdf.empty() || !(cdf.back() > 0.0)) {
    return fallback;
  }

  const Real direction_pick = rng.Uniform01() * cdf.back();
  const std::vector<Real>::const_iterator iter =
      std::lower_bound(cdf.begin(), cdf.end(), direction_pick);
  const std::size_t index = iter == cdf.end()
      ? cdf.size() - 1 : static_cast<std::size_t>(iter - cdf.begin());
  Sample result;
  result.psi = mDirections[index].Theta();
  result.zeta = mDirections[index].Phi();
  result.polarization = selected == GSS ? mPolarizations[index] : 0.0;
  result.output_type = (selected == GPP || selected == GSP)
      ? RAY_P : RAY_S;
  return result;
}

std::size_t AxisymmetricScatteringKernel::IncomingBinCount() const {
  return mBins.size();
}

std::size_t AxisymmetricScatteringKernel::DirectionNodeCount() const {
  return mDirections.size();
}

bool AxisymmetricScatteringKernel::CDFsAreValid(Real tolerance) const {
  for (const Bin & bin : mBins) {
    for (int conversion = 0; conversion < NUM_CONVERSIONS; ++conversion) {
      const std::vector<Real> & cdf = bin.cdf[conversion];
      if (cdf.empty() || !std::isfinite(cdf.back())) return false;
      if (std::abs(cdf.back() - bin.totals[conversion])
          > tolerance * std::max(static_cast<Real>(1.0),
                                 std::abs(bin.totals[conversion]))) {
        return false;
      }
      for (std::size_t i = 1; i < cdf.size(); ++i) {
        if (!std::isfinite(cdf[i]) || cdf[i] + tolerance < cdf[i-1]) {
          return false;
        }
      }
    }
    for (int intype = 0; intype < RAY_NBT; ++intype) {
      const int first = intype == 0 ? GPP : GSP;
      const int second = intype == 0 ? GPS : GSS;
      if (!(bin.totals[first] > 0.0) || !(bin.totals[second] > 0.0)) {
        return false;
      }
      if (!std::isfinite(bin.inverse_mfp[intype])
          || !std::isfinite(bin.dipole[intype])) {
        return false;
      }
    }
  }
  return true;
}
