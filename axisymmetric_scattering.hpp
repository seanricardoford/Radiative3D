// axisymmetric_scattering.hpp
//
// Cached axisymmetric anisotropic scattering kernel.  The kernel retains
// three-dimensional outgoing directions while reducing the tabulated
// scattering physics to incoming polar angle, outgoing deflection angle,
// and azimuth about the local vertical symmetry axis.
//
#ifndef AXISYMMETRIC_SCATTERING_H_
#define AXISYMMETRIC_SCATTERING_H_

#include <cstddef>
#include <vector>

#include "probability.hpp"
#include "raytype.hpp"
#include "scatparams.hpp"

class AxisymmetricScatteringKernel {
public:
  enum conversion_e { GPP, GPS, GSP, GSS, NUM_CONVERSIONS };

  struct Sample {
    Real psi;
    Real zeta;
    Real polarization;
    raytype output_type;
  };

  explicit AxisymmetricScatteringKernel(const ScatterParams & params);

  Real GetInverseMeanFreePath(raytype intype, Real mu_in) const;
  Real GetConversionWeight(conversion_e conversion, Real mu_in) const;
  Real GetAverageInverseMeanFreePath(raytype intype) const;
  Real GetAverageConversionWeight(conversion_e conversion) const;
  Real GetAverageDipole(raytype intype) const;
  Sample GetRandomSample(raytype intype, Real mu_in,
                         RandomEngine & rng) const;

  std::size_t IncomingBinCount() const;
  std::size_t DirectionNodeCount() const;
  bool CDFsAreValid(Real tolerance) const;

private:
  struct Bin {
    Real mu;
    Real inverse_mfp[RAY_NBT];
    Real totals[NUM_CONVERSIONS];
    Real dipole[RAY_NBT];
    std::vector<Real> cdf[NUM_CONVERSIONS];
  };

  std::vector<Bin> mBins;
  std::vector<S2::ThetaPhi> mDirections;
  std::vector<Real> mWeights;
  std::vector<Real> mPolarizations;

  void Bracket(Real mu_in, std::size_t & lower, std::size_t & upper,
               Real & fraction) const;
  static Real Interpolate(Real lower, Real upper, Real fraction);
};

#endif //#ifndef AXISYMMETRIC_SCATTERING_H_
