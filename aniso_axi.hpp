// aniso_axi.hpp
//
// Develops AxiProposal, the incidence-tailored proposal distribution
// used by Scatterer::AttemptAnisotropicScatter() (via
// Scatterer::sampler_e::SAMPLER_AXISYMMETRIC) for axisymmetric
// ScatterParams (ax == ay != az). See docs/superpowers/specs/
// 2026-10-16-anisotropic-scale-lengths-design.md, "Axisymmetric
// sampler", for the derivation.
//
// SUMMARY:
//
//   The general sampler (Task 4) draws proposals from the isotropic
//   (a=a_min) take-off-angle tables, which don't know anything about
//   incidence angle; for a strongly axisymmetric medium, most of those
//   proposals get rejected. AxiProposal factors the proposal into (1)
//   the same full-resolution isotropic tables, shared with the general
//   sampler, and (2) a small, coarse, per-incidence-bin weight table
//   that tilts the choice of *cell* (a contiguous block of take-off
//   angles) towards the directions the medium favors at that
//   incidence. A cell is then sampled down to a fine take-off angle by
//   a range-restricted draw from the shared isotropic table (see
//   ProbDist::GetRandomIndexInRange()). The resulting proposal is
//   finished off with an ordinary exact accept/reject test, so
//   sampling stays exact at full angular resolution no matter how
//   coarse the incidence binning is.
//
//   AxiProposal itself is stateless with respect to any particular
//   scattering event: it is built once per axisymmetric Scatterer (at
//   Scatterer construction) from that Scatterer's ScatterParams,
//   take-off-angle set, and isotropic probability tables, and is
//   read-only (all public methods const) thereafter. Everything it
//   stores is expressed in the *local ray frame* (the (S1,S2,S3) basis
//   Scatterer::AttemptAnisotropicScatter() builds for the incoming
//   phonon), never in global or fabric coordinates: the forward
//   (incoming-ray) direction is always (0,0,1) in this frame, and the
//   caller is responsible for projecting the local vertical direction
//   into it (as 'u_rel') before calling BinRate()/DrawProposal().
//
#ifndef ANISO_AXI_H_
#define ANISO_AXI_H_
//
#include <vector>
#include "typedefs.hpp"
#include "raytype.hpp"
#include "geom.hpp"
#include "scatparams.hpp"
#include "probability.hpp"

//////
// CLASS:   ::: AxiProposal :::
// ENCAPS:  Incidence-tailored proposal tables for the axisymmetric
//          sampler.
//
//   Read-only after construction. The 'conv' (conversion type) values
//   used throughout this class and its public interface are plain ints
//   using the same order as Scatterer's (private) out_types_e: GPP=0,
//   GPS=1, GSP=2, GSS=3.
//
class AxiProposal {
public:

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Constants  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  static const int kCoarseDegree = 3; // Cells are blocks of fine take-off
                                      // angles at this TesselSphere
                                      // degree (1,280 of them), or one
                                      // fine angle per cell below this
                                      // degree.
  static const int kRepDegree    = 6; // Bound-construction representative
                                      // directions are capped at this
                                      // TesselSphere degree (81,920 of
                                      // them); the run's own take-off-
                                      // angle set is used directly when
                                      // its degree is <= this.
  static const int kThetaBins    = 45; // Incidence theta bins, 2 degrees
                                       // each, folded to [0,90].
  static const int kChiBins      = 9;  // Incidence chi bins (S only), 10
                                       // degrees each, folded to [0,90].

  static constexpr Real kBoundMargin = 1e-6; // Relative slack applied to
                                             // each stored bound, to
                                             // cover rounding in the
                                             // symmetry folds.
  static constexpr Real kWeightFloor = 1e-9; // Minimum proposal weight
                                             // (relative to a bin's
                                             // total bound*S), so a cell
                                             // with a tiny-but-nonzero
                                             // isotropic weight is never
                                             // given exactly zero chance
                                             // of proposal.

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  AxiProposal(const ScatterParams & par, const S2::S2Set & toa,
              const std::vector<ProbDist> & pdists);

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Query Methods  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  Real MajorantRate(raytype intype) const;
                      // g*_max (1/length): the largest bin majorant
                      // rate for incoming ray type 'intype', over every
                      // incidence bin.

  Real BinRate(raytype intype, const R3::XYZ & u_rel) const;
                      // g*_bin (1/length): the majorant rate of the
                      // incidence bin 'u_rel' (the local vertical
                      // direction, expressed in the incoming ray's
                      // (S1,S2,S3) frame) folds into.

  struct Draw {
    int   Conv;       // out_types_e code (GPP=0,GPS=1,GSP=2,GSS=3) of
                      // the drawn conversion.
    Index Toa;        // Fine take-off-angle index drawn within the
                      // selected (post-mirror) cell.
    Real  Weight;     // The stored proposal weight of the table entry
                      // that was drawn (pre-mirror cell, same bin).
    Real  CellSum;    // S(conv,cell): the isotropic table's sum over
                      // the drawn (post-mirror, actual) cell's take-
                      // off-angle range, for 'Conv'. Freshly computed
                      // from 'pdists', not stored.
  };                  // Acceptance probability = AnisoRatio(m)*CellSum
                      // / Weight (see Scatterer::AttemptAnisotropic-
                      // Scatter()).

  Draw DrawProposal(raytype intype, const R3::XYZ & u_rel) const;
                      // Draws a (conv, cell) pair from u_rel's
                      // incidence bin's coarse cumulative weights, then
                      // a fine take-off-angle index within that cell
                      // from the shared isotropic table.

  Real AcceptanceBound(raytype intype, const R3::XYZ & u_rel,
                       int conv, Index cell) const;
                      // Weight/CellSum for the given (actual) 'conv'
                      // and 'cell', in u_rel's incidence bin: an upper
                      // bound on AnisoRatio(m) for every fine direction
                      // in 'cell' and every incidence state in the bin.
                      // Used by tests to validate bound coverage
                      // directly (DrawProposal() applies the same ratio
                      // internally).

  Index CellOf(Index toa) const;
                      // The cell index owning fine take-off-angle index
                      // 'toa'.

  Index NumCells() const {return mNumCells;}
                      // Number of cells (1,280 at or above
                      // kCoarseDegree; the run's own take-off-angle
                      // count below it).

private:

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Types  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  struct FoldResult {
    int  ThetaBin;
    int  ChiBin;     // Always 0 for intype==RAY_P (canonical frame).
    bool Mx;
    bool My;         // Always false for intype==RAY_P.
  };

  typedef std::vector<double> BinTable;  // Cumulative proposal weights
                                         // over (conv,cell), concatenated
                                         // as [conv0 cells][conv1 cells];
                                         // length 2*NumCells().

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Methods  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  void BuildMirrorPerms(int coarseDeg);

  void BuildBinTable(raytype intype, const R3::XYZ & u_c, Real delta,
                     const std::vector<R3::XYZ> & repCenters,
                     const std::vector<Real> & repRadius,
                     Index repsPerCell,
                     BinTable & table, Real & rate) const;

  Real ComputeBound(const R3::XYZ & t, Real r_i,
                    Real k_in, Real k_out,
                    const R3::XYZ & u_c, Real delta) const;

  FoldResult Fold(raytype intype, R3::XYZ u) const;

  Index MapCell(Index cell, bool mx, bool my) const;
                      // Applies Mx^mx(My^my(cell)); self-inverse (the
                      // mirror permutations are involutions and
                      // commute), so this same map converts a table
                      // cell to an actual cell, or vice versa.

  const BinTable & TableFor(raytype intype, int thetaBin, int chiBin) const;
  Real BinRateFor(raytype intype, int thetaBin, int chiBin) const;

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Member Data  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  const std::vector<ProbDist> & mPDists; // Scatterer's isotropic
                                         // (a=a_min) take-off-angle
                                         // tables, indexed by the same
                                         // out_types_e codes as 'Conv'.

  Index mToaSize;
  Index mNumCells;
  Index mCellBlock;    // Fine take-off-angle indices per cell.

  std::vector<Index> mMxPerm;  // cell -> x-mirrored cell
  std::vector<Index> mMyPerm;  // cell -> y-mirrored cell

  std::vector<BinTable> mPTable;   // [thetaBin]
  std::vector<BinTable> mSTable;   // [thetaBin*kChiBins + chiBin]

  std::vector<Real> mPBinRate;     // [thetaBin]
  std::vector<Real> mSBinRate;     // [thetaBin*kChiBins + chiBin]

  Real mPMajorant;
  Real mSMajorant;

  Real mM;       // ScatterParams::AnisoMajorantFactor()
  Real mAh;      // ax == ay
  Real mAz;
  Real mAMin;    // min(ah,az) == ScatterParams::GetA()
  Real mKappa;
  Real mKP;      // ScatterParams::GetKP()
  Real mKS;      // ScatterParams::GetKS()

}; // class AxiProposal
////

///
#endif //#ifndef ANISO_AXI_H_
//
