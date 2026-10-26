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
              const std::vector<ProbDist> & pdists, unsigned threads = 1);
                      // 'threads': the number of std::thread workers
                      // (ParallelFor(), parallel.hpp) used to build this
                      // AxiProposal: the geometry (representative radii
                      // and mirror-cell permutations, on a cache miss or
                      // fingerprint mismatch -- see
                      // GetOrBuildAxiGeometry(), aniso_axi.cpp), the
                      // per-representative BoundRep arrays, and the
                      // combined P+S bin table (Task 4b). Every result
                      // -- and so every query method of this class -- is
                      // bit-identical regardless of 'threads': each
                      // representative, mirror-cell, or bin is computed
                      // by exactly one task, writing only its own
                      // array/table/rate entry, and mPMajorant/
                      // mSMajorant are reduced afterward, serially, in
                      // ascending bin-index order. Default 1 (serial)
                      // keeps every pre-Task-3 call site unchanged.

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Query Methods  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  struct BinRef {
    int  ThetaBin;
    int  ChiBin;     // Always 0 for intype==RAY_P (canonical frame).
    bool Mx;
    bool My;         // Always false for intype==RAY_P.
  };                  // The incidence bin a local-vertical direction
                      // folds into, plus the two mirror flags needed to
                      // map a table cell back to the actual cell.
                      // Replaces the former private FoldResult: public
                      // so a caller that needs several queries for the
                      // same incidence (BinRate() then DrawProposal(),
                      // e.g. Scatterer::AttemptAxisymmetricScatter())
                      // can locate it once and reuse it, instead of
                      // relocating on every call.

  BinRef Locate(raytype intype, const R3::XYZ & u_rel) const;
                      // Locates 'u_rel' (the local vertical direction,
                      // expressed in the incoming ray's (S1,S2,S3)
                      // frame) into its incidence bin: the same fold
                      // the former private Fold() performed, but via
                      // binary search over precomputed edge
                      // cosines/sines instead of acos/atan2/sin -- see
                      // docs/superpowers/specs/2026-10-17-axisymmetric-
                      // sampler-performance-design.md, decision 4, and
                      // tests/anisotropy_tests.cpp's
                      // ReferenceLocate3db253a (ruling R-P1).

  Real MajorantRate(raytype intype) const;
                      // g*_max (1/length): the largest bin majorant
                      // rate for incoming ray type 'intype', over every
                      // incidence bin.

  Real BinRate(raytype intype, const R3::XYZ & u_rel) const;
                      // g*_bin (1/length): the majorant rate of the
                      // incidence bin 'u_rel' (the local vertical
                      // direction, expressed in the incoming ray's
                      // (S1,S2,S3) frame) folds into. Implemented via
                      // Locate().

  Real BinRate(raytype intype, const BinRef & bin) const;
                      // Same as above, for a bin already located by a
                      // prior Locate() call (avoids relocating).

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
                      // from the shared isotropic table. Implemented
                      // via Locate().

  Draw DrawProposal(raytype intype, const BinRef & bin) const;
                      // Same as above, for a bin already located by a
                      // prior Locate() call (avoids relocating).

  Real AcceptanceBound(raytype intype, const R3::XYZ & u_rel,
                       int conv, Index cell) const;
                      // Weight/CellSum for the given (actual) 'conv'
                      // and 'cell', in u_rel's incidence bin: an upper
                      // bound on AnisoRatio(m) for every fine direction
                      // in 'cell' and every incidence state in the bin.
                      // Used by tests to validate bound coverage
                      // directly (DrawProposal() applies the same ratio
                      // internally). Implemented via Locate().

  Index CellOf(Index toa) const;
                      // The cell index owning fine take-off-angle index
                      // 'toa'.

  Index NumCells() const {return mNumCells;}
                      // Number of cells (1,280 at or above
                      // kCoarseDegree; the run's own take-off-angle
                      // count below it).

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Bound-Kernel Types/Methods  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //   The per-(representative, k_in, k_out) and per-(incidence bin)
  //   halves of the old ComputeBound()'s rigorous AnisoRatio(m) bound
  //   (see docs/superpowers/specs/2026-10-16-anisotropic-scale-lengths-
  //   design.md and docs/superpowers/specs/2026-10-17-axisymmetric-
  //   sampler-performance-design.md, decision 1), split apart so the
  //   (representative, k_in, k_out) part -- everything but the per-bin
  //   cap-edge terms -- is computed once instead of once per
  //   (representative, bin). A cell's bound is
  //   'min(M, M*pow(maxRatio, kappa+1.5))', where maxRatio is the
  //   largest BoundRatio() over the cell's representatives; see
  //   BuildBinTable() in aniso_axi.cpp. Exposed publicly (and as static
  //   methods, since neither depends on instance state) so tests can
  //   exercise the kernel directly against a frozen reference
  //   (tests/anisotropy_tests.cpp's ReferenceBound3db253a, ruling R-P1).

  struct BoundRep {
    R3::XYZ MHat;  // m_rep.Unit(), where m_rep = k_out*t - k_in*zhat;
                   // the zero vector if Mr == 0 (m_rep.Unit() is never
                   // evaluated in that case -- see PrepareBoundRep()).
    Real Mr;       // |m_rep|
    Real Rho;      // k_out * 2*sin(r_i/2): converts representative
                   // 't's angular radius r_i into a bound on how far
                   // any fine direction's exchange wavevector can land
                   // from m_rep.
    Real M2Lo;     // max(Mr-Rho,0)^2 : lower bound on the represented
                   // directions' |m|^2.
    Real M2Hi;     // (Mr+Rho)^2      : upper bound on same.
    Real Numer;    // 1 + aMin^2*M2Hi : the bin-independent numerator of
                   // the AnisoRatio bound.
  };

  struct BinCap {
    R3::XYZ Uc;       // Incidence-bin centre direction.
    Real    Delta;    // Incidence-bin angular half-width (plus the
                      // chi-width contribution, for S bins).
    Real    CosDelta; // cos(Delta); unused until the trig-free cap
                      // kernel (design decision 3) -- computed here so
                      // that later task's BoundRatio() doesn't need a
                      // new MakeBinCap() call at every call site.
    Real    SinDelta; // sin(Delta); ditto.
  };

  static BoundRep PrepareBoundRep(const R3::XYZ & t, Real r_i,
                                  Real k_in, Real k_out, Real aMin);
                      // The (representative, k_in, k_out)-only part of
                      // the old ComputeBound(): m_rep, Mr, Rho, M2Lo,
                      // M2Hi, and the numerator Numer.

  static BinCap MakeBinCap(const R3::XYZ & u_c, Real delta);

  static Real BoundRatio(const BoundRep & rep, const BinCap & cap,
                         Real ah, Real az);
                      // (1 + aMin^2*m2_hi) / (1 + Q_lo): the quantity
                      // raised to (kappa+1.5) in the cell-bound formula
                      // above, evaluated with the trig-free cap-edge
                      // identities of design decision 3 (no acos, cos,
                      // or sin call) -- exactly the arithmetic of the
                      // ComputeBound() this replaces. This is the scalar
                      // reference implementation: no production code
                      // calls it any more (BuildBinTable() below calls
                      // only RatioRun()), it exists solely so
                      // RatioRun() -- and tests/anisotropy_tests.cpp's
                      // test_axi_ratio_run_matches_scalar and
                      // axi_kernel_matches_reference -- have a fixed,
                      // unvectorized target to check against.
                      // RatioRun() must match this function bit-for-bit
                      // (not merely within a tolerance): see its own doc
                      // comment below for why that's achievable despite
                      // vectorization, and why it originally wasn't.

  struct RepSoA {
    std::vector<Real> Mx, My, Mz, Mr, Rho, M2Lo, M2Hi, Numer;
    void Resize(std::size_t n);
    void AssignOne(std::size_t i, const BoundRep & rep);
    void Assign(const std::vector<BoundRep> & reps);
  };                  // Structure-of-arrays BoundRep: one array set per
                      // (k_in,k_out) pair, feeding RatioRun() below so
                      // its loop over representatives reads flat,
                      // contiguous Real arrays (letting GCC vectorize
                      // it) instead of an array of BoundRep structs.
                      // Resize()/AssignOne() let the AxiProposal
                      // constructor fill a RepSoA one representative at
                      // a time, directly from PrepareBoundRep()'s
                      // return value, with no intermediate
                      // vector<BoundRep> array ever allocated (ledger
                      // I2: see the constructor's own comment on why
                      // that matters for peak build memory). Assign()
                      // -- Resize() plus a loop of AssignOne(), kept for
                      // tests/anisotropy_tests.cpp's
                      // test_axi_ratio_run_matches_scalar, which builds
                      // a RepSoA directly from a std::vector<BoundRep>
                      // to exercise RatioRun() against BoundRatio() --
                      // and AssignOne() both copy column-wise with no
                      // arithmetic, so neither can itself introduce any
                      // rounding difference. Exposed publicly, like
                      // BoundRep/BinCap above, purely for that test;
                      // BuildBinTable() is still the only production
                      // reader of a RepSoA, and the AxiProposal
                      // constructor (aniso_axi.cpp) the only production
                      // writer.

  template <bool AhGeAz>
  static void RatioRun(const RepSoA & reps, Index lo, Index hi,
                       const BinCap & cap, Real ah, Real az, Real * ratio);
                      // Fills ratio[i-lo] = BoundRatio(rep_i, cap, ah, az)
                      // for every representative i in the half-open
                      // range [lo,hi) of 'reps' (rep_i reconstructed
                      // from reps' i'th column), using exactly
                      // BoundRatio()'s per-element arithmetic but
                      // written so the loop vectorizes: both cap-edge
                      // values are computed before selecting between
                      // them, std::max()/std::min() are spelled as
                      // ternaries, and the s2_lo branch test uses '&'
                      // rather than '&&' (see aniso_axi.cpp). The
                      // caller's output buffer need only be sized
                      // (hi-lo), not reps' full size -- BuildBinTable()
                      // calls this once per cell, over just that cell's
                      // contiguous representative block, specifically
                      // so its own scratch is O(repsPerCell) rather than
                      // O(numReps): seeing numReps (not a few dozen)
                      // concurrent calls (one per bin, up to 'threads'
                      // at a time, see BuildBinTable()'s own caller)
                      // each want their own buffer, the earlier
                      // (numReps-per-call) version scaled temporary
                      // memory with 'threads' and could exceed this
                      // class's ~32 MB/scatterer construction-memory
                      // budget at high thread counts. 'AhGeAz' selects
                      // BoundRatio()'s ah>=az branch (true) or its
                      // ah<az branch (false) at compile time -- valid
                      // because the caller already knows which holds
                      // for every representative in a single call (it
                      // depends only on 'ah'/'az', not on the
                      // representative), so the branch doesn't need to
                      // be re-tested per element. Defined, and
                      // explicitly instantiated for both 'AhGeAz'
                      // values, only in aniso_axi.cpp -- never in this
                      // header -- so every caller (BuildBinTable() here
                      // and tests/anisotropy_tests.cpp's
                      // test_axi_ratio_run_matches_scalar) links the
                      // one compiled copy, avoiding the cross-
                      // translation-unit rounding divergence documented
                      // in parallel.hpp's own IMPLEMENTATION NOTE.

private:

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Types  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  typedef std::vector<double> BinTable;  // Cumulative proposal weights
                                         // over (conv,cell), concatenated
                                         // as [conv0 cells][conv1 cells];
                                         // length 2*NumCells().

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Methods  (AxiProposal Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  void BuildBinTable(const R3::XYZ & u_c, Real delta,
                     const RepSoA & reps0,
                     const RepSoA & reps1,
                     Index repsPerCell,
                     const std::vector<Real> & S0,
                     const std::vector<Real> & S1,
                     BinTable & table, Real & rate) const;
                      // Builds one incidence bin's cumulative proposal-
                      // weight table (and its majorant rate) for the
                      // two conversions sharing this bin's u_c/delta:
                      // 'reps0'/'reps1' are RepSoA arrays (one entry
                      // per representative, already keyed to the right
                      // k_in/k_out for this bin type -- see the
                      // constructor), and 'S0'/'S1' are the matching
                      // conversions' per-cell isotropic sums (bin-
                      // independent, hoisted out by the constructor).

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

//////
// TEST SUPPORT (tests/anisotropy_tests.cpp only; ledger M3):
//
//   AxiProposal's take-off-angle-set geometry cache
//   (GetOrBuildAxiGeometry(), aniso_axi.cpp) is a function-local static,
//   so a test cannot reach it, or the exception a corrupted build might
//   throw partway through, directly -- and no take-off-angle set
//   reachable through AxiProposal's public constructor actually triggers
//   either of BuildAxiGeometry()'s own throws today: the non-icosahedral
//   guard (AxiProposal's constructor, "toa.size() is not a multiple of
//   mNumCells") runs *before* the geometry cache is ever touched, by
//   design, so it cannot demonstrate a cache-poisoning bug either way;
//   and the mirror-cell "no match" throw depends only on a canonical,
//   always-exactly-symmetric TesselSphere built from 'coarseDeg' alone
//   (never from the caller's actual take-off-angle data), so it cannot
//   be triggered by any take-off-angle set that survives the
//   non-icosahedral guard. This one-shot, self-disarming hook lets a
//   test simulate "BuildAxiGeometry() throws" regardless, right after its
//   Fingerprint field is set but before anything else is -- the same
//   partially-populated state a throw during the real (now unreachable)
//   mirror-match failure would have left, just earlier -- to check that
//   GetOrBuildAxiGeometry() never leaves a partial, fingerprint-matching
//   entry in its cache for a throw to have caused --
//   test_axi_geometry_cache_survives_build_exception
//   (tests/anisotropy_tests.cpp). Costs one atomic load per
//   geometry-cache-miss build (not per representative or cell) when
//   armed, and nothing when not: false by default, and production code
//   never calls the arming function below.
void AxiProposal_TestForceNextGeometryBuildToThrow();

///
#endif //#ifndef ANISO_AXI_H_
//
