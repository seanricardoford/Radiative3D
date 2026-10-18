// aniso_axi.cpp
//
#include <algorithm>    /* std::min, std::max */
#include <cmath>        /* sin, cos, acos, atan2, sqrt, pow, floor, log */
#include <map>          /* std::map: geometry cache, keyed by TOA count */
#include <mutex>        /* std::mutex, std::lock_guard: geometry cache */
#include <stdexcept>    /* std::logic_error */
#include "aniso_axi.hpp"
#include "random.hpp"

using namespace std;

// Out-of-class definitions for the static constexpr members, needed
// (C++11) whenever they are odr-used (e.g. bound to a const Real& in a
// std::min/std::max call). See scatterers.cpp's use of
// AxiProposal::kBoundMargin.
constexpr Real AxiProposal::kBoundMargin;
constexpr Real AxiProposal::kWeightFloor;


//////
// LOCAL HELPERS:
//

namespace {

  // Global out_types_e-equivalent conversion codes, duplicated here
  // because Scatterer::out_types_e is a private implementation detail
  // of Scatterer: GPP=0, GPS=1, GSP=2, GSS=3.
  const raytype kConvOutRT[4] = {RAY_P, RAY_S, RAY_P, RAY_S};

  inline Index Pow4(int n) {
    Index p = 1;
    for (int i = 0; i < n; i++) {p *= 4;}
    return p;
  }

  // ::: ComputeAxiShape(...) :::
  //
  //   The cell/block/representative sizing arithmetic that used to live
  //   directly in AxiProposal's constructor: a pure function of the
  //   take-off-angle count alone. Factored out into the single place
  //   both the constructor (which needs numCells/cellBlock up front, to
  //   run the non-icosahedral guard *before* touching the geometry
  //   cache below) and GetOrBuildAxiGeometry() (which needs the rest)
  //   call -- neither keeps its own copy of the formula, so the two
  //   genuinely cannot drift apart (not merely "shouldn't"). Every
  //   formula here is unchanged from the pre-Task-2 constructor.
  //
  void ComputeAxiShape(Index toaSize,
                       Index & numCells, Index & cellBlock,
                       int & coarseDeg, int & repDeg,
                       Index & numReps, Index & repBlock,
                       Index & repsPerCell, bool & repEqualsToa) {

    // D = log4(toa.size()/20):
    const int D = (int) ( log( (double) toaSize / 20.0 ) / log(4.0) + 0.5 );

    coarseDeg = std::min(D, (int) AxiProposal::kCoarseDegree);
    numCells  = 20 * Pow4(coarseDeg);
    cellBlock = toaSize / numCells;

    repDeg = std::min(D, (int) AxiProposal::kRepDegree);
    repEqualsToa = (repDeg == D);
    if (repEqualsToa) {
      numReps  = toaSize;
      repBlock = 1;
    } else {
      numReps  = 20 * Pow4(repDeg);
      repBlock = toaSize / numReps;
    }
    repsPerCell = numReps / numCells;
  }

  // ::: FingerprintOf(toa) :::
  //
  //   A cheap scalar fingerprint of a take-off-angle set: controller
  //   ruling R2's index-weighted sum
  //   Sum_k (k+1)*(x_k + 2*y_k + 3*z_k), k in index order, 0-based. This
  //   is only a sanity check guarding the geometry cache below against
  //   reuse across two distinct take-off-angle sets that happen to share
  //   a size (not a security hash) -- production builds a single
  //   TesselSphere per run and reuses it for every Scatterer (see
  //   PhononSource::Set_TOA_Array()), so a collision can't arise there.
  //
  //   An earlier version of this fingerprint was a plain (unweighted,
  //   equally-coefficiented) component sum, which a standalone check
  //   showed to be a poor discriminator: TESS_ICO's vertex directions
  //   sum to (numerically) zero by icosahedral symmetry, and since
  //   rotation is linear, a rotated copy's sum is *also* (numerically)
  //   zero -- two such fingerprints could differ by as little as
  //   ~1e-12 relative, no better than rounding noise. Per-index
  //   weighting breaks that cancellation (a rotation or reordering
  //   changes which direction sits at which index, and a different
  //   index carries a different weight) and the x/2y/3z coefficients
  //   break the remaining x<->y<->z degeneracy; a standalone check
  //   confirmed a 0.01 rad rotation about z of the degree-5 TESS_ICO
  //   set now changes this fingerprint by about 0.7% (relative) --
  //   orders of magnitude above floating-point noise -- and that
  //   reversing the set's order changes it by more than 100%. See
  //   test_axi_geometry_cache_mismatch and
  //   test_axi_geometry_cache_mismatch_rotation (tests/anisotropy_tests.cpp)
  //   for both as regression-tested mismatch cases.
  //
  Real FingerprintOf(const S2::S2Set & toa) {
    Real sum = 0.0;
    for (Index k = 0; k < (Index) toa.size(); k++) {
      const R3::XYZ p = toa[k].XYZ();
      sum += (Real) (k+1) * (p.x() + 2.0*p.y() + 3.0*p.z());
    }
    return sum;
  }

  // ::: struct CachedAxiGeometry ::: / ::: BuildAxiGeometry(...) ::: /
  // ::: GetOrBuildAxiGeometry(toa, scratch) :::
  //
  //   The take-off-angle-set-dependent, ScatterParams-independent part
  //   of AxiProposal's construction: representative centres/radii
  //   (formerly computed inline in the constructor) and the mirror-cell
  //   permutations (formerly AxiProposal::BuildMirrorPerms()). This is
  //   identical for every Scatterer built at a given take-off-angle
  //   degree -- see docs/superpowers/specs/2026-10-17-axisymmetric-
  //   sampler-performance-design.md, decision 1 -- so it is cached in a
  //   function-local static map keyed by take-off-angle count, guarded
  //   by a mutex. Construction of Scatterers (and so of AxiProposal) is
  //   single-threaded today; the mutex is defensive, and cheap enough
  //   (one lock per AxiProposal, not per representative or bin) not to
  //   matter either way. A fingerprint (FingerprintOf()) guards against
  //   serving a cached entry to a different take-off-angle set that
  //   happens to share a size: per the brief, a mismatch rebuilds the
  //   geometry for that one construction *without caching it* -- the
  //   result goes into the caller-owned 'scratch' instead of the map, so
  //   the entry already cached under that size (which may still be
  //   exactly right for whatever call eventually does come back with the
  //   set it was built from) is never evicted or overwritten by a
  //   same-size-but-different set passing through.
  //
  struct CachedAxiGeometry {
    Real                  Fingerprint;
    Index                 RepsPerCell;
    std::vector<R3::XYZ>  RepCenters;
    std::vector<Real>     RepRadius;
    std::vector<Index>    MxPerm;
    std::vector<Index>    MyPerm;
  };

  // ::: BuildAxiGeometry(toa, toaSize, fp, geom) :::
  //
  //   Does the actual (expensive) geometry construction, writing into
  //   'geom' (which the caller owns -- either a fresh map slot on a true
  //   cache miss, or a construction-local scratch object on a
  //   fingerprint mismatch; see GetOrBuildAxiGeometry() below). Split
  //   out from GetOrBuildAxiGeometry() so both call sites share one
  //   body instead of risking two copies drifting apart.
  //
  void BuildAxiGeometry(const S2::S2Set & toa, Index toaSize, Real fp,
                        CachedAxiGeometry & geom) {

    Index numCells, cellBlock, numReps, repBlock, repsPerCell;
    int coarseDeg, repDeg; bool repEqualsToa;
    ComputeAxiShape(toaSize, numCells, cellBlock, coarseDeg, repDeg,
                   numReps, repBlock, repsPerCell, repEqualsToa);

    geom.Fingerprint = fp;
    geom.RepsPerCell = repsPerCell;

    // Representative directions (for bound construction only):
    geom.RepCenters.resize(numReps);
    geom.RepRadius.assign(numReps, 0.0);
    if (repEqualsToa) {
      for (Index i = 0; i < numReps; i++) {
        geom.RepCenters[i] = toa[i].XYZ();
      }
      // RepRadius[] stays all-zero: each representative is its own
      // (and only) fine direction.
    } else {
      S2::TesselSphere repTess(S2::TESS_ICO, repDeg);
      for (Index i = 0; i < numReps; i++) {
        geom.RepCenters[i] = repTess[i].XYZ();
      }
      for (Index i = 0; i < numReps; i++) {
        Real maxAngle = 0.0;
        const Index base = i * repBlock;
        for (Index k = base; k < base + repBlock; k++) {
          Real c = geom.RepCenters[i].Dot(toa[k].XYZ());
          if (c > 1.0) c = 1.0;
          if (c < -1.0) c = -1.0;
          const Real ang = acos(c);
          if (ang > maxAngle) {maxAngle = ang;}
        }
        geom.RepRadius[i] = maxAngle;
      }
    }

    // Mirror-cell permutations (formerly AxiProposal::BuildMirrorPerms()):
    // for each coarse cell, the index of the cell whose centre is the x-
    // (resp. y-) mirror image of this cell's centre, matched within
    // 1e-9. Throws std::logic_error if any cell has no match (the
    // icosahedral tessellation used throughout this class is, by
    // construction, symmetric under negating x or y alone -- see
    // geom_s2.cpp's init_icosahedron() -- so a missing match signals a
    // real inconsistency, not just tessellation noise). Built from a
    // fresh canonical TesselSphere, not from 'toa' itself, so this part
    // never depends on 'toa's actual directions (only indirectly, via
    // 'coarseDeg', on its size).
    //
    // NOTE: building a fresh degree-coarseDeg TesselSphere always gives
    // exactly the same direction set as grouping 'toa' into numCells
    // cells would (whether coarseDeg==D or coarseDeg<D), because
    // TesselSphere's recursive quadrisection preserves each cell's
    // three corner nodes at every finer degree: see geom_s2.cpp's
    // Triangle::subdivide().
    std::vector<R3::XYZ> centers(numCells);
    S2::TesselSphere cellTess(S2::TESS_ICO, coarseDeg);
    for (Index j = 0; j < numCells; j++) {
      centers[j] = cellTess[j].XYZ();
    }

    geom.MxPerm.resize(numCells);
    geom.MyPerm.resize(numCells);

    const Real tol2 = 1e-9*1e-9;

    for (Index j = 0; j < numCells; j++) {

      const R3::XYZ mx_target(-centers[j].x(), centers[j].y(), centers[j].z());
      const R3::XYZ my_target( centers[j].x(),-centers[j].y(), centers[j].z());

      Index bestMx = 0, bestMy = 0;
      Real bestMxD2 = 1.0e300, bestMyD2 = 1.0e300;

      for (Index k = 0; k < numCells; k++) {
        const Real dMx = centers[k].DistFrom(mx_target);
        const Real dMy = centers[k].DistFrom(my_target);
        if (dMx*dMx < bestMxD2) {bestMxD2 = dMx*dMx; bestMx = k;}
        if (dMy*dMy < bestMyD2) {bestMyD2 = dMy*dMy; bestMy = k;}
      }

      if (bestMxD2 > tol2) {
        throw std::logic_error(
          "AxiProposal::BuildMirrorPerms: no x-mirror match for cell");
      }
      if (bestMyD2 > tol2) {
        throw std::logic_error(
          "AxiProposal::BuildMirrorPerms: no y-mirror match for cell");
      }

      geom.MxPerm[j] = bestMx;
      geom.MyPerm[j] = bestMy;
    }

  }//

  // ::: GetOrBuildAxiGeometry(toa, scratch) :::
  //
  //   Returns the geometry for 'toa': a cache hit returns a reference
  //   into the cache; a true miss (no entry yet for this size) builds
  //   directly into a fresh cache slot and returns a reference to it; a
  //   fingerprint mismatch (an entry exists for this size, but isn't
  //   'toa') builds into the caller-owned 'scratch' instead and returns
  //   a reference to *that*, leaving the cache untouched. 'scratch' must
  //   outlive the returned reference in every case -- the constructor
  //   keeps it as a local for exactly as long as it uses the result.
  //
  const CachedAxiGeometry & GetOrBuildAxiGeometry(const S2::S2Set & toa,
                                                  CachedAxiGeometry & scratch) {

    static std::mutex mtx;
    static std::map<Index, CachedAxiGeometry> cache;

    const Index toaSize = (Index) toa.size();
    const Real  fp       = FingerprintOf(toa);

    std::lock_guard<std::mutex> lock(mtx);

    std::map<Index, CachedAxiGeometry>::iterator it = cache.find(toaSize);
    if (it != cache.end()) {
      if (it->second.Fingerprint == fp) {
        return it->second;
      }
      BuildAxiGeometry(toa, toaSize, fp, scratch);
      return scratch;
    }

    CachedAxiGeometry & slot = cache[toaSize];
    BuildAxiGeometry(toa, toaSize, fp, slot);
    return slot;
  }

}//end anonymous namespace


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  AxiProposal                                         ****
// ****                                                              ****
//

//////
// CONSTRUCTOR:  AxiProposal()
//
AxiProposal::AxiProposal(const ScatterParams & par, const S2::S2Set & toa,
                         const std::vector<ProbDist> & pdists) :
  mPDists(pdists)
{

  mM     = par.AnisoMajorantFactor();
  mAh    = par.GetAx();   // == par.GetAy(), by axisymmetric assumption
  mAz    = par.GetAz();
  mAMin  = par.GetA();
  mKappa = par.GetKappa();
  mKP    = par.GetKP();
  mKS    = par.GetKS();

  mToaSize = (Index) toa.size();

  // Shape (cell/block/representative counts): obtained from
  // ComputeAxiShape() -- the one place this arithmetic lives -- rather
  // than re-derived here, so this and GetOrBuildAxiGeometry()'s own call
  // to the same function can never compute a different mNumCells/
  // mCellBlock for the same take-off-angle count.
  Index numCells, cellBlock, shapeNumReps, shapeRepBlock, shapeRepsPerCell;
  int coarseDeg, repDeg; bool repEqualsToa;
  ComputeAxiShape(mToaSize, numCells, cellBlock, coarseDeg, repDeg,
                 shapeNumReps, shapeRepBlock, shapeRepsPerCell, repEqualsToa);

  mNumCells  = numCells;
  mCellBlock = cellBlock;

  if (mToaSize != mNumCells * mCellBlock) {
    // The cell/block arithmetic throughout this class (hierarchical
    // TesselSphere ordering, contiguous per-cell blocks of fine
    // directions) assumes toa.size() is exactly 20*4^D for some D, as
    // TESS_ICO always produces. A non-icosahedral take-off-angle set
    // would silently drop the tail directions from every cell's range,
    // which is a silent bias rather than a loud failure. Production
    // always uses TESS_ICO; this is a belt-and-braces guard. This check
    // runs before the geometry cache below is ever touched, so a bad
    // take-off-angle set always throws fresh -- it is never masked by,
    // or allowed to poison, a cached entry.
    throw std::logic_error(
      "AxiProposal requires an icosahedral (TESS_ICO) take-off-angle "
      "set: toa.size() is not a multiple of mNumCells");
  }

  // Geometry (representative centres/radii, mirror permutations):
  // identical for every Scatterer at this take-off-angle count, so
  // fetched from (and, on first use, built into) a process-wide cache.
  // 'scratch' backs the returned reference only on a fingerprint
  // mismatch (a different, same-size take-off-angle set already
  // occupies this size's cache slot); it must stay alive for as long as
  // 'geom' is used below, which it does as a constructor-local. See
  // GetOrBuildAxiGeometry() above.
  CachedAxiGeometry scratch;
  const CachedAxiGeometry & geom = GetOrBuildAxiGeometry(toa, scratch);
  mMxPerm = geom.MxPerm;
  mMyPerm = geom.MyPerm;

  // Per-representative bound-kernel data (BoundRep), precomputed once
  // for each of the four (k_in,k_out) pairs the P and S bin tables need
  // -- P bins use (kP,kP) and (kP,kS); S bins use (kS,kP) and (kS,kS) --
  // instead of once per (representative, bin) as the old ComputeBound()
  // did. See PrepareBoundRep()/BoundRatio() and BuildBinTable() below.
  const Index numReps = (Index) geom.RepCenters.size();
  std::vector<BoundRep> repsPP(numReps), repsPS(numReps),
                        repsSP(numReps), repsSS(numReps);
  for (Index i = 0; i < numReps; i++) {
    const R3::XYZ & t = geom.RepCenters[i];
    const Real      r = geom.RepRadius[i];
    repsPP[i] = PrepareBoundRep(t, r, mKP, mKP, mAMin);
    repsPS[i] = PrepareBoundRep(t, r, mKP, mKS, mAMin);
    repsSP[i] = PrepareBoundRep(t, r, mKS, mKP, mAMin);
    repsSS[i] = PrepareBoundRep(t, r, mKS, mKS, mAMin);
  }

  // Per-cell isotropic sums S(conv,cell), hoisted out of the bin loop:
  // they depend only on conv and cell, never on the incidence bin, so
  // the old BuildBinTable() recomputed the same mNumCells values on
  // every one of the 45 P (or 405 S) calls that shared a conv group.
  std::vector<Real> S0_P(mNumCells), S1_P(mNumCells);
  std::vector<Real> S0_S(mNumCells), S1_S(mNumCells);
  for (Index c = 0; c < mNumCells; c++) {
    const Index lo = c*mCellBlock, hi = lo+mCellBlock;
    S0_P[c] = mPDists[0].GetRangeMagnitude(lo, hi);  // GPP
    S1_P[c] = mPDists[1].GetRangeMagnitude(lo, hi);  // GPS
    S0_S[c] = mPDists[2].GetRangeMagnitude(lo, hi);  // GSP
    S1_S[c] = mPDists[3].GetRangeMagnitude(lo, hi);  // GSS
  }

  // P incidence bins: theta only (canonical frame fixes chi == 0).
  mPTable.resize(kThetaBins);
  mPBinRate.resize(kThetaBins);
  for (int b = 0; b < kThetaBins; b++) {
    const Real theta_lo = (2.0*b)     * Geometry::DtoR;
    const Real theta_hi = (2.0*b+2.0) * Geometry::DtoR;
    const Real theta_c  = 0.5*(theta_lo+theta_hi);
    const Real delta    = 0.5*(theta_hi-theta_lo);
    const R3::XYZ u_c(-sin(theta_c), 0.0, cos(theta_c));
    BuildBinTable(u_c, delta, repsPP, repsPS, geom.RepsPerCell,
                 S0_P, S1_P, mPTable[b], mPBinRate[b]);
  }
  mPMajorant = 0.0;
  for (int b = 0; b < kThetaBins; b++) {
    mPMajorant = std::max(mPMajorant, mPBinRate[b]);
  }

  // S incidence bins: theta and chi.
  mSTable.resize(kThetaBins*kChiBins);
  mSBinRate.resize(kThetaBins*kChiBins);
  for (int b = 0; b < kThetaBins; b++) {
    const Real theta_lo = (2.0*b)     * Geometry::DtoR;
    const Real theta_hi = (2.0*b+2.0) * Geometry::DtoR;
    const Real theta_c  = 0.5*(theta_lo+theta_hi);
    for (int c = 0; c < kChiBins; c++) {
      const Real chi_lo = (10.0*c)      * Geometry::DtoR;
      const Real chi_hi = (10.0*c+10.0) * Geometry::DtoR;
      const Real chi_c  = 0.5*(chi_lo+chi_hi);
      const Real delta = 0.5*(theta_hi-theta_lo)
                        + sin(theta_hi)*0.5*(chi_hi-chi_lo);
      const R3::XYZ u_c(-sin(theta_c)*cos(chi_c),
                        sin(theta_c)*sin(chi_c),
                        cos(theta_c));
      const int idx = b*kChiBins + c;
      BuildBinTable(u_c, delta, repsSP, repsSS, geom.RepsPerCell,
                   S0_S, S1_S, mSTable[idx], mSBinRate[idx]);
    }
  }
  mSMajorant = 0.0;
  for (std::size_t i = 0; i < mSBinRate.size(); i++) {
    mSMajorant = std::max(mSMajorant, mSBinRate[i]);
  }

}//
//


//////
// METHOD:  AxiProposal :: PrepareBoundRep()                   [static]
//
//   The (representative, k_in, k_out)-only part of the bound kernel
//   that used to be recomputed, per bin, inside ComputeBound(): m_rep's
//   direction and magnitude, the chord-length radius rho, the |m|^2
//   bracket (M2Lo, M2Hi), and the bin-independent numerator Numer. See
//   docs/superpowers/specs/2026-10-16-anisotropic-scale-lengths-
//   design.md and docs/superpowers/specs/2026-10-17-axisymmetric-
//   sampler-performance-design.md, decision 1.
//
AxiProposal::BoundRep
AxiProposal::PrepareBoundRep(const R3::XYZ & t, Real r_i,
                             Real k_in, Real k_out, Real aMin) {

  const R3::XYZ zhat(0.0, 0.0, 1.0);   // Forward (incoming-ray) direction,
                                       // in the local ray frame.
  const R3::XYZ m_rep = t.ScaledBy(k_out) + zhat.ScaledBy(-k_in);

  BoundRep rep;
  rep.Mr   = m_rep.Mag();
  rep.MHat = (rep.Mr > 0.0) ? m_rep.Unit() : R3::XYZ(0.0, 0.0, 0.0);
  rep.Rho  = k_out * 2.0 * sin(r_i/2.0);

  const Real mlo = std::max(rep.Mr-rep.Rho, (Real)0.0);
  rep.M2Lo = mlo*mlo;
  const Real mhi = rep.Mr+rep.Rho;
  rep.M2Hi = mhi*mhi;

  rep.Numer = 1.0 + aMin*aMin*rep.M2Hi;

  return rep;

}//
//


//////
// METHOD:  AxiProposal :: MakeBinCap()                        [static]
//
AxiProposal::BinCap AxiProposal::MakeBinCap(const R3::XYZ & u_c, Real delta) {

  BinCap cap;
  cap.Uc       = u_c;
  cap.Delta    = delta;
  cap.CosDelta = cos(delta);
  cap.SinDelta = sin(delta);
  return cap;

}//
//


//////
// METHOD:  AxiProposal :: BoundRatio()                        [static]
//
//   Rigorous upper bound on AnisoRatio(m) (see ScatterParams::
//   AnisoRatio()), divided by M, over every fine direction within angle
//   'rep's originating r_i of its representative direction, and every
//   local-vertical direction within angle 'cap.Delta' of 'cap.Uc'. This
//   is exactly today's ComputeBound() (frozen, for ruling R-P1, as
//   tests/anisotropy_tests.cpp's ReferenceBound3db253a), split into the
//   per-representative part (PrepareBoundRep(), above) and this
//   per-(representative,bin) remainder: one acos and two cos calls, and
//   the Q_lo algebra. The caller applies 'pow(.,kappa+1.5)' and the
//   M/min(M,.) wrapper once per cell, on this function's largest return
//   value over the cell's representatives -- bit-identical to applying
//   them per-representative-then-maxing, since both pow(.,kappa+1.5)
//   (kappa+1.5 > 0) and the M-scaled min are monotone nondecreasing.
//
Real AxiProposal::BoundRatio(const BoundRep & rep, const BinCap & cap,
                             Real ah, Real az) {

  Real alpha_c = 0.0;
  if (rep.Mr > 0.0) {
    Real cosang = rep.MHat.Dot(cap.Uc);
    if (cosang > 1.0)  cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;
    alpha_c = acos(cosang);
  }
  const Real alpha_lo = std::max((Real)0.0, alpha_c - cap.Delta);
  const Real alpha_hi = std::min((Real)Geometry::Pi, alpha_c + cap.Delta);

  const Real s_lo = rep.Mr*cos(alpha_hi) - rep.Rho;
  const Real s_hi = rep.Mr*cos(alpha_lo) + rep.Rho;

  const Real s2_hi = std::max(s_lo*s_lo, s_hi*s_hi);
  const Real s2_lo = (s_lo <= 0.0 && 0.0 <= s_hi) ?
                      (Real)0.0 : std::min(s_lo*s_lo, s_hi*s_hi);

  Real Q_lo;
  if (ah >= az) {
    Q_lo = std::max(ah*ah*rep.M2Lo - (ah*ah-az*az)*std::min(s2_hi,rep.M2Hi),
                    az*az*rep.M2Lo);
  } else {
    Q_lo = ah*ah*rep.M2Lo + (az*az-ah*ah)*s2_lo;
  }

  return rep.Numer / (1.0 + Q_lo);

}//
//


//////
// METHOD:  AxiProposal :: BuildBinTable()                     [const]
//
//   Builds one incidence bin's cumulative proposal-weight table (and
//   its majorant rate) for the two conversions sharing this bin's
//   u_c/delta ('reps0'/'reps1', 'S0'/'S1' -- see aniso_axi.hpp).
//
void AxiProposal::BuildBinTable(const R3::XYZ & u_c, Real delta,
                                const std::vector<BoundRep> & reps0,
                                const std::vector<BoundRep> & reps1,
                                Index repsPerCell,
                                const std::vector<Real> & S0,
                                const std::vector<Real> & S1,
                                BinTable & table, Real & rate) const {

  const BinCap cap = MakeBinCap(u_c, delta);

  std::vector<Real> maxRatio0(mNumCells, 0.0), maxRatio1(mNumCells, 0.0);

  const Index numReps = (Index) reps0.size();
  for (Index i = 0; i < numReps; i++) {
    const Index cell = i / repsPerCell;
    const Real r0 = BoundRatio(reps0[i], cap, mAh, mAz);
    if (r0 > maxRatio0[cell]) {maxRatio0[cell] = r0;}
    const Real r1 = BoundRatio(reps1[i], cap, mAh, mAz);
    if (r1 > maxRatio1[cell]) {maxRatio1[cell] = r1;}
  }

  std::vector<Real> bound0(mNumCells), bound1(mNumCells);
  Real sumBinBoundS = 0.0;
  for (Index c = 0; c < mNumCells; c++) {
    bound0[c] = std::min(mM, mM * pow(maxRatio0[c], mKappa+1.5));
    bound1[c] = std::min(mM, mM * pow(maxRatio1[c], mKappa+1.5));
    sumBinBoundS += bound0[c]*S0[c] + bound1[c]*S1[c];
  }

  table.assign(2*mNumCells, 0.0);
  Real running = 0.0;
  for (Index c = 0; c < mNumCells; c++) {
    const Real w = (S0[c] > 0.0) ?
      std::max(bound0[c]*(1.0+kBoundMargin)*S0[c], kWeightFloor*sumBinBoundS)
      : 0.0;
    running += w;
    table[c] = running;
  }
  for (Index c = 0; c < mNumCells; c++) {
    const Real w = (S1[c] > 0.0) ?
      std::max(bound1[c]*(1.0+kBoundMargin)*S1[c], kWeightFloor*sumBinBoundS)
      : 0.0;
    running += w;
    table[mNumCells+c] = running;
  }

  rate = running / (Real) mToaSize;

}//
//


//////
// METHOD:  AxiProposal :: Fold()                              [const]
//
//   Folds the local-vertical direction 'u' into a canonical incidence
//   bin (theta, chi for S; theta only for P), and records the mirror
//   flags (Mx, My) needed to map a table cell back to the actual cell.
//
AxiProposal::FoldResult AxiProposal::Fold(raytype intype, R3::XYZ u) const {

  if (u.z() < 0.0) {u = u.Negative();}

  const bool mx = (u.x() > 0.0);
  const Real ux = mx ? -u.x() : u.x();

  bool my = false;
  Real uy = 0.0;
  if (intype == RAY_S) {
    my = (u.y() < 0.0);
    uy = my ? -u.y() : u.y();
  }

  Real uz = u.z();
  if (uz > 1.0)  uz = 1.0;
  if (uz < -1.0) uz = -1.0;
  const Real theta = acos(uz);
  const Real sinTheta = sin(theta);

  Real chi = 0.0;
  if (intype == RAY_S && sinTheta > 1e-12) {
    chi = atan2(uy, -ux);
  }

  int thetaBin = (int) floor(theta*Geometry::RtoD/2.0);
  if (thetaBin < 0) thetaBin = 0;
  if (thetaBin > kThetaBins-1) thetaBin = kThetaBins-1;

  int chiBin = 0;
  if (intype == RAY_S) {
    chiBin = (int) floor(chi*Geometry::RtoD/10.0);
    if (chiBin < 0) chiBin = 0;
    if (chiBin > kChiBins-1) chiBin = kChiBins-1;
  }

  FoldResult r;
  r.ThetaBin = thetaBin;
  r.ChiBin   = chiBin;
  r.Mx       = mx;
  r.My       = my;
  return r;

}//
//


//////
// METHOD:  AxiProposal :: MapCell()                           [const]
//
Index AxiProposal::MapCell(Index cell, bool mx, bool my) const {
  Index c = cell;
  if (my) {c = mMyPerm[c];}
  if (mx) {c = mMxPerm[c];}
  return c;
}//
//


//////
// METHOD:  AxiProposal :: TableFor() / BinRateFor()           [const]
//
const AxiProposal::BinTable &
AxiProposal::TableFor(raytype intype, int thetaBin, int chiBin) const {
  return (intype == RAY_P) ? mPTable[thetaBin] : mSTable[thetaBin*kChiBins+chiBin];
}

Real AxiProposal::BinRateFor(raytype intype, int thetaBin, int chiBin) const {
  return (intype == RAY_P) ? mPBinRate[thetaBin] : mSBinRate[thetaBin*kChiBins+chiBin];
}


//////
// METHOD:  AxiProposal :: MajorantRate()                      [const]
//
Real AxiProposal::MajorantRate(raytype intype) const {
  return (intype == RAY_P) ? mPMajorant : mSMajorant;
}


//////
// METHOD:  AxiProposal :: BinRate()                           [const]
//
Real AxiProposal::BinRate(raytype intype, const R3::XYZ & u_rel) const {
  const FoldResult f = Fold(intype, u_rel);
  return BinRateFor(intype, f.ThetaBin, f.ChiBin);
}


//////
// METHOD:  AxiProposal :: DrawProposal()                      [const]
//
AxiProposal::Draw
AxiProposal::DrawProposal(raytype intype, const R3::XYZ & u_rel) const {

  const FoldResult f = Fold(intype, u_rel);
  const BinTable & table = TableFor(intype, f.ThetaBin, f.ChiBin);

  const Real total = table.back();
  const Real r = RandomUniform() * total;

  Index k1 = 0, k2 = (Index) table.size() - 1, k;
  while (k1 != k2) {
    k = (k1+k2) >> 1;
    if (r <= table[k]) {k2 = k;} else {k1 = k+1;}
  }

  const Index tableConv = k2 / mNumCells;   // 0 or 1
  const Index tableCell = k2 % mNumCells;
  const Real weight = table[k2] - ((k2 > 0) ? table[k2-1] : 0.0);

  const int g0 = (intype == RAY_P) ? 0 : 2;
  const int conv = g0 + (int) tableConv;

  const Index actualCell = MapCell(tableCell, f.Mx, f.My);
  const Index lo = actualCell * mCellBlock, hi = lo + mCellBlock;

  Draw d;
  d.Conv    = conv;
  d.Toa     = mPDists[conv].GetRandomIndexInRange(lo, hi);
  d.Weight  = weight;
  d.CellSum = mPDists[conv].GetRangeMagnitude(lo, hi);
  return d;

}//
//


//////
// METHOD:  AxiProposal :: AcceptanceBound()                   [const]
//
Real AxiProposal::AcceptanceBound(raytype intype, const R3::XYZ & u_rel,
                                  int conv, Index cell) const {

  const FoldResult f = Fold(intype, u_rel);
  const BinTable & table = TableFor(intype, f.ThetaBin, f.ChiBin);

  const int g0 = (intype == RAY_P) ? 0 : 2;
  const Index tableConv = (Index) (conv - g0);
  const Index tableCell = MapCell(cell, f.Mx, f.My);  // self-inverse

  const Index idx = tableConv*mNumCells + tableCell;
  const Real weight = table[idx] - ((idx > 0) ? table[idx-1] : 0.0);

  const Index lo = cell * mCellBlock, hi = lo + mCellBlock;
  const Real cellSum = mPDists[conv].GetRangeMagnitude(lo, hi);

  if (cellSum <= 0.0) {return 0.0;}
  return weight / cellSum;

}//
//


//////
// METHOD:  AxiProposal :: CellOf()                            [const]
//
Index AxiProposal::CellOf(Index toa) const {
  return toa / mCellBlock;
}
