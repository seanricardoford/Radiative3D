// aniso_axi.cpp
//
#include <algorithm>    /* std::min, std::max */
#include <cmath>        /* sin, cos, acos, atan2, sqrt, pow, floor, log */
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

  // D = log4(toa.size()/20):
  const int D = (int) ( log( (double) mToaSize / 20.0 ) / log(4.0) + 0.5 );

  const int coarseDeg = std::min(D, kCoarseDegree);
  mNumCells  = 20 * Pow4(coarseDeg);
  mCellBlock = mToaSize / mNumCells;

  if (mToaSize != mNumCells * mCellBlock) {
    // The cell/block arithmetic throughout this class (hierarchical
    // TesselSphere ordering, contiguous per-cell blocks of fine
    // directions) assumes toa.size() is exactly 20*4^D for some D, as
    // TESS_ICO always produces. A non-icosahedral take-off-angle set
    // would silently drop the tail directions from every cell's range,
    // which is a silent bias rather than a loud failure. Production
    // always uses TESS_ICO; this is a belt-and-braces guard.
    throw std::logic_error(
      "AxiProposal requires an icosahedral (TESS_ICO) take-off-angle "
      "set: toa.size() is not a multiple of mNumCells");
  }

  BuildMirrorPerms(coarseDeg);

  // Representative directions (for bound construction only -- these
  // don't need to survive past the constructor):
  const int repDeg = std::min(D, kRepDegree);
  Index numReps, repBlock;
  if (repDeg == D) {
    numReps  = mToaSize;
    repBlock = 1;
  } else {
    numReps  = 20 * Pow4(repDeg);
    repBlock = mToaSize / numReps;
  }
  const Index repsPerCell = numReps / mNumCells;

  std::vector<R3::XYZ> repCenters(numReps);
  std::vector<Real>    repRadius(numReps, 0.0);

  if (repDeg == D) {
    for (Index i = 0; i < numReps; i++) {
      repCenters[i] = toa[i].XYZ();
    }
    // repRadius[] stays all-zero: each representative is its own (and
    // only) fine direction.
  } else {
    S2::TesselSphere repTess(S2::TESS_ICO, repDeg);
    for (Index i = 0; i < numReps; i++) {
      repCenters[i] = repTess[i].XYZ();
    }
    for (Index i = 0; i < numReps; i++) {
      Real maxAngle = 0.0;
      const Index base = i * repBlock;
      for (Index k = base; k < base + repBlock; k++) {
        Real c = repCenters[i].Dot(toa[k].XYZ());
        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;
        const Real ang = acos(c);
        if (ang > maxAngle) {maxAngle = ang;}
      }
      repRadius[i] = maxAngle;
    }
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
    BuildBinTable(RAY_P, u_c, delta, repCenters, repRadius, repsPerCell,
                 mPTable[b], mPBinRate[b]);
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
      BuildBinTable(RAY_S, u_c, delta, repCenters, repRadius, repsPerCell,
                   mSTable[idx], mSBinRate[idx]);
    }
  }
  mSMajorant = 0.0;
  for (std::size_t i = 0; i < mSBinRate.size(); i++) {
    mSMajorant = std::max(mSMajorant, mSBinRate[i]);
  }

}//
//


//////
// METHOD:  AxiProposal :: BuildMirrorPerms()         [Helper Function]
//
//   Populates mMxPerm/mMyPerm: for each cell, the index of the cell
//   whose center is the x- (resp. y-) mirror image of this cell's
//   center, matched within 1e-9. Throws std::logic_error if any cell
//   has no match (the icosahedral tessellation used throughout this
//   class is, by construction, symmetric under negating x or y alone --
//   see geom_s2.cpp's init_icosahedron() -- so a missing match signals
//   a real inconsistency, not just tessellation noise).
//
void AxiProposal::BuildMirrorPerms(int coarseDeg) {

  std::vector<R3::XYZ> centers(mNumCells);

  // NOTE: Building a fresh degree-coarseDeg TesselSphere always gives
  // exactly the same direction set as grouping 'toa' into mNumCells
  // cells would (whether coarseDeg==D or coarseDeg<D), because
  // TesselSphere's recursive quadrisection preserves each cell's three
  // corner nodes at every finer degree: see geom_s2.cpp's
  // Triangle::subdivide().
  S2::TesselSphere cellTess(S2::TESS_ICO, coarseDeg);
  for (Index j = 0; j < mNumCells; j++) {
    centers[j] = cellTess[j].XYZ();
  }

  mMxPerm.resize(mNumCells);
  mMyPerm.resize(mNumCells);

  const Real tol2 = 1e-9*1e-9;

  for (Index j = 0; j < mNumCells; j++) {

    const R3::XYZ mx_target(-centers[j].x(), centers[j].y(), centers[j].z());
    const R3::XYZ my_target( centers[j].x(),-centers[j].y(), centers[j].z());

    Index bestMx = 0, bestMy = 0;
    Real bestMxD2 = 1.0e300, bestMyD2 = 1.0e300;

    for (Index k = 0; k < mNumCells; k++) {
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

    mMxPerm[j] = bestMx;
    mMyPerm[j] = bestMy;
  }

}//
//


//////
// METHOD:  AxiProposal :: ComputeBound()                      [const]
//
//   Rigorous upper bound on AnisoRatio(m) (see ScatterParams::
  // AnisoRatio()) over every fine direction within angle 'r_i' of
//   representative direction 't' (scaled by 'k_out'), and every local-
//   vertical direction within angle 'delta' of 'u_c', for a conversion
//   with the given k_in/k_out. See docs/superpowers/specs/
//   2026-10-16-anisotropic-scale-lengths-design.md and the task-7
//   brief for the derivation.
//
Real AxiProposal::ComputeBound(const R3::XYZ & t, Real r_i,
                               Real k_in, Real k_out,
                               const R3::XYZ & u_c, Real delta) const {

  const R3::XYZ zhat(0.0, 0.0, 1.0);   // Forward (incoming-ray) direction,
                                       // in the local ray frame.
  const R3::XYZ m_rep = t.ScaledBy(k_out) + zhat.ScaledBy(-k_in);
  const Real mr  = m_rep.Mag();
  const Real rho = k_out * 2.0 * sin(r_i/2.0);

  Real alpha_c = 0.0;
  if (mr > 0.0) {
    Real cosang = m_rep.Unit().Dot(u_c);
    if (cosang > 1.0)  cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;
    alpha_c = acos(cosang);
  }
  const Real alpha_lo = std::max((Real)0.0, alpha_c - delta);
  const Real alpha_hi = std::min((Real)Geometry::Pi, alpha_c + delta);

  const Real s_lo = mr*cos(alpha_hi) - rho;
  const Real s_hi = mr*cos(alpha_lo) + rho;

  const Real s2_hi = std::max(s_lo*s_lo, s_hi*s_hi);
  const Real s2_lo = (s_lo <= 0.0 && 0.0 <= s_hi) ?
                      (Real)0.0 : std::min(s_lo*s_lo, s_hi*s_hi);

  const Real mlo   = std::max(mr-rho, (Real)0.0);
  const Real m2_lo = mlo*mlo;
  const Real mhi   = mr+rho;
  const Real m2_hi = mhi*mhi;

  Real Q_lo;
  if (mAh >= mAz) {
    Q_lo = std::max(mAh*mAh*m2_lo - (mAh*mAh-mAz*mAz)*std::min(s2_hi,m2_hi),
                    mAz*mAz*m2_lo);
  } else {
    Q_lo = mAh*mAh*m2_lo + (mAz*mAz-mAh*mAh)*s2_lo;
  }

  const Real ratioBase = (1.0 + mAMin*mAMin*m2_hi) / (1.0 + Q_lo);
  const Real bound = mM * pow(ratioBase, mKappa+1.5);

  return std::min(mM, bound);

}//
//


//////
// METHOD:  AxiProposal :: BuildBinTable()                     [const]
//
//   Builds one incidence bin's cumulative proposal-weight table (and
//   its majorant rate) for the two conversions available to 'intype'.
//
void AxiProposal::BuildBinTable(raytype intype, const R3::XYZ & u_c,
                                Real delta,
                                const std::vector<R3::XYZ> & repCenters,
                                const std::vector<Real> & repRadius,
                                Index repsPerCell,
                                BinTable & table, Real & rate) const {

  const Real k_in = (intype == RAY_P) ? mKP : mKS;
  // conv0's output is always P (GPP for intype P, GSP for intype S);
  // conv1's output is always S (GPS or GSS):
  const Real k_out0 = mKP;
  const Real k_out1 = mKS;
  const int g0 = (intype == RAY_P) ? 0 : 2;
  const int g1 = g0 + 1;

  std::vector<Real> bound0(mNumCells, 0.0), bound1(mNumCells, 0.0);

  const Index numReps = (Index) repCenters.size();
  for (Index i = 0; i < numReps; i++) {
    const Index cell = i / repsPerCell;
    const Real b0 = ComputeBound(repCenters[i], repRadius[i], k_in, k_out0,
                                 u_c, delta);
    if (b0 > bound0[cell]) {bound0[cell] = b0;}
    const Real b1 = ComputeBound(repCenters[i], repRadius[i], k_in, k_out1,
                                 u_c, delta);
    if (b1 > bound1[cell]) {bound1[cell] = b1;}
  }

  std::vector<Real> S0(mNumCells), S1(mNumCells);
  Real sumBinBoundS = 0.0;
  for (Index c = 0; c < mNumCells; c++) {
    const Index lo = c*mCellBlock, hi = lo+mCellBlock;
    S0[c] = mPDists[g0].GetRangeMagnitude(lo, hi);
    S1[c] = mPDists[g1].GetRangeMagnitude(lo, hi);
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
