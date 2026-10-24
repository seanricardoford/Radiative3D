// aniso_axi.cpp
//
#include <algorithm>    /* std::min, std::max */
#include <cmath>        /* sin, cos, acos, atan2, sqrt, pow, floor, log */
#include <cstdint>      /* std::uint64_t: FingerprintOf() chunking */
#include <map>          /* std::map: geometry cache, keyed by TOA count */
#include <mutex>        /* std::mutex, std::lock_guard: geometry cache */
#include <stdexcept>    /* std::logic_error */
#include "aniso_axi.hpp"
#include "parallel.hpp"
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

  // ::: ThetaEdgeCosines() / ChiEdgeCosines() / ChiEdgeSines() :::
  //
  //   Precomputed edge-angle cosines/sines used by AxiProposal::Locate()
  //   to fold a local-vertical direction into its incidence bin by
  //   comparison alone -- no acos, atan2, or sin call on the hot path
  //   (design decision 4). Each is built once (on first use; C++11
  //   guarantees thread-safe function-local static initialization) and
  //   reused for the life of the program: these depend only on the fixed
  //   bin widths (2 degrees for theta, 10 degrees for chi), never on any
  //   ScatterParams or take-off-angle set.
  //
  //   ThetaEdgeCosines()[b] == cos(2*b degrees), b in [0,kThetaBins]
  //   (46 entries): theta bin b covers the (acos-monotonicity) interval
  //   cos(2*(b+1) deg) < u_z <= cos(2*b deg).
  //
  //   ChiEdgeCosines()[c] == cos(10*c degrees),
  //   ChiEdgeSines()[c]   == sin(10*c degrees), c in [0,kChiBins-1] (9
  //   entries each): chi bin c is the largest c with
  //   v*cos(10*c deg) >= w*sin(10*c deg) (see Locate()).
  //
  const std::vector<Real> & ThetaEdgeCosines() {
    static const std::vector<Real> edges = [] {
      std::vector<Real> e(AxiProposal::kThetaBins + 1);
      for (int b = 0; b <= AxiProposal::kThetaBins; ++b) {
        e[b] = cos(2.0 * (Real) b * Geometry::DtoR);
      }
      return e;
    }();
    return edges;
  }

  const std::vector<Real> & ChiEdgeCosines() {
    static const std::vector<Real> c = [] {
      std::vector<Real> v(AxiProposal::kChiBins);
      for (int i = 0; i < AxiProposal::kChiBins; ++i) {
        v[i] = cos(10.0 * (Real) i * Geometry::DtoR);
      }
      return v;
    }();
    return c;
  }

  const std::vector<Real> & ChiEdgeSines() {
    static const std::vector<Real> s = [] {
      std::vector<Real> v(AxiProposal::kChiBins);
      for (int i = 0; i < AxiProposal::kChiBins; ++i) {
        v[i] = sin(10.0 * (Real) i * Geometry::DtoR);
      }
      return v;
    }();
    return s;
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

  // ::: FingerprintOf(toa, threads) :::
  //
  //   A cheap scalar fingerprint of a take-off-angle set: controller
  //   ruling R6 (amending R2) -- the index-weighted sum
  //   Sum_k (k+1)*(theta_k + 2*phi_k), k in index order, 0-based, over
  //   the take-off-angle set's *stored* theta/phi (S2::ThetaPhi::Theta()/
  //   Phi()), with no trigonometry. This is only a sanity check guarding
  //   the geometry cache below against reuse across two distinct
  //   take-off-angle sets that happen to share a size (not a security
  //   hash) -- production builds a single TesselSphere per run and
  //   reuses it for every Scatterer (see PhononSource::Set_TOA_Array()),
  //   so a collision can't arise there.
  //
  //   R2's version of this fingerprint, Sum_k (k+1)*(x_k+2*y_k+3*z_k)
  //   over toa[k].XYZ(), replaced an even earlier plain (unweighted,
  //   equally-coefficiented) component sum that a standalone check
  //   showed to be a poor discriminator: TESS_ICO's vertex directions
  //   sum to (numerically) zero by icosahedral symmetry, and since
  //   rotation is linear, a rotated copy's sum is *also* (numerically)
  //   zero. Per-index weighting breaks that cancellation (a rotation or
  //   reordering changes which direction sits at which index, and a
  //   different index carries a different weight); R2's x/2y/3z
  //   coefficients additionally broke the remaining x<->y<->z
  //   degeneracy, which R6's theta/2*phi analogue doesn't need to worry
  //   about in the first place -- theta and phi aren't interchangeable
  //   the way x/y/z are, so there's no equivalent degeneracy to break,
  //   and summing raw polar angles (every theta_k >= 0) carries none of
  //   x_k/y_k/z_k's own signed-component cancellation either. R6's sole
  //   purpose (per the controller) is dropping every sin/cos call
  //   XYZ() made -- profiling (task-5b-brief.md) found this loop, called
  //   fresh for every one of a run's axisymmetric Scatterers regardless
  //   of geometry-cache hits (see GetOrBuildAxiGeometry() below), costing
  //   about 1.7s of wall time alone at a degree-9 take-off-angle count.
  //   A standalone check confirmed R6 stays just as discriminating as
  //   R2 for this task's two regression cases: a 0.01 rad rotation about
  //   z of the degree-5 TESS_ICO set (test_axi_geometry_cache_mismatch_
  //   rotation) shifts every direction's phi by the same constant, so
  //   the fingerprint changes by exactly 2*0.01*Sum_k(k+1) -- far above
  //   floating-point noise, with none of R2's x/y/z-cancellation risk to
  //   begin with; the per-index phi perturbation used by
  //   test_axi_geometry_cache_mismatch is similarly non-cancelling.
  //
  //   Computed as 64 fixed-size chunks (chunk c covers
  //   [toa.size()*c/64, toa.size()*(c+1)/64), the division done in
  //   64-bit so it can't overflow for any toa.size() that itself fits in
  //   an Index), each summed by one ParallelFor(0,64,threads,...) task
  //   into its own part[c], combined by a final serial loop over chunk
  //   index 0..63 -- so the chunk boundaries (fixed at 64, independent of
  //   'threads') and the fixed combination order make the result
  //   independent of 'threads', exactly as every other construction-time
  //   table in this class is (see the class-level doc comment on
  //   AxiProposal's constructor).
  //
  Real FingerprintOf(const S2::S2Set & toa, unsigned threads) {
    const Index n = (Index) toa.size();
    const Index kChunks = 64;
    std::vector<Real> part(kChunks, 0.0);
    ParallelFor(0, kChunks, threads, [&](Index c) {
      const Index lo = (Index) ((std::uint64_t) n * c / kChunks);
      const Index hi = (Index) ((std::uint64_t) n * (c+1) / kChunks);
      Real s = 0.0;
      for (Index k = lo; k < hi; k++) {
        s += (Real) (k+1) * (toa[k].Theta() + 2.0*toa[k].Phi());
      }
      part[c] = s;
    });
    Real sum = 0.0;
    for (Index c = 0; c < kChunks; c++) {sum += part[c];}
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

  // ::: BuildAxiGeometry(toa, toaSize, fp, threads, geom) :::
  //
  //   Does the actual (expensive) geometry construction, writing into
  //   'geom' (which the caller owns -- either a fresh map slot on a true
  //   cache miss, or a construction-local scratch object on a
  //   fingerprint mismatch; see GetOrBuildAxiGeometry() below). Split
  //   out from GetOrBuildAxiGeometry() so both call sites share one
  //   body instead of risking two copies drifting apart.
  //
  //   'threads': the number of std::thread workers (ParallelFor(),
  //   parallel.hpp) used to build the representative-radius array and
  //   the mirror-cell permutations below (Task 4b). 'geom' comes out
  //   bit-identical regardless of 'threads': each representative or
  //   cell is computed by exactly one task, writing only its own
  //   RepRadius[i]/MxPerm[j]/MyPerm[j] entry, with no cross-task
  //   reduction (see docs/superpowers/specs/2026-10-17-axisymmetric-
  //   sampler-performance-design.md, decision 2).
  //
  void BuildAxiGeometry(const S2::S2Set & toa, Index toaSize, Real fp,
                        unsigned threads, CachedAxiGeometry & geom) {

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
      // Each representative's angular radius -- the largest angle from
      // its centre to any of its 'repBlock' fine take-off-angle
      // directions (one acos per (representative, fine-direction) pair)
      // -- is independent across representatives: each task below
      // writes only its own geom.RepRadius[i].
      ParallelFor(0, numReps, threads, [&](Index i) {
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
      });
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

    // Brute-force nearest-centre match for each cell's x- and y-mirror
    // image, independent across cells: each task below writes only its
    // own geom.MxPerm[j]/MyPerm[j], with no cross-task reduction. A "no
    // match" throw is preserved -- ParallelFor rethrows the exception
    // from the numerically lowest-indexed failing task (parallel.hpp).
    ParallelFor(0, numCells, threads, [&](Index j) {

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
    });

  }//

  // ::: GetOrBuildAxiGeometry(toa, threads, scratch) :::
  //
  //   Returns the geometry for 'toa': a cache hit returns a reference
  //   into the cache (built, possibly long ago, under whatever thread
  //   count that first construction used -- see the AxiProposal
  //   constructor's own doc comment on why that's fine); a true miss (no
  //   entry yet for this size) builds directly into a fresh cache slot,
  //   using 'threads' workers, and returns a reference to it; a
  //   fingerprint mismatch (an entry exists for this size, but isn't
  //   'toa') builds into the caller-owned 'scratch' instead, also using
  //   'threads' workers, and returns a reference to *that*, leaving the
  //   cache untouched. 'scratch' must outlive the returned reference in
  //   every case -- the constructor keeps it as a local for exactly as
  //   long as it uses the result.
  //
  const CachedAxiGeometry & GetOrBuildAxiGeometry(const S2::S2Set & toa,
                                                  unsigned threads,
                                                  CachedAxiGeometry & scratch) {

    static std::mutex mtx;
    static std::map<Index, CachedAxiGeometry> cache;

    const Index toaSize = (Index) toa.size();
    const Real  fp       = FingerprintOf(toa, threads);

    std::lock_guard<std::mutex> lock(mtx);

    std::map<Index, CachedAxiGeometry>::iterator it = cache.find(toaSize);
    if (it != cache.end()) {
      if (it->second.Fingerprint == fp) {
        return it->second;
      }
      BuildAxiGeometry(toa, toaSize, fp, threads, scratch);
      return scratch;
    }

    CachedAxiGeometry & slot = cache[toaSize];
    BuildAxiGeometry(toa, toaSize, fp, threads, slot);
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
                         const std::vector<ProbDist> & pdists,
                         unsigned threads) :
  mPDists(pdists)
{

  // Force Locate()'s support tables (ThetaEdgeCosines()/ChiEdgeCosines()/
  // ChiEdgeSines(), above) to build now, during construction, rather than
  // lazily on this process's first Locate() call. Locate() isn't called
  // by anything in this constructor -- only by BinRate()/DrawProposal()/
  // AcceptanceBound(), all per-scattering-event queries -- so left alone,
  // the three functions' function-local statics would otherwise first
  // build during propagation (the first tentative collision any
  // axisymmetric Scatterer handles), a one-time write to process-wide
  // shared state that "Scatterers are read-only once constructed.
  // Propagation writes no shared state." rules out, even though C++11
  // guarantees the write itself is thread-safe. These three tables
  // depend only on the fixed bin widths (2 degrees theta, 10 degrees
  // chi), never on 'par'/'toa'/'pdists', so touching them here changes
  // neither their contents nor which call ends up building them the
  // *first* time overall (still whichever AxiProposal, of however many a
  // run constructs, happens to run first) -- it only guarantees that
  // first time is during some Scatterer's construction, never during
  // propagation. After the first call, each is just a guard-variable
  // check -- negligible cost on every later AxiProposal built.
  (void) ThetaEdgeCosines();
  (void) ChiEdgeCosines();
  (void) ChiEdgeSines();

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
  const CachedAxiGeometry & geom = GetOrBuildAxiGeometry(toa, threads, scratch);
  mMxPerm = geom.MxPerm;
  mMyPerm = geom.MyPerm;

  // Per-representative bound-kernel data (BoundRep), precomputed once
  // for each of the four (k_in,k_out) pairs the P and S bin tables need
  // -- P bins use (kP,kP) and (kP,kS); S bins use (kS,kP) and (kS,kS) --
  // instead of once per (representative, bin) as the old ComputeBound()
  // did. Independent across representatives (each task writes only its
  // own repsPP[i]/repsPS[i]/repsSP[i]/repsSS[i]), so this runs via
  // ParallelFor on up to 'threads' workers (Task 4b). See
  // PrepareBoundRep()/BoundRatio() and BuildBinTable() below.
  const Index numReps = (Index) geom.RepCenters.size();
  RepSoA soaPP, soaPS, soaSP, soaSS;
  {
    // Scoped so the four vector<BoundRep> temporaries below (~5 MB each
    // at kRepDegree's 81,920 representatives) are freed as soon as their
    // data is copied into the RepSoA arrays above (also ~5 MB each),
    // rather than staying alive alongside them for the rest of the
    // constructor: keeping both sets live through the 450-bin loop
    // further down would cost ~42 MB of temporaries, over the ~32 MB/
    // scatterer budget (docs/superpowers/specs/2026-10-17-axisymmetric-
    // sampler-performance-design.md); freed here, only one ~21 MB set is
    // ever live at a time.
    std::vector<BoundRep> repsPP(numReps), repsPS(numReps),
                          repsSP(numReps), repsSS(numReps);
    ParallelFor(0, numReps, threads, [&](Index i) {
      const R3::XYZ & t = geom.RepCenters[i];
      const Real      r = geom.RepRadius[i];
      repsPP[i] = PrepareBoundRep(t, r, mKP, mKP, mAMin);
      repsPS[i] = PrepareBoundRep(t, r, mKP, mKS, mAMin);
      repsSP[i] = PrepareBoundRep(t, r, mKS, mKP, mAMin);
      repsSS[i] = PrepareBoundRep(t, r, mKS, mKS, mAMin);
    });

    // Structure-of-arrays copies of the four BoundRep vectors above, fed
    // to BuildBinTable()/RatioRun() instead of the vector<BoundRep>
    // arrays themselves so that function's hot loop reads flat,
    // contiguous Real arrays and can vectorize (see aniso_axi.hpp's
    // RepSoA/RatioRun doc comments). Assign() is a plain column-wise
    // copy, so this adds no arithmetic and so no rounding difference.
    soaPP.Assign(repsPP); soaPS.Assign(repsPS);
    soaSP.Assign(repsSP); soaSS.Assign(repsSS);
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

  // Incidence bins: P (theta only; canonical frame fixes chi == 0), then
  // S (theta and chi, flattened to a single index), over one combined
  // index space of kThetaBins + kThetaBins*kChiBins = 450 bins, P first
  // then S in today's per-table index order (Task 4b: previously two
  // separate ParallelFor calls). Each bin's table/rate (mPTable[b]/
  // mPBinRate[b], or the mSTable/mSBinRate equivalent) is written by
  // exactly one task and depends only on its own bin index and the
  // read-only data prepared above, so every bin runs via a single
  // ParallelFor on up to 'threads' workers -- halving the per-
  // AxiProposal thread spawn/join overhead from two calls to one (docs/
  // superpowers/specs/2026-10-17-axisymmetric-sampler-performance-
  // design.md, decision 2). mPMajorant/mSMajorant are then each reduced
  // serially, in ascending bin-index order, exactly as before, so both
  // tables are bit-identical for any 'threads'.
  mPTable.resize(kThetaBins);
  mPBinRate.resize(kThetaBins);
  mSTable.resize(kThetaBins*kChiBins);
  mSBinRate.resize(kThetaBins*kChiBins);

  const Index numPBins = (Index) kThetaBins;
  const Index numBins  = numPBins + (Index) (kThetaBins*kChiBins);

  ParallelFor(0, numBins, threads, [&](Index idx) {
    if (idx < numPBins) {
      const int b = (int) idx;
      const Real theta_lo = (2.0*b)     * Geometry::DtoR;
      const Real theta_hi = (2.0*b+2.0) * Geometry::DtoR;
      const Real theta_c  = 0.5*(theta_lo+theta_hi);
      const Real delta    = 0.5*(theta_hi-theta_lo);
      const R3::XYZ u_c(-sin(theta_c), 0.0, cos(theta_c));
      BuildBinTable(u_c, delta, soaPP, soaPS, geom.RepsPerCell,
                   S0_P, S1_P, mPTable[b], mPBinRate[b]);
    } else {
      const Index sidx = idx - numPBins;
      const int b = (int) (sidx / kChiBins);
      const int c = (int) (sidx % kChiBins);
      const Real theta_lo = (2.0*b)     * Geometry::DtoR;
      const Real theta_hi = (2.0*b+2.0) * Geometry::DtoR;
      const Real theta_c  = 0.5*(theta_lo+theta_hi);
      const Real chi_lo = (10.0*c)      * Geometry::DtoR;
      const Real chi_hi = (10.0*c+10.0) * Geometry::DtoR;
      const Real chi_c  = 0.5*(chi_lo+chi_hi);
      const Real delta = 0.5*(theta_hi-theta_lo)
                        + sin(theta_hi)*0.5*(chi_hi-chi_lo);
      const R3::XYZ u_c(-sin(theta_c)*cos(chi_c),
                        sin(theta_c)*sin(chi_c),
                        cos(theta_c));
      BuildBinTable(u_c, delta, soaSP, soaSS, geom.RepsPerCell,
                   S0_S, S1_S, mSTable[sidx], mSBinRate[sidx]);
    }
  });

  mPMajorant = 0.0;
  for (int b = 0; b < kThetaBins; b++) {
    mPMajorant = std::max(mPMajorant, mPBinRate[b]);
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
//   per-(representative,bin) remainder: the Q_lo algebra, now fed by a
//   trig-free evaluation of cos(alpha_lo)/cos(alpha_hi) (design decision
//   3) instead of one acos and two cos calls. The caller applies
//   'pow(.,kappa+1.5)' and the M/min(M,.) wrapper once per cell, on this
//   function's largest return value over the cell's representatives --
//   bit-identical to applying them per-representative-then-maxing, since
//   both pow(.,kappa+1.5) (kappa+1.5 > 0) and the M-scaled min are
//   monotone nondecreasing.
//
//   With alpha_c = angle(MHat,Uc) (0 when Mr==0, matching the acos
//   kernel's 'alpha_c=0.0' default -- m_rep.Unit() is never evaluated in
//   that case), alpha_lo = max(0,alpha_c-Delta), alpha_hi =
//   min(pi,alpha_c+Delta):
//     cos(alpha_lo) == 1          whenever alpha_c <= Delta,
//                                 i.e. c=cos(alpha_c) >= cos(Delta);
//                   == cos(alpha_c-Delta) = c*CosDelta + s*SinDelta
//                                 otherwise;
//     cos(alpha_hi) == -1         whenever alpha_c+Delta >= pi,
//                                 i.e. c <= -cos(Delta);
//                   == cos(alpha_c+Delta) = c*CosDelta - s*SinDelta
//                                 otherwise;
//   using s = sin(alpha_c) = sqrt(1-c^2) >= 0 (alpha_c in [0,pi]). No
//   acos, cos or sin call remains in this function.
//
// This is the scalar reference implementation: BuildBinTable() (below)
// calls only RatioRun(), never this function directly, so BoundRatio()
// has no production caller any more -- it exists solely as the fixed,
// unvectorized target tests/anisotropy_tests.cpp's
// test_axi_ratio_run_matches_scalar and axi_kernel_matches_reference
// check against. RatioRun() must reproduce this function's result
// bit-for-bit (not merely within a tolerance) for every representative;
// the rest of this comment block explains how that's achieved despite
// vectorization, and the bug (an out-of-line Dot() call compiled under
// different ambient flags) that first broke it.
//
// This function and RatioRun() (below) carry a matching GCC
// optimize("fp-contract=off") attribute so neither ever contracts a
// multiply-add into a fused instruction: GCC's default for this file's
// build flags is -ffp-contract=fast (confirmed by -Q --help=optimizers),
// and disassembly showed BOTH this scalar function and RatioRun's
// vectorized loop already contain real fmadd/fmsub/fnmsub(/NEON
// fmla/fmls) instructions under that default -- just not, for at least
// one representative/bin found by test_axi_ratio_run_matches_scalar
// (tests/anisotropy_tests.cpp), the *same* ones, so the two rounded
// differently for that input. Forcing both to fully expand every
// multiply-add (never fusing) removes the ambiguity: per IEEE 754,
// scalar and SIMD add/sub/mul/div/sqrt are each required to round
// identically lane-for-lane, so two functions evaluating the same
// source expressions with no fused contraction anywhere are
// bit-identical by construction, independent of whether either one
// vectorizes. See the exactness note on RatioRun() below.
#if defined(__GNUC__) && !defined(__clang__)
#define R3D_EXACT_FP_ATTR \
  __attribute__((optimize("no-math-errno","no-trapping-math","fp-contract=off")))
#else
#define R3D_EXACT_FP_ATTR
#endif

R3D_EXACT_FP_ATTR
Real AxiProposal::BoundRatio(const BoundRep & rep, const BinCap & cap,
                             Real ah, Real az) {

  Real c = 1.0;              // cos(alpha_c); the Mr==0 default (see above).
  if (rep.Mr > 0.0) {
    // Written out component-wise, rather than rep.MHat.Dot(cap.Uc),
    // because R3::XYZ::Dot() (geom_r3.hpp) is an ordinary inline method
    // with its own separately-emitted out-of-line copy for call sites
    // GCC doesn't inline -- compiled under this translation unit's
    // *default* flags (-ffp-contract=fast; confirmed by disassembly:
    // that copy contains two fmadd instructions), not under this
    // function's R3D_EXACT_FP_ATTR. A BoundRatio() call that happened to
    // reach that copy would therefore round this dot product using
    // fused multiply-add while RatioRun()'s own (always manually
    // inlined, always attributed) version of the same sum never does --
    // precisely the mismatch test_axi_ratio_run_matches_scalar
    // (tests/anisotropy_tests.cpp) caught. Dot() itself is untouched
    // (every other call site across the codebase keeps today's
    // behavior); only this one call is expanded so it's governed by the
    // same attribute as the rest of this function.
    c = rep.MHat.x()*cap.Uc.x() + rep.MHat.y()*cap.Uc.y()
      + rep.MHat.z()*cap.Uc.z();
    if (c > 1.0)  c = 1.0;
    if (c < -1.0) c = -1.0;
  }
  const Real s = sqrt(std::max((Real)0.0, 1.0 - c*c));  // sin(alpha_c) >= 0.

  const Real cos_alpha_lo = (c >= cap.CosDelta) ?
    (Real)1.0 : c*cap.CosDelta + s*cap.SinDelta;
  const Real cos_alpha_hi = (c <= -cap.CosDelta) ?
    (Real)-1.0 : c*cap.CosDelta - s*cap.SinDelta;

  const Real s_lo = rep.Mr*cos_alpha_hi - rep.Rho;
  const Real s_hi = rep.Mr*cos_alpha_lo + rep.Rho;

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
// METHOD:  AxiProposal :: RepSoA :: Assign()
//
//   Column-wise copy of 'reps' into the eight parallel arrays: no
//   arithmetic, so this cannot itself be a source of any rounding
//   difference between RatioRun() and BoundRatio().
//
void AxiProposal::RepSoA::Assign(const std::vector<BoundRep> & reps) {
  const std::size_t n = reps.size();
  Mx.resize(n); My.resize(n); Mz.resize(n); Mr.resize(n);
  Rho.resize(n); M2Lo.resize(n); M2Hi.resize(n); Numer.resize(n);
  for (std::size_t i = 0; i < n; i++) {
    Mx[i]    = reps[i].MHat.x();
    My[i]    = reps[i].MHat.y();
    Mz[i]    = reps[i].MHat.z();
    Mr[i]    = reps[i].Mr;
    Rho[i]   = reps[i].Rho;
    M2Lo[i]  = reps[i].M2Lo;
    M2Hi[i]  = reps[i].M2Hi;
    Numer[i] = reps[i].Numer;
  }
}//
//


// GCC's auto-vectorizer needs permission to reorder the libm-adjacent
// comparisons/selects below without worrying about errno or trapping FP
// exceptions (neither of which this codebase ever inspects); combined,
// here, with the same fp-contract=off as BoundRatio() above, for the
// exactness reasons given there. Granted locally, via a function
// attribute, rather than globally via a command-line flag, so it
// affects only these two definitions. Clang (and anything else) gets
// the unannotated definition: the attribute is a GCC extension, and
// this loop vectorizes under Clang's defaults anyway.
//
//////
// METHOD:  AxiProposal :: RatioRun()                          [static]
//
//   BoundRatio(), run over representatives [lo,hi) of 'reps' -- same
//   per-element arithmetic, rewritten so the loop vectorizes:
//     - 'c' (cos(alpha_c)) is computed unconditionally as the clamped
//       dot product, then replaced with the Mr==0 default (1.0) by a
//       final select, instead of branching around the Dot() call --
//       correctness-preserving because the discarded dot-product value
//       is exactly 0 when Mr==0 (MHat is the zero vector then) and is
//       never used in that case either way.
//     - Both cap-edge candidate values (what BoundRatio() calls
//       cos_alpha_lo/cos_alpha_hi's non-clamped branch) are computed
//       unconditionally, before the clamp select, rather than only in
//       an 'else' branch.
//     - std::max()/std::min() are spelled as explicit ternaries (with
//       the same operand order and tie-break as the std:: versions --
//       see the exactness note below), and the s2_lo guard uses '&'
//       (bitwise, both sides always evaluated) rather than '&&'.
//   None of this changes any individual result: every selection below
//   is a choice between two already-computed exact values (a
//   comparison, never an arithmetic operation), so it introduces no
//   rounding beyond what BoundRatio() itself performs, regardless of
//   which branch the scalar function would have taken. The ternaries
//   match std::max(a,b)/std::min(a,b)'s own tie-break (returns 'a' when
//   a==b) exactly, so the two are interchangeable bit-for-bit -- not
//   merely mathematically equivalent -- for every input, including
//   ties.
//
template <bool AhGeAz>
R3D_EXACT_FP_ATTR
void AxiProposal::RatioRun(const RepSoA & reps, Index lo, Index hi,
                           const BinCap & cap, Real ah, Real az,
                           Real * ratio) {

  const Real * const mx    = &reps.Mx[0];
  const Real * const my    = &reps.My[0];
  const Real * const mz    = &reps.Mz[0];
  const Real * const mr    = &reps.Mr[0];
  const Real * const rho   = &reps.Rho[0];
  const Real * const m2lo  = &reps.M2Lo[0];
  const Real * const m2hi  = &reps.M2Hi[0];
  const Real * const numer = &reps.Numer[0];
  Real * const out = ratio;

  const Real ux = cap.Uc.x(), uy = cap.Uc.y(), uz = cap.Uc.z();
  const Real cosD = cap.CosDelta, sinD = cap.SinDelta;
  const Real ah2 = ah*ah, az2 = az*az;
  const Real d1 = ah2 - az2;   // AhGeAz==true  branch coefficient.
  const Real d2 = az2 - ah2;   // AhGeAz==false branch coefficient.

  for (Index i = lo; i < hi; i++) {

    Real d = ux*mx[i] + uy*my[i] + uz*mz[i];
    d = (d > 1.0)  ? (Real) 1.0  : d;
    d = (d < -1.0) ? (Real) -1.0 : d;
    const Real c = (mr[i] > 0.0) ? d : (Real) 1.0;

    const Real omc = 1.0 - c*c;
    const Real s = sqrt((omc > 0.0) ? omc : (Real) 0.0);

    const Real tlo = c*cosD + s*sinD;
    const Real thi = c*cosD - s*sinD;
    const Real cos_alpha_lo = (c >= cosD)  ? (Real) 1.0  : tlo;
    const Real cos_alpha_hi = (c <= -cosD) ? (Real) -1.0 : thi;

    const Real s_lo = mr[i]*cos_alpha_hi - rho[i];
    const Real s_hi = mr[i]*cos_alpha_lo + rho[i];
    const Real slo2 = s_lo*s_lo, shi2 = s_hi*s_hi;
    const Real s2_hi = (slo2 < shi2) ? shi2 : slo2;

    Real Q_lo;
    if (AhGeAz) {
      const Real mn = (m2hi[i] < s2_hi) ? m2hi[i] : s2_hi;
      const Real qa = ah2*m2lo[i] - d1*mn;
      const Real qb = az2*m2lo[i];
      Q_lo = (qa < qb) ? qb : qa;
    } else {
      const Real s2_lo = ((s_lo <= 0.0) & (0.0 <= s_hi)) ?
                         (Real) 0.0 : ((shi2 < slo2) ? shi2 : slo2);
      Q_lo = ah2*m2lo[i] + d2*s2_lo;
    }

    out[i-lo] = numer[i] / (1.0 + Q_lo);
  }

}//
//

#undef R3D_EXACT_FP_ATTR

// Explicit instantiation: RatioRun<true>/<false>'s only compiled copies
// live here, in aniso_axi.o -- both 'main' and 'anisotropy_tests' link
// that one object file (Makefile), never recompiling this template from
// the header, so the vectorized loop above can't end up compiled
// differently (e.g. under different ambient flags) between the
// production call site (BuildBinTable(), below) and
// tests/anisotropy_tests.cpp's test_axi_ratio_run_matches_scalar -- the
// same cross-translation-unit-divergence risk parallel.hpp's own
// IMPLEMENTATION NOTE documents for body(i).
template void AxiProposal::RatioRun<true>(
  const RepSoA &, Index, Index, const BinCap &, Real, Real, Real *);
template void AxiProposal::RatioRun<false>(
  const RepSoA &, Index, Index, const BinCap &, Real, Real, Real *);


//////
// METHOD:  AxiProposal :: BuildBinTable()                     [const]
//
//   Builds one incidence bin's cumulative proposal-weight table (and
//   its majorant rate) for the two conversions sharing this bin's
//   u_c/delta ('reps0'/'reps1', 'S0'/'S1' -- see aniso_axi.hpp).
//
void AxiProposal::BuildBinTable(const R3::XYZ & u_c, Real delta,
                                const RepSoA & reps0,
                                const RepSoA & reps1,
                                Index repsPerCell,
                                const std::vector<Real> & S0,
                                const std::vector<Real> & S1,
                                BinTable & table, Real & rate) const {

  const BinCap cap = MakeBinCap(u_c, delta);

  // ratio[j] = BoundRatio(rep_i, cap, mAh, mAz) for a cell's own
  // repsPerCell representatives (i = c*repsPerCell+j), computed by one
  // RatioRun() call (vectorized) over just that cell's contiguous block
  // -- not the whole group -- so this function's own scratch is
  // O(repsPerCell) (a few hundred bytes at degree 9's repsPerCell==64),
  // not O(numReps) (~0.63 MB). This runs from inside a 450-bin
  // ParallelFor (the constructor, below), so up to 'threads' of these
  // can be live at once; an O(numReps) buffer (an earlier version of
  // this function had one) would scale total temporary memory with
  // thread count and, at this class's own ~32 MB/scatterer construction
  // budget and the design doc's 24-thread reference configuration,
  // pushed it over (caught in review, not by any test). Each cell's max
  // is taken over its own small buffer immediately after that cell's
  // RatioRun() call. Max is exact (never rounds: it always returns one
  // of its inputs' own bit patterns) and doesn't care about grouping or
  // order, so processing one cell at a time like this is bit-identical
  // to the original "update both cells' running max inline, one
  // representative at a time" loop, for both groups.
  std::vector<Real> ratio(repsPerCell);

  std::vector<Real> maxRatio0(mNumCells, 0.0), maxRatio1(mNumCells, 0.0);
  const RepSoA      * const groups[2]    = {&reps0, &reps1};
  std::vector<Real> * const maxRatios[2] = {&maxRatio0, &maxRatio1};

  for (int g = 0; g < 2; g++) {
    const RepSoA & reps = *groups[g];
    std::vector<Real> & maxR = *maxRatios[g];
    for (Index c = 0; c < mNumCells; c++) {
      const Index lo = c*repsPerCell, hi = lo+repsPerCell;
      if (mAh >= mAz) {
        RatioRun<true>(reps, lo, hi, cap, mAh, mAz, &ratio[0]);
      } else {
        RatioRun<false>(reps, lo, hi, cap, mAh, mAz, &ratio[0]);
      }
      Real m = maxR[c];
      for (Index j = 0; j < repsPerCell; j++) {
        if (ratio[j] > m) {m = ratio[j];}
      }
      maxR[c] = m;
    }
  }

  // maxRatio[c] >= 1.0 sets the cell's bound to M outright, skipping
  // pow(): exact, since pow(x>=1, kappa+1.5>0) >= 1 makes
  // min(M, M*pow(x,kappa+1.5)) equal M either way -- this just avoids
  // computing (and rounding) a pow() call whose result is thrown away.
  std::vector<Real> bound0(mNumCells), bound1(mNumCells);
  Real sumBinBoundS = 0.0;
  for (Index c = 0; c < mNumCells; c++) {
    bound0[c] = (maxRatio0[c] >= 1.0) ?
      mM : std::min(mM, mM * pow(maxRatio0[c], mKappa+1.5));
    bound1[c] = (maxRatio1[c] >= 1.0) ?
      mM : std::min(mM, mM * pow(maxRatio1[c], mKappa+1.5));
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
// METHOD:  AxiProposal :: Locate()                             [const]
//
//   Folds the local-vertical direction 'u' into a canonical incidence
//   bin (theta, chi for S; theta only for P), and records the mirror
//   flags (Mx, My) needed to map a table cell back to the actual cell.
//   Trig-free (design decision 4): no acos, atan2, or sin call. This
//   folds exactly as the former (now test-only) Fold() did -- compared
//   directly against a frozen transcription of it,
//   tests/anisotropy_tests.cpp's ReferenceLocate3db253a, under ruling
//   R-P1.
//
//   Theta: after the z>=0 fold, u.z() is compared directly against the
//   precomputed edge cosines (ThetaEdgeCosines()) by binary search: bin
//   b is the one with ThetaEdgeCosines()[b+1] < u.z() <=
//   ThetaEdgeCosines()[b], the same interval the old acos(u.z())-based
//   floor(theta_deg/2) assigned (acos is strictly decreasing on
//   [-1,1]), clamped to [0,kThetaBins-1] for free by the search range
//   (no out-of-range u.z() clamp needed, unlike the old acos-domain
//   clamp).
//
//   Chi (S only): w = |u.x()|, v = |u.y()| (post z-fold; the mx/my signs
//   already record which quadrant u.x()/u.y() actually came from). Since
//   chi = atan2(v,w) in [0,90deg] once folded, and tan is increasing
//   there, v*cos(10c deg) >= w*sin(10c deg) iff chi >= 10c deg (dividing
//   both sides by w*cos(10c deg), when positive); chi's bin is the
//   largest c satisfying it, found by binary search over
//   ChiEdgeCosines()/ChiEdgeSines(). w^2+v^2 == sin(theta)^2 exactly (u
//   is a unit vector), so 'w^2+v^2 <= 1e-24' is exactly the old
//   'sinTheta > 1e-12' guard, without the sqrt: near either pole, chi is
//   undefined, and the old code's default (chi=0, i.e. chiBin=0) is used
//   directly.
//
AxiProposal::BinRef AxiProposal::Locate(raytype intype,
                                        const R3::XYZ & u_rel) const {

  R3::XYZ u = u_rel;
  if (u.z() < 0.0) {u = u.Negative();}

  const bool mx = (u.x() > 0.0);
  const Real w  = std::fabs(u.x());   // == -(folded u.x()), >= 0.

  bool my = false;
  Real v = 0.0;
  if (intype == RAY_S) {
    my = (u.y() < 0.0);
    v  = std::fabs(u.y());            // == folded u.y(), >= 0.
  }

  const std::vector<Real> & thetaEdges = ThetaEdgeCosines();
  const Real uz = u.z();
  int thetaLo = 0, thetaHi = kThetaBins - 1;
  while (thetaLo < thetaHi) {
    const int mid = (thetaLo + thetaHi) / 2;
    if (thetaEdges[mid+1] < uz) {thetaHi = mid;} else {thetaLo = mid+1;}
  }

  int chiBin = 0;
  if (intype == RAY_S && (w*w + v*v) > 1e-24) {
    const std::vector<Real> & chiCos = ChiEdgeCosines();
    const std::vector<Real> & chiSin = ChiEdgeSines();
    int chiLo = 0, chiHi = kChiBins - 1;
    while (chiLo < chiHi) {
      const int mid = (chiLo + chiHi + 1) / 2;   // round up: searching
                                                 // for the largest c
                                                 // satisfying the test.
      if (v*chiCos[mid] >= w*chiSin[mid]) {chiLo = mid;} else {chiHi = mid-1;}
    }
    chiBin = chiLo;
  }

  BinRef r;
  r.ThetaBin = thetaLo;
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
  return BinRate(intype, Locate(intype, u_rel));
}

Real AxiProposal::BinRate(raytype intype, const BinRef & bin) const {
  return BinRateFor(intype, bin.ThetaBin, bin.ChiBin);
}


//////
// METHOD:  AxiProposal :: DrawProposal()                      [const]
//
AxiProposal::Draw
AxiProposal::DrawProposal(raytype intype, const R3::XYZ & u_rel) const {
  return DrawProposal(intype, Locate(intype, u_rel));
}

AxiProposal::Draw
AxiProposal::DrawProposal(raytype intype, const BinRef & bin) const {

  const BinTable & table = TableFor(intype, bin.ThetaBin, bin.ChiBin);

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

  const Index actualCell = MapCell(tableCell, bin.Mx, bin.My);
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

  const BinRef f = Locate(intype, u_rel);
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
