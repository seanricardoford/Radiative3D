// scatterers.cpp
//
#include <algorithm>    /* std::max */
#include <cmath>        /* log(), sqrt() */
#include <iomanip>      /* setw() */
#include <cstdlib>
#include <cstddef>      /* std::size_t */
#include <sstream>      /* std::ostringstream */
#include <stdexcept>    /* std::logic_error */
#include "random.hpp"
#include "parallel.hpp"
#include "scatterers.hpp"
#include "phonons.hpp"
#include "ecs.hpp"      /* Global ECS object: GetEast/GetNorth/GetUp */

                        /* Assume srand() has been called elsewhere */

using namespace std;

// in this file:
// CLASS IMPLEMENTATIONS FOR:
//
//   o  Class Scatterer
//


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  Scatterer                                           ****
// ****                                                              ****
//

//
// Static Member Initialization:  (Scatterer Class)
//

Scatterer * Scatterer::cm_ll_first = 0; // Linked list pointers
Scatterer * Scatterer::cm_ll_last  = 0; //
Count Scatterer::cm_ll_count = 0;  // Linked list size
bool Scatterer::cm_MFPOverride_b = false;
Real Scatterer::cm_MFPOverrides[] = {0,0};
bool Scatterer::cm_NoDeflect_b = false;
unsigned Scatterer::cm_BuildThreads = 1;

const raytype Scatterer::kOutTypes[4] = {RAY_P, RAY_S, RAY_P, RAY_S};
                                // Maps conversion types (GPP, GPS, GSP,
                                // GSS) to output raytypes (RAY_P or
                                // RAY_S).

//////
// STATIC METHOD:  Scatterer :: GetScattererMatchingParams()
//
//   Return pointer to an existing Scatterer if one can be found that
//   suitable matches the provided ScatterParams object, otherwise
//   create a new Scatterer and return pointer if no such match can be
//   found.
//
Scatterer *
Scatterer::GetScattererMatchingParams(ScatterParams par) {

  if (cm_MFPOverride_b && cm_NoDeflect_b) { // Check for overrides
    par = ScatterParams(Elastic::HSneak(1.0,0.0,1.0,1.0),1.0,1.0);
    // Squash parameters to dummy value, so that one and only one
    // object will ever be built.
  } //

  Scatterer * result = 0;
  Scatterer * fromlist = cm_ll_first;

  while (fromlist != 0) { // First scan linked list of existing
                          // Scatterers for possible match:
                          //
    Real diff = par.CompareRoughly(fromlist->mParams);
                          // A very rough comparison metric...
    if (diff <= 0) {      // A looser cuttoff may be in order...
      result = fromlist;
      break;
    } else {              // else keep walking the list
      fromlist = fromlist->mpllNext;
    }
  }

  if (result != 0) {      // If we found a match, return it
    return result;        //
  }                       // else continue as follows:

  result = new Scatterer(par);      // Allocate new scatterer
  cm_ll_count++;                    //

  if (cm_ll_first == 0) {           // And append to linked list...
    cm_ll_first = result;           // (set pointer to first if needed)
  }

  result->mpllPrev = cm_ll_last;    // link new object to previous oject;
  result->mpllNext = 0;             // mark new object as last in list;
  cm_ll_last = result;              // update pointer to last;

  if (result->mpllPrev != 0) {      // and link previous object (if any)
    result->mpllPrev->mpllNext = result;  // to new object;
  }

  return result;

}


//////
// CONSTRUCTOR:  Scatterer()
//
Scatterer::Scatterer(ScatterParams par, sampler_e choice) :
  PhononSource(NUM_INTYPES, NUM_OUTTYPES),
                        // Init base class for two input raytypes and
                        // four output raytypes.
  mParams(par),         // Record params for later use if needed
  mAnisotropic( (choice != SAMPLER_AUTO) || !par.IsIsotropic() ),
                        // AUTO defers to the spec: isotropic specs
                        // keep today's exact path. Any explicit,
                        // non-AUTO choice forces an anisotropic
                        // sampling path even on an isotropic spec.
  mUseAxi( choice == SAMPLER_AXISYMMETRIC ||
          (choice == SAMPLER_AUTO && par.IsAxisymmetric()) )
                        // AUTO selects the axisymmetric sampler for
                        // axisymmetric specs; SAMPLER_AXISYMMETRIC
                        // forces it; SAMPLER_GENERAL never selects it.
{

  if (choice == SAMPLER_AXISYMMETRIC && par.GetAx() != par.GetAy()) {
    // AxiProposal always builds with ah=GetAx(), silently ignoring ay.
    // On a spec where ax!=ay that bound is invalid for the actual
    // spectrum, which would later surface as a bound-violation throw in
    // AttemptAxisymmetricScatter() (or worse, a wrong-but-undetected
    // bound). SAMPLER_AXISYMMETRIC is test-only API (production always
    // goes through SAMPLER_AUTO, which never selects it for a
    // non-axisymmetric spec), so reject the misuse here instead.
    throw std::logic_error(
      "SAMPLER_AXISYMMETRIC requires ax == ay");
  }

  if (cm_MFPOverride_b && mAnisotropic && mParams.GetEps() == 0.0) {
    // Controller ruling R11: with eps==0, every GSATO()/mPDists value
    // (and, for the axisymmetric sampler, every AxiProposal bin weight)
    // is identically 0 -- see ScatterParams::PSATO()'s eps^2 prefactor
    // -- so --overridemfp's skip-the-thinning-test loop (both samplers;
    // see AttemptAnisotropicScatter()/AttemptAxisymmetricScatter())
    // would otherwise either compute 0/0 (axisymmetric: NaN never
    // satisfies p>1+kBoundMargin OR RandomUniform()<p, hanging forever)
    // or spin on a degenerate, meaningless all-mass-at-index-0
    // distribution (general). Neither is useful, so fail loudly at
    // construction instead of silently hanging or scattering nonsense.
    throw std::logic_error(
      "--overridemfp requires eps > 0 in anisotropic regions");
  }

  PopulateProbDists(par);   // Populate probability distributions
  PopulateWholeProbs();     //

  if (mAnisotropic && mUseAxi) {   // Build the axisymmetric sampler's
                                   // proposal tables. mPDists is
                                   // already in integrated form (every
                                   // entry's GetMagnitude() was called
                                   // just above, in PopulateWhole-
                                   // Probs()), as AxiProposal's
                                   // GetRangeMagnitude()/GetRandom-
                                   // IndexInRange() calls require.
    mpAxi.reset(new AxiProposal(par, *pTOA, mPDists, GetBuildThreads()));
  }

  if (!cm_MFPOverride_b) {  // Populate Mean Free Path values
    if (mAnisotropic && mUseAxi) {
      mMeanFreeP[RAY_P] = 1.0 / mpAxi->MajorantRate(RAY_P);
      mMeanFreeP[RAY_S] = 1.0 / mpAxi->MajorantRate(RAY_S);
                            // (+inf when the bin majorants are all 0,
                            // e.g. eps==0 -- ordinary IEEE division,
                            // no special-casing needed.)
    } else {
      ComputeMFPs();          //
      if (mAnisotropic) {     // Anisotropic (general-sampler) scatterers
                              // use the isotropic (a=a_min) MFP as the
                              // basis for the majorant proposal
                              // distribution, shortened by the factor
                              // M = ax*ay*az/a_min^3 by which the
                              // anisotropic spectrum can exceed the
                              // isotropic one at the same |m| (see
                              // ScatterParams::AnisoMajorantFactor()).
        Real M = mParams.AnisoMajorantFactor();
        mMeanFreeP[RAY_P] /= M;
        mMeanFreeP[RAY_S] /= M;
      }
    }
  } else {
    mMeanFreeP[RAY_P] = cm_MFPOverrides[RAY_P];   // Override values
    mMeanFreeP[RAY_S] = cm_MFPOverrides[RAY_S];   // if requested
  }

  if (!cm_NoDeflect_b) {    // Populate Dipole Moments
    ComputeDipoles();       //
  } else {
    mDipoles[RAY_P] = 1.0;  // Note: these are informational values, over-
    mDipoles[RAY_S] = 1.0;  // riding them does not prevent deflection. A
  }                         // check in GRSPh() does that.

  if (mAnisotropic) {       // Populate the reporting-only directional
    ComputeReportStats();   // MFPs/DMs that PrintStats() shows for an
  }                         // anisotropic Scatterer. Independent of the
                            // override flags above: these values
                            // describe the physical scattering pattern,
                            // not the (possibly overridden) sampling
                            // behavior.

  //
}//
//


//////
// METHOD:  Scatterer :: PopulateProbDists()         [Helper Function]
//
//   Computes and populates the mPDists members.  Called by the
//   constructor, and whenever the (memory-hungry) mPDists member
//   might need to be re-computed (e.g. after a cache-miss if the
//   distribution has been flushed).
//
//   GSATO() is evaluated over contiguous direction CHUNKS via
//   ParallelFor() (parallel.hpp): a small, GetBuildThreads()-sized index
//   range, one chunk per claim, each processed by an ordinary serial
//   for-loop over its own [chunkBegin,chunkEnd) range of take-off-angle
//   indices. This -- not a ParallelFor spanning all nTOA directions
//   one-at-a-time -- is deliberate, found necessary by profiling this
//   task's first (per-direction) attempt at degree 9 (5,242,880
//   directions): claiming one shared atomic<Index> counter 5.2 million
//   times, split across ~24 threads, spends most of the run's CPU time
//   (confirmed with `sample`) contending for that one cache line, not
//   computing GSATO() -- and it barely changed wall-clock time at all
//   versus threads=1, despite burning roughly 15x the CPU-seconds.
//   Chunking reduces the claim count to GetBuildThreads(), eliminating
//   that contention, and (as a side effect) confines false-sharing on
//   mPDists[*]/m_spol writes to the single element at each chunk
//   boundary instead of every element.
//
//   Each worker writes mPDists[*][k] and m_spol[k] directly for the k's
//   in its own chunk -- not into a separate temporary array -- because a
//   per-direction temporary here (4 Reals x nTOA) would itself cost
//   ~167 MB at degree 9, far past the ~32 MB/scatterer temporary-build-
//   memory budget (docs/superpowers/specs/2026-10-16-anisotropic-scale-
//   lengths-design.md). This is safe: probability.hpp's
//   SetRelativeProb() lazily resizes/allocates its shared mDist vector
//   (and flips mbAllocated) on a distribution's *first* call, which
//   would race if that first call came from more than one thread at
//   once -- so each of the four distributions is allocated by one
//   dummy, serial SetRelativeProb(0, ...) call below, *before*
//   ParallelFor spawns any workers (std::thread's constructor
//   establishes the needed happens-before relationship). From then on,
//   SetRelativeProb(k, ·) only ever assigns mDist[k] -- a plain,
//   non-cumulative store to an already-allocated std::vector<Real>,
//   never touching any other element or shared bookkeeping -- so
//   concurrent calls from different chunks (which never touch the same
//   k) are data-race-free. m_spol is likewise already resized (below)
//   before the parallel section starts.
//
//   Because every (distribution, k) pair is written at most once, by
//   whichever chunk's range contains k, with a value computed by a
//   single, deterministic, side-effect-free GSATO() call, the final
//   contents of mPDists[*]/m_spol are bit-identical to a fully serial
//   build for any thread count or chunking -- see docs/superpowers/
//   specs/2026-10-17-axisymmetric-sampler-performance-design.md,
//   decision 8.
//
//   AREA-WEIGHTED QUADRATURE (Task 3): Get_TOA_Quadrature() is read once,
//   into 'areaMode' below, before ParallelFor spawns any worker -- not
//   re-read per chunk or per direction. In TOA_AREA mode, each of the
//   four GSATO() values (never m_spol[k]; see Scatterer::IntegrateRows()'s
//   identical rule and the Global Constraints) is multiplied by
//   TOAWeight(k) before SetRelativeProb(). The TOA_EQUAL branch is its
//   own textually separate loop, storing exactly the values it always
//   has: this file is compiled with GCC's default -ffp-contract=fast,
//   which could otherwise fuse a weight multiply -- even one behind an
//   'if' that's always false in TOA_EQUAL mode, even one that would
//   evaluate to exactly 1.0 -- into a different FMA than before, shifting
//   the default mode's output by a rounding ulp and breaking byte-
//   identity with commit 319c8dd (see IntegrateRows()'s own precaution in
//   this file). Both branches still write every (distribution, k) pair
//   at most once from a single chunk, so the thread-invariance argument
//   above holds in either mode.
//
void Scatterer::PopulateProbDists(ScatterParams par) {

  S2::S2Set & toa = (*pTOA);    // (Convenience alias)

  m_spol.clear();       // spol is an array of polarization angles for
  m_spol.resize(nTOA);  // the S->S conversionevents.  Technically,
                        // there are probabilities of going to either
                        // SH (spol=90) or SV (spol=0), but we only
                        // treat the probability of going to S
                        // generally, and code a polarization angle
                        // that represents an appropriate linear
                        // combination of 0 and 90, (based on the
                        // individual SV/SH probabilities, which are
                        // known to the GSATO function.)

  if (nTOA > 0) {   // Force each distribution's one-time lazy allocation
                    // now, serially -- see header comment. (Guards the
                    // degenerate nTOA==0 case, which SetRelativeProb(0,
                    // ...) itself could not: PopulateWholeProbs()'s
                    // later GetMagnitude() call already turns that case
                    // into today's controlled AckSpit() failure, same as
                    // before this change.)
    mPDists[GPP].SetRelativeProb(0, 0.0);  // Dummy value: whichever
    mPDists[GPS].SetRelativeProb(0, 0.0);  // chunk's range covers k==0
    mPDists[GSP].SetRelativeProb(0, 0.0);  // below overwrites it with
    mPDists[GSS].SetRelativeProb(0, 0.0);  // the real GSATO() result.
  }

  //
  // RADIATION PATTERNS:
  //   (Populate PF arrays)
  const Index n = (Index) nTOA;
  const unsigned threads = GetBuildThreads();
  const Index nChunks = (threads == 0) ? 1 : (Index) threads;
  const bool areaMode = (Get_TOA_Quadrature() == TOA_AREA);  // Read once;
                              // see header comment above.

  ParallelFor(0, nChunks, threads, [&](Index c) {
    const Index chunkBegin = (Index) (((BigIndex) n * c)       / nChunks);
    const Index chunkEnd   = (Index) (((BigIndex) n * (c + 1)) / nChunks);
    if (areaMode) {
      for (Index k = chunkBegin; k < chunkEnd; k++) {
        Real gpp, gps, gsp, gss, spolv;
        par.GSATO(toa[k], gpp, gps, gsp, gss, spolv); // Get g-values and
        const Real w = TOAWeight(k);                  // spol value
        mPDists[GPP].SetRelativeProb(k, w * gpp);
        mPDists[GPS].SetRelativeProb(k, w * gps);
        mPDists[GSP].SetRelativeProb(k, w * gsp);
        mPDists[GSS].SetRelativeProb(k, w * gss);
        m_spol[k]  = spolv;        // Never weighted.
      }
    } else {
      for (Index k = chunkBegin; k < chunkEnd; k++) {
        Real gpp, gps, gsp, gss, spolv;
        par.GSATO(toa[k], gpp, gps, gsp, gss, spolv); // Get g-values and
        mPDists[GPP].SetRelativeProb(k, gpp);         // spol value
        mPDists[GPS].SetRelativeProb(k, gps);
        mPDists[GSP].SetRelativeProb(k, gsp);
        mPDists[GSS].SetRelativeProb(k, gss);
        m_spol[k]  = spolv;
      }
    }
  });

}//
//


//////
// METHOD:  Scatterer :: PopulateWholeProbs()        [Helper Function]
//
//   Computes and populates the WholeProbs probability arrays.
//
//   Assumes mPDists[] has already been populated.
//
void Scatterer::PopulateWholeProbs() {

  mWholeProbs[IN_P].SetRelativeProb(GPP, mPDists[GPP].GetMagnitude());
  mWholeProbs[IN_P].SetRelativeProb(GPS, mPDists[GPS].GetMagnitude());
  mWholeProbs[IN_P].SetRelativeProb(GSP, 0); // (S->P precluded for IN_P)
  mWholeProbs[IN_P].SetRelativeProb(GSS, 0); // (S->S precluded for IN_P)

  mWholeProbs[IN_S].SetRelativeProb(GPP, 0); // (P->P precluded for IN_S)
  mWholeProbs[IN_S].SetRelativeProb(GPS, 0); // (P->S precluded for IN_S)
  mWholeProbs[IN_S].SetRelativeProb(GSP, mPDists[GSP].GetMagnitude());
  mWholeProbs[IN_S].SetRelativeProb(GSS, mPDists[GSS].GetMagnitude());

}//
//


//////
// METHOD:  Scatterer :: ComputeMFPs()               [Helper Function]
//
//   Computes and stores Mean Free Paths.
//
//   Assumes mWholeProbs[] has already been populated.
//
void Scatterer::ComputeMFPs() {

  //
  // Get the *inverse* MFPs  (Basically the G0 values)
  //
  Real imfp_p = mWholeProbs[IN_P].GetMagnitude();
  Real imfp_s = mWholeProbs[IN_S].GetMagnitude();
  //
  // Compute mean-free-paths:
  //
  //   (The imfp's currently contain the sums of the G-functions -- each
  //   term already weighted by TOAWeight(k) when Get_TOA_Quadrature() ==
  //   TOA_AREA, or left unweighted (today's TOA_EQUAL behavior)
  //   otherwise; see PopulateProbDists(). The "probability of scattering
  //   per unit length" is the "surface average" of the G-functions, so
  //   instead of dividing by 4*Pi we divide by nTOA; that stays correct
  //   in EITHER mode, because Sum_k TOAWeight(k) == nTOA always --
  //   TOA_EQUAL sums nTOA terms of 1.0, and TOA_AREA's terms
  //   SolidAngles()[k]*(nTOA/(4*Pi)) sum to nTOA precisely because the
  //   solid angles themselves sum to 4*Pi (tessel_solid_angles,
  //   tests/anisotropy_tests.cpp). See docs/toa-quadrature-bias.md for
  //   the bias this quadrature choice resolves.)
  //
  imfp_p /= nTOA;
  imfp_s /= nTOA;
  mMeanFreeP[IN_P] = 1.0 / imfp_p;
  mMeanFreeP[IN_S] = 1.0 / imfp_s;

}//
//


//////
// METHOD:   Scatterer :: ComputeDipoles()
//
//   Computes the dipole moments of the scattering shapes. This
//   provides a useful way to characterize the scattering shapes on a
//   spectrum between "forward scattering" and "back scattering."
//
//   Returns a real-valued number between -1.0 and 1.0 where numbers
//   approaching 1.0 mean dominantly forward-scattering (minimal path
//   deviation), -1.0 means dominant back-scattering (~180*
//   reversals), and 0.0 means either equal contribution from both
//   forward and backwards or else the dominant scattering direction
//   is side-deflection.
//
//   Computation is achieved by multiplying each differential
//   probability by cos(theta) where theta is the deflection angle
//   (co-lattitude from forward axis).
//
//   Assumes mPDists[] and mWholeProbs[] have already been computed
//   and are available for analysis.
//
//   Inherits whichever quadrature PopulateProbDists() applied to
//   mPDists (TOA_EQUAL or TOA_AREA) with no code change here:
//   GetDiffProb() normalizes by GetMagnitude(), so each moments[*] term
//   and the magnitude it's eventually scaled by (mWholeProbs[*]'s own
//   GetDiffProb(), computed the same way from PopulateWholeProbs()'s
//   already-weighted mPDists[*].GetMagnitude() sums) carry the same
//   per-direction TOAWeight() factors throughout -- an area-mode dipole
//   is therefore the properly area-weighted moment, not the equal-
//   weighted approximation.
//
void Scatterer::ComputeDipoles() {

  const S2::S2Set & toa = (*pTOA);          // Alias (for convenience)
  Real moments[NUM_OUTTYPES] = {0,0,0,0};   // Four moments

  for (int idx = 0; idx < nTOA; idx++) {
    Real costh = toa[idx].z();  // (Unit vec, so z is cos theta)
    moments[GPP] += costh * mPDists[GPP].GetDiffProb(idx);
    moments[GPS] += costh * mPDists[GPS].GetDiffProb(idx);
    moments[GSP] += costh * mPDists[GSP].GetDiffProb(idx);
    moments[GSS] += costh * mPDists[GSS].GetDiffProb(idx);
  }

  mDipoles[IN_P] = moments[GPP] * mWholeProbs[IN_P].GetDiffProb(GPP)
                 + moments[GPS] * mWholeProbs[IN_P].GetDiffProb(GPS);
  mDipoles[IN_S] = moments[GSP] * mWholeProbs[IN_S].GetDiffProb(GSP)
                 + moments[GSS] * mWholeProbs[IN_S].GetDiffProb(GSS);

}//
//


//////
// METHOD:  Scatterer :: SetMeanFreePathsPS()
//
//   Though the MFPs are computed automatically by the constructor,
//   this method can be used to change them after the fact.  This is
//   mainly for diagnostic and development/debugging purposes, or for
//   creating very specialized and not-necessarily-physical test
//   models. It should not be part of a normal grid-based model
//   building process.
//
//   NOTE: C.f. this function with class static OverrideMFP().  This function
//   overrides the MFP values for a particular object.  OverrideMFP()
//   overrides all objects at construction time.
//
void Scatterer::SetMeanFreePathsPS(Real mfpP, Real mfpS) {
  mMeanFreeP[RAY_P] = mfpP;
  mMeanFreeP[RAY_S] = mfpS;
}//
//


//////
// METHOD:  Scatterer :: GetRandomPathLength()
//
// PURPOSE: Return a randomly-generated path length to a scattering
//          event.  Lengths are generated in such a way that the
//          average, or "mean free path," of many iterations of this
//          function will tend towards the values stored in
//          MeanFreePath[].
//
Real Scatterer::GetRandomPathLength(raytype intype) {

  Real r = RandomUniform();
                // A double in the range [0,1)
  r = 1.0 - r;  // A double in the range (0,1]
                // TODO: This is a little ugly. Maybe compute r like
                // we do in rtcoef.cpp instead. (~CJS 20140514)

  return -log(r) * mMeanFreeP[intype];

}


//////
// METHOD:  Scatterer :: GetRandomScatteredRelativePhonon()
//
// PURPOSE: Generate and return a "relative" phonon at random from the
//          scattering probability distributions.  The "relative"
//          phonon can be used to "bend and rotate" the path of the
//          incoming phonon.
//
Phonon Scatterer::GetRandomScatteredRelativePhonon(raytype intype) {
  // Guard against anisotropic scattering via the isotropic path.
  // Anisotropic scatterers must use AttemptAnisotropicScatter instead.
  if (!mParams.IsIsotropic()) {
    throw std::logic_error("anisotropic scatterers must be sampled through AttemptAnisotropicScatter");
  }

  S2::S2Set & toa = (*pTOA);             // Alias

  if (cm_NoDeflect_b) {     // First, check if deflection is overridden:
    Phonon nodeflect = Phonon(S2::ThetaPhi(0,0),intype);
    nodeflect.SetPolarization(0); // And return a non-deflected
    return nodeflect;             // result if so.
  }                               // Otherwise, continue:

  // Get output conversion type and corresponding output ray type:
  out_types_e conv = (out_types_e) mWholeProbs[intype].GetRandomIndex();
  raytype ort = kOutTypes[conv];

  // Get take-off-angle index for output ray:
  Index toa_index = mPDists[conv].GetRandomIndex();

  // Ray polariztion depends on conversion type:
  Real pol = 0;
  switch (conv) {
  case GPP : // Polarization meaningless for P-waves; default to zero.
    pol = 0;
    break;
  case GSP :
    pol = 0;
    break;
  case GPS : // This conversion always makes theta-hat polarization.
    pol = 0; // (0 for theta-hat, Pi/2 for phi-hat)
    break;
  case GSS : // This conversion looks up polarization in spol array.
    pol = m_spol[toa_index];
    break;
  default:   // (Shouldn't end up here - this is just a #pragma
    break;   //  for NUM_OUTTYPES - avoids compiler warning.)
  }

  // Generate the Phonon:
  Phonon phon = Phonon(toa[toa_index], ort);
  phon.SetPolarization(pol);

  return phon;

}


//////
// METHOD:  Scatterer :: ToFabric()                            [static]
//
// PURPOSE: Express vector 'v' in the local fabric (East, North, Up)
//          frame at 'loc', using the global ECS object to determine
//          the local East/North/Up directions. The fabric frame is
//          what ScatterParams::AnisoRatio()/AnisoAcceptance() expect
//          their wavevector argument to be expressed in.
//
R3::XYZ Scatterer::ToFabric(const R3::XYZ & loc, const R3::XYZ & v) {
  return R3::XYZ(v.Dot(ECS.GetEast(loc)),
                 v.Dot(ECS.GetNorth(loc)),
                 v.Dot(ECS.GetUp(loc)));
}


//////
// METHOD:  Scatterer :: RotationForS1()                       [static]
//
// PURPOSE: Return the rotation angle 'rot' such that
//          R3::OrthoAxes(dir.Theta(),dir.Phi(),rot).S1() is parallel to
//          the component of 's1' normal to 'dir' (i.e., with any
//          component of 's1' along 'dir' removed). If that normal
//          component is (numerically) zero -- e.g. because 's1' is
//          parallel to 'dir', as happens when a caller has no
//          meaningful polarization direction to offer (incoming
//          P-waves) -- returns 0, matching OrthoAxes' own convention
//          that rot=0 aligns S1 with E1 (theta-hat).
//
Real Scatterer::RotationForS1(const S2::ThetaPhi & dir, const R3::XYZ & s1) {

  const R3::XYZ dhat = dir.XYZ();
  const R3::XYZ proj = s1 + dhat.ScaledBy(-s1.Dot(dhat));
                              // s1, with its component along dhat
                              // removed.

  if (proj.MagSquared() < 1.e-18) {
    return 0;                // Degenerate: no meaningful normal
                             // component to align S1 with.
  }

  const R3::OrthoAxes AA0(dir.Theta(), dir.Phi(), 0);
                              // At rot=0, AA0.S1()==AA0.E1() and
                              // AA0.S2()==AA0.E2() (see OrthoAxes'
                              // constructor); E1,E2 don't depend on
                              // rot, so we can read them straight off
                              // this reference frame.
  const R3::XYZ proj_hat = proj.Unit();
  const Real cosrot = proj_hat.Dot(AA0.E1());
  const Real sinrot = proj_hat.Dot(AA0.E2());

  return atan2(sinrot, cosrot);

}


//////
// METHOD:  Scatterer :: CanonicalPPolarization()              [static]
//
// PURPOSE: Returns the polarization-rotation angle that the
//          axisymmetric sampler (SAMPLER_AXISYMMETRIC) uses in place of
//          a real polarization for an incoming P phonon (which has
//          none): the angle that puts the local vertical 'up' into the
//          (S1,S3) plane, with a non-positive S1 component. This is
//          RotationForS1(dir, (up.d)d - up), i.e. S1 is aligned with
//          the normal component of -up (equivalently, of (up.d)d-up,
//          the negative of up's component perpendicular to d). With
//          this convention the local vertical always sits at chi==0
//          (see AxiProposal's fold convention), so the P incidence
//          tables need only a theta dimension.
//
Real Scatterer::CanonicalPPolarization(const S2::ThetaPhi & dir,
                                       const R3::XYZ & up) {
  const R3::XYZ d = dir.XYZ();
  const R3::XYZ s1target = d.ScaledBy(up.Dot(d)) + up.ScaledBy(-1.0);
  return RotationForS1(dir, s1target);
}


//////
// METHOD:  Scatterer :: BuildRelativePhonon()                 [const]
//
// PURPOSE: Builds the "relative" Phonon for an accepted anisotropic-
//          scatter event (general or axisymmetric sampler alike):
//          looks up the output raytype via kOutTypes[conv] and the
//          polarization (m_spol[toa_index] for a GSS conversion, else
//          0), or -- if Scatterer::SetNoDeflect() is in effect -- a
//          non-deflected dummy Phonon of raytype 'intype' instead (the
//          event rate stays anisotropic; only the reported deflection
//          is suppressed). Factored out of AttemptAnisotropicScatter()
//          and AttemptAxisymmetricScatter(), which were otherwise
//          identical here.
//
Phonon Scatterer::BuildRelativePhonon(out_types_e conv, Index toa_index,
                                      raytype intype) const {

  if (cm_NoDeflect_b) {
    Phonon rel(S2::ThetaPhi(0,0), intype);
    rel.SetPolarization(0);
    return rel;
  }

  const raytype ort = kOutTypes[conv];

  Real relpol = 0;
  switch (conv) {
  case GSS:
    relpol = m_spol[toa_index];
    break;
  default:
    relpol = 0;
    break;
  }

  Phonon rel((*pTOA)[toa_index], ort);
  rel.SetPolarization(relpol);
  return rel;

}



//////
// METHOD:  Scatterer :: AttemptAnisotropicScatter()
//
// PURPOSE: Exact (thinning / delta-tracking) sampler for anisotropic
//          (ax,ay,az not all equal) ScatterParams. Proposals are drawn
//          from the same isotropic (a=a_min) tables that
//          GetRandomScatteredRelativePhonon() uses, then accepted with
//          probability ScatterParams::AnisoAcceptance() evaluated at
//          the proposal's exchange wavevector, expressed in the local
//          fabric frame. See docs/superpowers/specs/
//          2026-10-16-anisotropic-scale-lengths-design.md, "Sampling
//          method", for the derivation of why this is exact.
//
// INPUTS:  intype - Incoming ray type (RAY_P or RAY_S).
//          loc    - Location of the tentative collision (used only to
//                   look up the local fabric frame via ECS).
//          dir    - Incoming phonon direction, exactly as
//                   Phonon::Transform() would read it (mDir).
//          pol    - Incoming polarization-rotation angle, exactly as
//                   Phonon::Transform() would read it (mPol).
//
Scatterer::AnisoAttempt
Scatterer::AttemptAnisotropicScatter(raytype intype, const R3::XYZ & loc,
                                     const S2::ThetaPhi & dir, Real pol) {

  if (mUseAxi) {
    return AttemptAxisymmetricScatter(intype, loc, dir, pol);
  }

  const Real k_in = (intype == RAY_P) ? mParams.GetKP() : mParams.GetKS();

  // Incoming P has no physical polarization-rotation angle, but the
  // equal-weight take-off-angle quadrature is not rotation-invariant, so
  // the frame AA is built in determines the realized P rate (see
  // docs/toa-quadrature-bias.md). Use the canonical frame -- the
  // highest-rate orientation, same as AttemptAxisymmetricScatter()'s
  // u_rel closed form -- instead of the carried 'pol', so the two
  // samplers agree for already-scattered P phonons. S incidence keeps
  // the real carried polarization unchanged.
  const Real pol_used = (intype == RAY_P) ?
                        CanonicalPPolarization(dir, ECS.GetUp(loc)) : pol;

  const R3::OrthoAxes AA(dir.Theta(), dir.Phi(), pol_used);  // Same
                              // construction as Phonon::Transform()
                              // would use with pol_used as mPol.
  const R3::XYZ dirXYZ = dir.XYZ();

  out_types_e conv = GPP;
  Index toa_index = 0;
  int draws = 0;
  bool accepted = false;

  do {
    draws++;

    // 1. Draw a proposal (conversion type, then take-off angle) from
    //    the isotropic (a=a_min) tables:
    conv = (out_types_e) mWholeProbs[intype].GetRandomIndex();
    toa_index = mPDists[conv].GetRandomIndex();

    // 2. Compute the exchange wavevector m = k_out*e_out - k_in*dir:
    const Real k_out = (kOutTypes[conv] == RAY_P) ?
                        mParams.GetKP() : mParams.GetKS();
    const R3::XYZ e_out = AA.Express((*pTOA)[toa_index].XYZ());
    const R3::XYZ m = e_out.ScaledBy(k_out) + dirXYZ.ScaledBy(-k_in);

    // 3. Accept/reject using the anisotropic/isotropic spectrum ratio
    //    at 'm', expressed in the local fabric frame:
    const Real accept_prob = mParams.AnisoAcceptance(ToFabric(loc, m));
    accepted = (RandomUniform() < accept_prob);

  } while (!accepted && cm_MFPOverride_b);
                              // Without an MFP override, a single
                              // proposal/accept trial corresponds to
                              // one majorant-distance draw: a
                              // rejection is a null collision, and the
                              // caller draws the next tentative
                              // distance and tries again. With an MFP
                              // override the path length is fixed, so
                              // we must retry internally until a real
                              // event is found.

  AnisoAttempt result;
  result.InPol = pol;        // Non-real attempts carry the incoming pol
                              // as-is, matching AttemptAxisymmetric-
                              // Scatter(): Propagate() reads InPol only
                              // when IsReal (phonons.cpp).
  result.Conv = (int) conv;
  result.Toa = toa_index;
  result.Draws = draws;
  result.IsReal = accepted;

  if (accepted) {
    if (intype == RAY_P) result.InPol = pol_used;  // Canonical frame,
                              // so a later scattering event off this
                              // (now real) P phonon also starts from
                              // the canonical orientation.
    result.Rel = BuildRelativePhonon(conv, toa_index, intype);
  }

  return result;

}


//////
// METHOD:  Scatterer :: AttemptAxisymmetricScatter()         [protected]
//
// PURPOSE: The SAMPLER_AXISYMMETRIC implementation of
//          AttemptAnisotropicScatter() (which dispatches here when
//          mUseAxi): draws proposals from the incidence-tailored
//          AxiProposal tables (mpAxi) instead of the general sampler's
//          flat isotropic tables, then finishes with the same kind of
//          exact accept/reject test. See docs/superpowers/specs/
//          2026-10-16-anisotropic-scale-lengths-design.md, "Axisymmetric
//          sampler".
//
// INPUTS:  Same as AttemptAnisotropicScatter().
//
Scatterer::AnisoAttempt
Scatterer::AttemptAxisymmetricScatter(raytype intype, const R3::XYZ & loc,
                                      const S2::ThetaPhi & dir, Real pol) {

  const R3::XYZ up = ECS.GetUp(loc);
  const R3::XYZ dirXYZ = dir.XYZ();

  // u_rel: the local vertical 'up', expressed in the incoming ray's
  // (S1,S2,S3) frame. For P there is no real polarization to build that
  // frame from; instead of building it via CanonicalPPolarization()
  // (which costs an atan2 through RotationForS1()) just to read off
  // up's components, use the closed form that frame is guaranteed to
  // produce -- see test_axi_p_urel_matches_canonical and
  // docs/superpowers/specs/2026-10-17-axisymmetric-sampler-performance-
  // design.md, decision 4. For S, 'pol' is a real polarization-rotation
  // angle, so the frame is built exactly as before.
  R3::XYZ u_rel;
  if (intype == RAY_P) {
    Real cz = up.Dot(dirXYZ);
    if (cz > 1.0)  cz = 1.0;
    if (cz < -1.0) cz = -1.0;
    const Real sz = sqrt(std::max((Real)0.0, 1.0 - cz*cz));
    u_rel = R3::XYZ(-sz, 0.0, cz);
  } else {
    const R3::OrthoAxes AA(dir.Theta(), dir.Phi(), pol);  // Same
                                // construction as Phonon::Transform().
    u_rel = R3::XYZ(up.Dot(AA.S1()), up.Dot(AA.S2()), up.Dot(AA.S3()));
  }

  const Real k_in = (intype == RAY_P) ? mParams.GetKP() : mParams.GetKS();

  AnisoAttempt result;
  result.InPol = pol;        // Non-real attempts carry the incoming pol
                              // as-is: Propagate() reads InPol only when
                              // IsReal (phonons.cpp). The canonical P
                              // polarization is computed below only for
                              // an accepted event.
  result.Draws = 0;

  const AxiProposal::BinRef bin = mpAxi->Locate(intype, u_rel);  // Once
                              // per tentative collision -- trig-free
                              // (design decision 4), unlike the two
                              // Fold() calls (one per BinRate()/
                              // DrawProposal()) this replaces.

  if (!cm_MFPOverride_b) {
    // 1. Thinning test, costing no table draw (step 4 of the brief):
    const Real majorant = mpAxi->MajorantRate(intype);
    const Real binrate  = mpAxi->BinRate(intype, bin);
    if (RandomUniform() * majorant >= binrate) {
      result.IsReal = false;
      return result;             // Null collision: Draws stays 0.
    }
  }

  out_types_e conv = GPP;
  Index toa_index = 0;
  bool accepted = false;

  do {
    result.Draws++;

    // 2. Draw a (conversion, take-off-angle) proposal from the
    //    incidence-tailored tables:
    const AxiProposal::Draw d = mpAxi->DrawProposal(intype, bin);
    conv = (out_types_e) d.Conv;
    toa_index = d.Toa;

    // 3. |m|^2 and (m.up)^2, straight in the ray frame -- no AA.Express()
    //    or ToFabric() transform needed. With t = the take-off direction
    //    (expressed in the (S1,S2,S3) frame, whose third axis is the
    //    incoming direction) and e_out = AA.Express(t): e_out.Dot(up) ==
    //    t.Dot(u_rel) and e_out.Dot(dirXYZ) == t.z() (AA.S3()==dirXYZ),
    //    since AA is an orthonormal change of basis; |e_out|==|dirXYZ|==1
    //    gives m2 directly. See docs/superpowers/specs/2026-10-17-
    //    axisymmetric-sampler-performance-design.md, decision 4.
    const Real k_out = (kOutTypes[conv] == RAY_P) ?
                        mParams.GetKP() : mParams.GetKS();
    const R3::XYZ t = (*pTOA)[toa_index].XYZ();
    const Real mu = k_out*t.Dot(u_rel) - k_in*u_rel.z();
    const Real m2 = k_out*k_out + k_in*k_in - 2.0*k_out*k_in*t.z();

    // 4. Exact accept/reject using the stored proposal weight and a
    //    freshly-computed cell sum (never the stored one -- see
    //    AxiProposal::Draw's header comment). A non-positive Weight
    //    would make p a 0/0 (or negative) NaN/garbage value that
    //    neither the bound-violation check below nor the acceptance
    //    test can catch (NaN compares false to everything), so fail
    //    loudly instead: DrawProposal() must never actually hand back
    //    such an entry (the Scatterer constructor's eps==0+
    //    --overridemfp guard rules out the one known way to reach an
    //    identically-zero table), so this is a belt-and-braces check.
    if (d.Weight <= 0.0) {
      std::ostringstream oss;
      oss << "AxiProposal drew a non-positive Weight: ax=" << mParams.GetAx()
          << " ay=" << mParams.GetAy() << " az=" << mParams.GetAz()
          << " kappa=" << mParams.GetKappa()
          << " u_rel=(" << u_rel.x() << "," << u_rel.y() << ","
          << u_rel.z() << ")"
          << " conv=" << (int)conv << " toa=" << toa_index
          << " Weight=" << d.Weight;
      throw std::logic_error(oss.str());
    }

    const Real ratio = mParams.AxiRatio(m2, mu*mu);
    const Real p = ratio * d.CellSum / d.Weight;

    if (p > 1.0 + AxiProposal::kBoundMargin) {
      std::ostringstream oss;
      oss << "AxiProposal bound violated: ax=" << mParams.GetAx()
          << " ay=" << mParams.GetAy() << " az=" << mParams.GetAz()
          << " kappa=" << mParams.GetKappa()
          << " u_rel=(" << u_rel.x() << "," << u_rel.y() << ","
          << u_rel.z() << ")"
          << " conv=" << (int)conv << " toa=" << toa_index
          << " p=" << p;
      throw std::logic_error(oss.str());
    }

    accepted = (RandomUniform() < p);

  } while (!accepted && cm_MFPOverride_b);
                              // Without an MFP override, a single
                              // proposal/accept trial corresponds to
                              // one majorant-distance draw (see
                              // AttemptAnisotropicScatter()'s comment
                              // on the general sampler's identical
                              // loop). With an MFP override the path
                              // length is fixed, so we must retry
                              // internally until a real event is found.

  result.Conv = (int) conv;
  result.Toa = toa_index;
  result.IsReal = accepted;

  if (accepted) {
    result.InPol = (intype == RAY_P) ?
                   CanonicalPPolarization(dir, up) : pol;
    result.Rel = BuildRelativePhonon(conv, toa_index, intype);
  }

  return result;

}


//////
// METHOD:  Scatterer :: IntegratedRate()                     [const]
//
// PURPOSE: Exact total scattering rate (1/mean-free-path) for an
//          incoming ray of type 'intype' travelling in direction 'dir'
//          with polarization direction 'pol' (both already expressed
//          in fabric coordinates), used to validate
//          AttemptAnisotropicScatter()'s realized event rate.
//
//          Computed by averaging, over the relative-direction grid
//          'dirs' (nominally the class-level take-off-angle set), the
//          isotropic (a=a_min) GSATO() value times AnisoRatio() at the
//          corresponding exchange wavevector, summed over the two
//          conversions available to 'intype' (GPP+GPS for RAY_P,
//          GSP+GSS for RAY_S). This mirrors ComputeMFPs(), but with
//          each term weighted by the anisotropic/isotropic spectrum
//          ratio instead of assuming it is 1.
//
//          'solidAngles', if non-null, replaces that equal-weighted
//          average with an area-weighted one -- see IntegrateRows() --
//          and must be the same size as 'dirs' or this throws
//          std::logic_error. A thin wrapper around
//          IntegrateRateAndDipole() that discards the dipole.
//
Real Scatterer::IntegratedRate(raytype intype, const R3::XYZ & dir,
                               const R3::XYZ & pol,
                               const S2::S2Set & dirs,
                               const std::vector<Real> * solidAngles) const {

  Real rate, dipole;
  IntegrateRateAndDipole(intype, dir, pol, dirs, rate, dipole, solidAngles);
  return rate;

}


//////
// METHOD:  Scatterer :: BuildReportRows()                     [const]
//
// PURPOSE: Fills 'rows' (resized to dirs.size()) with the GSATO()
//          values and cos(psi) for every direction in 'dirs', via
//          ParallelFor() (parallel.hpp) over direction indices run on
//          up to 'threads' workers. Each row is a pure function of
//          dirs[k] and mParams (ScatterParams::GSATO() touches no
//          mutable state), so 'rows' comes out bit-identical to a
//          serial build regardless of 'threads' or how indices are
//          scheduled across workers -- see docs/superpowers/specs/
//          2026-10-17-axisymmetric-sampler-performance-design.md,
//          decision 2.
//
//          Shared by IntegrateRateAndDipole() (builds rows for its own
//          'dirs' argument) and ComputeReportStats() (builds rows once
//          for ReportDirs(), reused by every one of its 9 or 17
//          states), so GSATO() is never evaluated more than once per
//          (Scatterer,direction) pair.
//
void Scatterer::BuildReportRows(const S2::S2Set & dirs, unsigned threads,
                                std::vector<GRow> & rows) const {

  const Index n = (Index) dirs.size();
  rows.resize(n);

  ParallelFor(0, n, threads, [&](Index k) {
    Real gpp, gps, gsp, gss, spol;
    mParams.GSATO(dirs[k], gpp, gps, gsp, gss, spol);
    rows[k].G[GPP] = gpp;
    rows[k].G[GPS] = gps;
    rows[k].G[GSP] = gsp;
    rows[k].G[GSS] = gss;
    rows[k].CosPsi = dirs[k].z();   // (Unit vec, so z is cos(psi).)
  });

}


//////
// METHOD:  Scatterer :: IntegrateRows()                       [const]
//
// PURPOSE: The row-based integrator shared by IntegrateRateAndDipole()
//          and ComputeReportStats(): computes the same rate/dipole sum
//          IntegrateRateAndDipole()'s header comment describes, but
//          reads each direction's GSATO() values and cos(psi) from the
//          precomputed 'rows' (BuildReportRows()) instead of calling
//          GSATO() itself. The summation stays serial, in k-ascending
//          order -- identical term-for-term to today's single combined
//          loop -- so a call's result never depends on how (or on how
//          many threads) 'rows' was built.
//
//          'solidAngles' selects the quadrature: null (TOA_EQUAL) runs
//          today's loop, completely unchanged -- same expressions, same
//          k-ascending summation order as before this method took a
//          'solidAngles' argument. This file is compiled with GCC's
//          default -ffp-contract=fast, which is free to fuse any
//          multiply immediately followed by an add into a single FMA
//          instruction; inserting a weight factor into that loop, even
//          behind an 'if' that is always false in TOA_EQUAL mode, could
//          still change which multiplications the compiler chooses to
//          fuse and shift the default mode's output by a rounding
//          ulp -- enough to break byte-identity with 319c8dd. The
//          TOA_AREA case below is therefore a separate loop, never the
//          one above, reusing only the shared setup (c1/c2, k_in, AA)
//          both cases need.
//
//          Non-null 'solidAngles' (TOA_AREA, or any other explicit
//          per-direction quadrature over 'dirs'): its size must equal
//          dirs.size() (std::logic_error if not -- a silent out-of-
//          bounds read would otherwise follow), and each direction k's
//          term is scaled by w_k = (*solidAngles)[k] * n/(4*Pi), 'n'
//          here being THIS call's dirs.size() -- not the class-level
//          nTOA, since 'dirs' is often ReportDirs(), a fixed degree-6
//          set independent of the run's --toa-degree:
//
//              rate   = Sum_k w_k*Sum_conv(G*AnisoRatio) / n
//              dipole = Sum w_k*G*AnisoRatio*cos(psi)
//                       / Sum w_k*G*AnisoRatio   (0 when that sum is 0)
//
void Scatterer::IntegrateRows(raytype intype, const R3::XYZ & dir,
                              const R3::XYZ & pol, const S2::S2Set & dirs,
                              const std::vector<GRow> & rows,
                              const std::vector<Real> * solidAngles,
                              Real & rate, Real & dipole) const {

  out_types_e c1, c2;
  if (intype == RAY_P) { c1 = GPP; c2 = GPS; }
  else                 { c1 = GSP; c2 = GSS; }
  const out_types_e convs[2] = {c1, c2};

  const Real k_in = (intype == RAY_P) ? mParams.GetKP() : mParams.GetKS();

  const S2::ThetaPhi dirTP(dir.Theta(), dir.Phi());
  const Real rot = RotationForS1(dirTP, pol);
  const R3::OrthoAxes AA(dirTP.Theta(), dirTP.Phi(), rot);  // Same
                              // construction as Phonon::Transform().

  const std::size_t n = dirs.size();

  if (solidAngles == nullptr) {

    Real sumRate = 0;
    Real sumDip  = 0;
    for (std::size_t k = 0; k < n; k++) {

      const Real costh = rows[k].CosPsi;
      const R3::XYZ e_out = AA.Express(dirs[k].XYZ());

      for (int ci = 0; ci < 2; ci++) {
        const out_types_e conv = convs[ci];
        const Real k_out = (kOutTypes[conv] == RAY_P) ?
                            mParams.GetKP() : mParams.GetKS();
        const R3::XYZ m = e_out.ScaledBy(k_out) + dir.ScaledBy(-k_in);
        const Real w = rows[k].G[conv] * mParams.AnisoRatio(m);
        sumRate += w;
        sumDip  += w * costh;
      }

    }

    rate = sumRate / (Real) n;
    dipole = (sumRate > 0) ? (sumDip / sumRate) : 0.0;

  } else {

    if (solidAngles->size() != n) {
      throw std::logic_error(
        "Scatterer::IntegrateRows: solidAngles size does not match dirs "
        "size");
    }
    const Real scale = (Real) n / (4.0 * Geometry::Pi);

    Real sumRate = 0;
    Real sumDip  = 0;
    for (std::size_t k = 0; k < n; k++) {

      const Real wk = (*solidAngles)[k] * scale;
      const Real costh = rows[k].CosPsi;
      const R3::XYZ e_out = AA.Express(dirs[k].XYZ());

      for (int ci = 0; ci < 2; ci++) {
        const out_types_e conv = convs[ci];
        const Real k_out = (kOutTypes[conv] == RAY_P) ?
                            mParams.GetKP() : mParams.GetKS();
        const R3::XYZ m = e_out.ScaledBy(k_out) + dir.ScaledBy(-k_in);
        const Real w = wk * rows[k].G[conv] * mParams.AnisoRatio(m);
        sumRate += w;
        sumDip  += w * costh;
      }

    }

    rate = sumRate / (Real) n;
    dipole = (sumRate > 0) ? (sumDip / sumRate) : 0.0;

  }

}


//////
// METHOD:  Scatterer :: IntegrateRateAndDipole()              [const]
//
// PURPOSE: Computes the exact total scattering rate (1/mean-free-path)
//          for an incoming ray of type 'intype' travelling in direction
//          'dir' with polarization direction 'pol' (both already
//          expressed in fabric coordinates), and the dipole moment of
//          the corresponding integrated scattering pattern:
//
//              dipole = Sum( cos(psi) * G ) / Sum( G )
//
//          where, for each direction dirs[k] and each of the two
//          conversions available to 'intype', G is the isotropic
//          (a=a_min) GSATO() value times AnisoRatio() at the
//          corresponding exchange wavevector, and psi is the deflection
//          angle between dirs[k] and the forward (incoming-ray)
//          direction, i.e. dirs[k].z() in the take-off-angle frame
//          (dirs[] is conventionally expressed with its own +z as the
//          forward direction -- the same convention ComputeDipoles()
//          relies on for the isotropic case). When 'solidAngles' is
//          non-null, both sums are area-weighted instead -- see
//          IntegrateRows() -- and 'solidAngles' must be the same size
//          as 'dirs' or IntegrateRows() throws std::logic_error.
//
//          Builds one GRow per direction of 'dirs' (BuildReportRows(),
//          using GetBuildThreads() workers) and hands them to the same
//          IntegrateRows() ComputeReportStats() uses, so the
//          integration loop itself is never duplicated. Used by tests
//          with their own 'dirs' sets; IntegratedRate() is a thin
//          wrapper that discards 'dipole'.
//
void Scatterer::IntegrateRateAndDipole(raytype intype, const R3::XYZ & dir,
                                       const R3::XYZ & pol,
                                       const S2::S2Set & dirs,
                                       Real & rate, Real & dipole,
                                       const std::vector<Real> *
                                         solidAngles) const {

  std::vector<GRow> rows;
  BuildReportRows(dirs, GetBuildThreads(), rows);
  IntegrateRows(intype, dir, pol, dirs, rows, solidAngles, rate, dipole);

}


//////
// METHOD:  Scatterer :: ReportDirs()                         [static]
//
// PURPOSE: Returns the fixed direction set ComputeReportStats() uses
//          to integrate the directional mean free paths and dipoles
//          shown by PrintStats() for an anisotropic Scatterer: a
//          degree-6 icosahedral tessellation (81,920 directions),
//          independent of a run's --toa-degree. Built lazily the first
//          time an anisotropic Scatterer is constructed, during the
//          single-threaded model build; see docs/superpowers/specs/
//          2026-10-16-anisotropic-scale-lengths-design.md, "Output and
//          reporting".
//
//          Always built with solid angles (TesselSphere's third
//          argument), in both quadrature modes -- about 0.66 MB and a
//          few ms, once per process -- so ComputeReportStats() can
//          integrate it area-weighted whenever Get_TOA_Quadrature() ==
//          TOA_AREA, without this lazily-built set needing to know the
//          mode in advance.
//
const S2::TesselSphere & Scatterer::ReportDirs() {
  static S2::TesselSphere dirs(S2::TESS_ICO, 6, true);
  return dirs;
}


//////
// METHOD:  Scatterer :: ComputeReportStats()         [Helper Function]
//
// PURPOSE: Populates mReportMeanFreeP, mReportDipoles, and
//          mReportMFPs, the reporting-only values PrintStats() shows
//          for an anisotropic Scatterer. Called by the constructor
//          whenever IsAnisotropic() is true.
//
//          Every anisotropic Scatterer evaluates the same 9 axis
//          states over ReportDirs() -- P travelling along East, North,
//          and Up; S travelling along each axis, polarized along each
//          of the other two axes. The legacy-style "MFP (P)" is
//          3/(sum of the 3 P rates); "MFP (S)" is 6/(sum of the 6 S
//          rates); the legacy "DM (P)"/"DM (S)" are the plain mean of
//          the corresponding dipoles.
//
//          An axisymmetric spec (IsAxisymmetric()) additionally
//          evaluates 8 incidence-angle states referenced to the az
//          axis (fabric Up): P at 0, 45, and 90 degrees from vertical;
//          S at 0 degrees (SV and SH coincide there by the medium's
//          axial symmetry, so only one value is computed), and SV, SH
//          separately at 45 and 90 degrees. GetReportMFPs() returns
//          these 8 values (instead of the 9 axis values, which still
//          feed the legacy MFP/DM values above) for an axisymmetric
//          spec.
//
//          GSATO() is evaluated once per ReportDirs() direction
//          (BuildReportRows(), GetBuildThreads() workers), not once per
//          state as before; the 9 axis states (and, for an axisymmetric
//          spec, the 8 angle states) are then integrated by
//          IntegrateRows() over that shared row table via a single
//          ParallelFor spanning every state (axis states first, then
//          angle states -- Task 4b merged what were two separate
//          ParallelFor calls into one) -- every state writes only its
//          own rate[i]/dip[i] (or r8[j]/d8[j]) slot, so there is no data
//          race. Each state's own summation stays serial and in today's
//          order, and the combination below runs only after the
//          ParallelFor has joined, serially, in the same fixed index
//          order as before: the results are bit-identical to the
//          fully-serial computation, for any thread count (see docs/
//          superpowers/specs/2026-10-17-axisymmetric-sampler-
//          performance-design.md, decision 2).
//
//          Every IntegrateRows() call below is handed ReportDirs()'s
//          SolidAngles() when Get_TOA_Quadrature() == TOA_AREA, or
//          nullptr (today's equal-weighted integration) otherwise --
//          ReportDirs() always carries solid angles (see its own header
//          comment) regardless of mode, so this is a simple mode check,
//          not a rebuild.
//
void Scatterer::ComputeReportStats() {

  const S2::TesselSphere & dirs = ReportDirs();
  const unsigned threads = GetBuildThreads();

  std::vector<GRow> rows;
  BuildReportRows(dirs, threads, rows);

  const std::vector<Real> * solidAngles =
    (Get_TOA_Quadrature() == TOA_AREA) ? &dirs.SolidAngles() : nullptr;

  const R3::XYZ E(1,0,0), N(0,1,0), U(0,0,1);

  struct AxisState {raytype in; R3::XYZ dir; R3::XYZ pol;};
  const AxisState axis[9] = {
    {RAY_P, E, E}, {RAY_P, N, N}, {RAY_P, U, U},
    {RAY_S, E, N}, {RAY_S, E, U},
    {RAY_S, N, E}, {RAY_S, N, U},
    {RAY_S, U, E}, {RAY_S, U, N},
  };

  const bool isAxi = mParams.IsAxisymmetric();

  // Angle states (axisymmetric specs only): declared/filled here, ahead
  // of the single ParallelFor below, so that ParallelFor can cover both
  // the axis and angle states in one call. Left zero-initialized (never
  // read) when !isAxi, since 'numStates' then keeps the ParallelFor from
  // ever reaching an angle-state index.
  struct AngleState {raytype in; R3::XYZ dir; R3::XYZ pol;};
  AngleState angle[8] = {};

  if (isAxi) {

    const S2::ThetaPhi tp0 (0.0,                  0.0);
    const S2::ThetaPhi tp45(45.0*Geometry::DtoR,  0.0);
    const S2::ThetaPhi tp90(90.0*Geometry::DtoR,  0.0);

    const R3::XYZ d0  = tp0.XYZ();
    const R3::XYZ d45 = tp45.XYZ();
    const R3::XYZ d90 = tp90.XYZ();

    const AngleState init[8] = {
      {RAY_P, d0,  d0},
      {RAY_P, d45, d45},
      {RAY_P, d90, d90},
      {RAY_S, d0,  tp0.ThetaHat().XYZ()},
      {RAY_S, d45, tp45.ThetaHat().XYZ()},
      {RAY_S, d45, tp45.PhiHat().XYZ()},
      {RAY_S, d90, tp90.ThetaHat().XYZ()},
      {RAY_S, d90, tp90.PhiHat().XYZ()},
    };
    for (int i = 0; i < 8; i++) {angle[i] = init[i];}

  }

  Real rate[9], dip[9];
  Real r8[8] = {}, d8[8] = {};

  const Index numStates = isAxi ? (Index) 17 : (Index) 9;
  ParallelFor(0, numStates, threads, [&](Index i) {
    if (i < 9) {
      IntegrateRows(axis[i].in, axis[i].dir, axis[i].pol, dirs, rows,
                   solidAngles, rate[i], dip[i]);
    } else {
      const Index j = i - 9;
      IntegrateRows(angle[j].in, angle[j].dir, angle[j].pol, dirs, rows,
                   solidAngles, r8[j], d8[j]);
    }
  });

  const Real sumRateP = rate[0] + rate[1] + rate[2];
  const Real sumRateS = rate[3] + rate[4] + rate[5]
                       + rate[6] + rate[7] + rate[8];
  mReportMeanFreeP[RAY_P] = 3.0 / sumRateP;
  mReportMeanFreeP[RAY_S] = 6.0 / sumRateS;
  mReportDipoles[RAY_P] = (dip[0] + dip[1] + dip[2]) / 3.0;
  mReportDipoles[RAY_S] = (dip[3] + dip[4] + dip[5]
                          + dip[6] + dip[7] + dip[8]) / 6.0;

  mReportMFPs.clear();

  if (isAxi) {

    mReportMFPs.push_back(1.0/r8[0]);                         // P0
    mReportMFPs.push_back(1.0/r8[1]);                         // P45
    mReportMFPs.push_back(1.0/r8[2]);                         // P90
    mReportMFPs.push_back(1.0/r8[3]);                         // S0
    mReportMFPs.push_back(1.0/r8[4]);                         // S45SV
    mReportMFPs.push_back(1.0/r8[5]);                         // S45SH
    mReportMFPs.push_back(1.0/r8[6]);                         // S90SV
    mReportMFPs.push_back(1.0/r8[7]);                         // S90SH

  } else {

    for (int i = 0; i < 9; i++) {
      mReportMFPs.push_back(1.0 / rate[i]);
    }

  }

}


//////
// METHOD:  Scatterer :: test_random_rayset()
//
// PURPOSE: Print on stdout a randomly-generated set of rays,
// formatted as take-off angles and symbol codes suitable for GMT
// plotting.  This is for testing/validating of scatter patterns.
//
void Scatterer::test_random_rayset(raytype intype, int count) {
  S2::S2Set & toa = (*pTOA); 
  vector< vector<int> > bins;
  bins.clear();
  bins.resize(3);               // For raytypes RAY_P, RAY_SH, RAY_SV
  for (int k=0; k<3; k++) {
    bins[k].clear();
    bins[k].resize(nTOA,0);
  }
  //
  // FILL BINS:
  //  
  for (int ctr=0; ctr<count; ctr++) {
    Phonon phon = GetRandomScatteredRelativePhonon(intype);
    raytype rt = phon.GetFullRaytype();
    int toa_index = toa.GetBestIndexFromPoint(phon.GetDirection());
    bins[rt][toa_index]++;
  }
  //
  // OUTPUT:
  //
  const char * symbol[3];
  symbol[RAY_P]  = "c";    // GMT circle code
  symbol[RAY_SH] = "-";    // GMT horiz. bar code
  symbol[RAY_SV] = "y";    // GMT vert. bar code
  for (int irt=0; irt<3; irt++) {
    for (int itoa=0; itoa<nTOA; itoa++) {
      int num = bins[irt][itoa];
      if (num != 0) {
        Real val = (Real) num / count;
        val = sqrt(val * 3.0 * nTOA) * 0.2;
        cout << setw(16) << toa[itoa].Lon()
             << setw(16) << toa[itoa].Lat()
             << setw(16) << val << " "
             << symbol[irt] << endl;
      }
    }
  }

}

//////
// METHOD:  Scatterer :: PrintAllScatteringStats()            [static]
//
//   Calls the PrintStats() method on each Scatterer object in the
//   class-maintained linked list of Scatterer objects.
//
void Scatterer::PrintAllScatteringStats() {     // STATIC METHOD

  Scatterer * next = cm_ll_first;   // First object in linked list

  bool anyAniso = false;            // Determines whether the header
  for (Scatterer * s = cm_ll_first; s != 0; s = s->mpllNext) {
    if (s->IsAnisotropic()) {anyAniso = true; break;}
  }                                 // gets the appended columns below.

  cout << "#  BEGIN SCATTERER DUMP:\n";
  cout << "#  Overrides: "
       << ((cm_MFPOverride_b) ?
           ((cm_NoDeflect_b) ? "MFPs, Deflections" : "MFPs") :
           ((cm_NoDeflect_b) ? "Deflections" : "None"))
       << std::endl;
  if (Get_TOA_Quadrature() == TOA_AREA) {
    cout << "#  TOA quadrature: area" << std::endl;
  }
  cout << "#  "
       << "    nu    eps        a    kappa         el       gam0    "
       << "     address      MFP (P)    MFP (S)      DM (P)    DM (S)";
  if (anyAniso) {
    cout << " ay az mode directional-MFPs";
  }
  cout << "\n";
  cout << "#  "
       << "====== ====== ======== ======== ========== ==========    "
       << "============    =========  =========    ========  ========\n";

  while (next != 0) {
    next->PrintStats();
    next = next->mpllNext;
  }

  cout << "#  "
       << "====== ====== ======== ======== ========== ==========    "
       << "============    =========  =========    ========  ========\n";
  cout << "#  END SCATTERERS\n";

}


//////
// METHOD:  Scatterer :: PrintStats()
//
//   Prints to stdout a breakdown of various stats of the scattering
//   object, including the parameters it was constructed on, and
//   various derived quantities, such as mean-free-paths.
//
void Scatterer::PrintStats() {

  const int width1 = 8;
  const int width2 = 12;

  cout << "   "
       << setprecision(6)
       << setw(width1-2) << mParams.GetNu() << " "
       << setw(width1-2) << mParams.GetEps() << " "
       << setw(width1)   << mParams.GetAx() << " "
       << setw(width1)   << mParams.GetKappa() << " "
       << setw(width1+2) << mParams.GetL() << " "
       << setw(width1+2) << mParams.GetGam0() << "    "
       << setw(width2)   << this << "    ";

  if (mAnisotropic) {          // The integrated legacy values from
                                // Task 6, not the majorant MFP held in
                                // mMeanFreeP[] (see GetMajorant-
                                // MeanFreePath()).
    cout << setw(width1+1) << mReportMeanFreeP[RAY_P] << "  "
         << setw(width1+1) << mReportMeanFreeP[RAY_S] << "    "
         << setprecision(4)
         << setw(width1)   << mReportDipoles[RAY_P] << "  "
         << setw(width1)   << mReportDipoles[RAY_S];
  } else {
    cout << setw(width1+1) << mMeanFreeP[RAY_P] << "  "
         << setw(width1+1) << mMeanFreeP[RAY_S] << "    "
         << setprecision(4)
         << setw(width1)   << mDipoles[RAY_P] << "  "
         << setw(width1)   << mDipoles[RAY_S];
  }

  if (mAnisotropic) {          // Append ay, az, mode, and the
                                // directional mean free paths. Mode is
                                // a property of the spec (Ruling R2),
                                // not of which sampler is in use.
    cout << "  "
         << setprecision(6)
         << setw(width1) << mParams.GetAy() << "  "
         << setw(width1) << mParams.GetAz() << "  "
         << setw(4) << (mParams.IsAxisymmetric() ? "AXI" : "ANI");
    for (std::size_t i = 0; i < mReportMFPs.size(); i++) {
      cout << "  " << setw(width1+1) << mReportMFPs[i];
    }
  }

  cout << "\n";

}
