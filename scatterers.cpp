// scatterers.cpp
//
#include <cmath>        /* log(), sqrt() */
#include <iomanip>      /* setw() */
#include <cstdlib>
#include <cstddef>      /* std::size_t */
#include <sstream>      /* std::ostringstream */
#include <stdexcept>    /* std::logic_error */
#include "random.hpp"
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
    mpAxi.reset(new AxiProposal(par, *pTOA, mPDists));
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
  //
  // RADIATION PATTERNS:
  //   (Populate PF arrays)
  Real gpp, gps, gsp, gss, spolv;
  for (int k = 0; k < nTOA; k++) {
    par.GSATO(toa[k], gpp, gps, gsp, gss, spolv); // Get g-values and
    mPDists[GPP].SetRelativeProb(k, gpp);         // spol value
    mPDists[GPS].SetRelativeProb(k, gps);
    mPDists[GSP].SetRelativeProb(k, gsp);
    mPDists[GSS].SetRelativeProb(k, gss);
    m_spol[k]  = spolv;
  }

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
  //   (The imfp's currently contain the sums of the G-functions.  The
  //   "probability of scattering per unit length" is the "surface
  //   average" of the G-functions. Note that the individual G-values
  //   have NOT been scaled by the size of their area element (this
  //   remains on the TODO: list - for now we make the approximation
  //   that each ToA represents an equal area), thus instead of
  //   dividing by 4Pi to get the surface average, we actually need to
  //   divide by nTOA.)  //TODO: TesselSphere needs to code area
  //   elements. See docs/toa-quadrature-bias.md for the bias this
  //   introduces and a future-fix plan.
  //
  imfp_p /= nTOA;  // TODO: This will change when above TODO fixed
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

  const R3::OrthoAxes AA(dir.Theta(), dir.Phi(), pol);  // Same
                              // construction as Phonon::Transform().
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
  result.InPol = pol;
  result.Conv = (int) conv;
  result.Toa = toa_index;
  result.Draws = draws;
  result.IsReal = accepted;

  if (accepted) {
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
  const Real pol_used = (intype == RAY_P) ?
                        CanonicalPPolarization(dir, up) : pol;

  const R3::OrthoAxes AA(dir.Theta(), dir.Phi(), pol_used);  // Same
                              // construction as Phonon::Transform().
  const R3::XYZ dirXYZ = dir.XYZ();
  const R3::XYZ u_rel(up.Dot(AA.S1()), up.Dot(AA.S2()), up.Dot(AA.S3()));

  const Real k_in = (intype == RAY_P) ? mParams.GetKP() : mParams.GetKS();

  AnisoAttempt result;
  result.InPol = pol_used;
  result.Draws = 0;

  if (!cm_MFPOverride_b) {
    // 1. Thinning test, costing no table draw (step 4 of the brief):
    const Real majorant = mpAxi->MajorantRate(intype);
    const Real binrate  = mpAxi->BinRate(intype, u_rel);
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
    const AxiProposal::Draw d = mpAxi->DrawProposal(intype, u_rel);
    conv = (out_types_e) d.Conv;
    toa_index = d.Toa;

    // 3. Compute the exchange wavevector m = k_out*e_out - k_in*dir:
    const Real k_out = (kOutTypes[conv] == RAY_P) ?
                        mParams.GetKP() : mParams.GetKS();
    const R3::XYZ e_out = AA.Express((*pTOA)[toa_index].XYZ());
    const R3::XYZ m = e_out.ScaledBy(k_out) + dirXYZ.ScaledBy(-k_in);

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

    const Real ratio = mParams.AnisoRatio(ToFabric(loc, m));
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
Real Scatterer::IntegratedRate(raytype intype, const R3::XYZ & dir,
                               const R3::XYZ & pol,
                               const S2::S2Set & dirs) const {

  Real rate, dipole;
  IntegrateRateAndDipole(intype, dir, pol, dirs, rate, dipole);
  return rate;

}


//////
// METHOD:  Scatterer :: IntegrateRateAndDipole()              [const]
//
// PURPOSE: Computes the same exact total scattering rate as
//          IntegratedRate() (see that method's header comment), and
//          additionally the dipole moment of the corresponding
//          integrated scattering pattern:
//
//              dipole = Sum( cos(psi) * G ) / Sum( G )
//
//          where, for each direction dirs[k] and each of the two
//          conversions available to 'intype', G is the isotropic
//          (a=a_min) GSATO() value times AnisoRatio() at the
//          corresponding exchange wavevector (the same weight
//          IntegratedRate() sums), and psi is the deflection angle
//          between dirs[k] and the forward (incoming-ray) direction,
//          i.e. dirs[k].z() in the take-off-angle frame (dirs[] is
//          conventionally expressed with its own +z as the forward
//          direction -- the same convention ComputeDipoles() relies on
//          for the isotropic case).
//
//          IntegratedRate() is a thin wrapper around this method that
//          discards 'dipole'.
//
void Scatterer::IntegrateRateAndDipole(raytype intype, const R3::XYZ & dir,
                                       const R3::XYZ & pol,
                                       const S2::S2Set & dirs,
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

  Real sumRate = 0;
  Real sumDip  = 0;
  const std::size_t n = dirs.size();
  for (std::size_t k = 0; k < n; k++) {

    Real gpp, gps, gsp, gss, spol;
    mParams.GSATO(dirs[k], gpp, gps, gsp, gss, spol);
    const Real gvals[NUM_OUTTYPES] = {gpp, gps, gsp, gss};
    const Real costh = dirs[k].z();   // (Unit vec, so z is cos(psi))

    const R3::XYZ e_out = AA.Express(dirs[k].XYZ());

    for (int ci = 0; ci < 2; ci++) {
      const out_types_e conv = convs[ci];
      const Real k_out = (kOutTypes[conv] == RAY_P) ?
                          mParams.GetKP() : mParams.GetKS();
      const R3::XYZ m = e_out.ScaledBy(k_out) + dir.ScaledBy(-k_in);
      const Real w = gvals[conv] * mParams.AnisoRatio(m);
      sumRate += w;
      sumDip  += w * costh;
    }

  }

  rate = sumRate / (Real) n;
  dipole = (sumRate > 0) ? (sumDip / sumRate) : 0.0;

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
const S2::TesselSphere & Scatterer::ReportDirs() {
  static S2::TesselSphere dirs(S2::TESS_ICO, 6);
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
void Scatterer::ComputeReportStats() {

  const S2::TesselSphere & dirs = ReportDirs();

  const R3::XYZ E(1,0,0), N(0,1,0), U(0,0,1);

  struct AxisState {raytype in; R3::XYZ dir; R3::XYZ pol;};
  const AxisState axis[9] = {
    {RAY_P, E, E}, {RAY_P, N, N}, {RAY_P, U, U},
    {RAY_S, E, N}, {RAY_S, E, U},
    {RAY_S, N, E}, {RAY_S, N, U},
    {RAY_S, U, E}, {RAY_S, U, N},
  };

  Real rate[9], dip[9];
  for (int i = 0; i < 9; i++) {
    IntegrateRateAndDipole(axis[i].in, axis[i].dir, axis[i].pol, dirs,
                           rate[i], dip[i]);
  }

  const Real sumRateP = rate[0] + rate[1] + rate[2];
  const Real sumRateS = rate[3] + rate[4] + rate[5]
                       + rate[6] + rate[7] + rate[8];
  mReportMeanFreeP[RAY_P] = 3.0 / sumRateP;
  mReportMeanFreeP[RAY_S] = 6.0 / sumRateS;
  mReportDipoles[RAY_P] = (dip[0] + dip[1] + dip[2]) / 3.0;
  mReportDipoles[RAY_S] = (dip[3] + dip[4] + dip[5]
                          + dip[6] + dip[7] + dip[8]) / 6.0;

  mReportMFPs.clear();

  if (mParams.IsAxisymmetric()) {

    const S2::ThetaPhi tp0 (0.0,                  0.0);
    const S2::ThetaPhi tp45(45.0*Geometry::DtoR,  0.0);
    const S2::ThetaPhi tp90(90.0*Geometry::DtoR,  0.0);

    const R3::XYZ d0  = tp0.XYZ();
    const R3::XYZ d45 = tp45.XYZ();
    const R3::XYZ d90 = tp90.XYZ();

    Real r, d;

    IntegrateRateAndDipole(RAY_P, d0,  d0,  dirs, r, d);
    mReportMFPs.push_back(1.0/r);                             // P0
    IntegrateRateAndDipole(RAY_P, d45, d45, dirs, r, d);
    mReportMFPs.push_back(1.0/r);                             // P45
    IntegrateRateAndDipole(RAY_P, d90, d90, dirs, r, d);
    mReportMFPs.push_back(1.0/r);                             // P90

    IntegrateRateAndDipole(RAY_S, d0, tp0.ThetaHat().XYZ(), dirs, r, d);
    mReportMFPs.push_back(1.0/r);                             // S0

    IntegrateRateAndDipole(RAY_S, d45, tp45.ThetaHat().XYZ(), dirs, r, d);
    mReportMFPs.push_back(1.0/r);                             // S45SV
    IntegrateRateAndDipole(RAY_S, d45, tp45.PhiHat().XYZ(), dirs, r, d);
    mReportMFPs.push_back(1.0/r);                             // S45SH

    IntegrateRateAndDipole(RAY_S, d90, tp90.ThetaHat().XYZ(), dirs, r, d);
    mReportMFPs.push_back(1.0/r);                             // S90SV
    IntegrateRateAndDipole(RAY_S, d90, tp90.PhiHat().XYZ(), dirs, r, d);
    mReportMFPs.push_back(1.0/r);                             // S90SH

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
