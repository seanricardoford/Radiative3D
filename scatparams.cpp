// scatparams.cpp
//
#include <iostream>
#include <cmath>      /* for pow(), sin(), cos(), tgamma() */
#include "scatparams.hpp"

using namespace std;

// in this file:
// CLASS IMPLEMENTATIONS FOR:
//
//   o  Class ScatterParams
//


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  ScatterParams                                       ****
// ****                                                              ****
//

//
// STATIC MEMBERS:
//

Real ScatterParams::cm_omega = 1.0;         // Meaningless default - Actual
bool ScatterParams::cm_omega_known = false; // value is to be set in the
                                            // Model constructor via a call
                                            // to SetFrequencyHertz()

//////
// METHOD:  ScatterParams :: CompareRoughly()
//
//   Very rough comparrison between two ScatterParams objects,
//   returned as a sum of squared differences in parameter values.
//   This is very unsophisticated.
//
Real ScatterParams::CompareRoughly(const ScatterParams & other) const {
  Real dnu = (other.nu - nu);
  Real deps = (other.eps - eps);
  Real da = (other.a - a);
  Real dax = (other.ax - ax);
  Real day = (other.ay - ay);
  Real daz = (other.az - az);
  Real dkappa = (other.kappa - kappa);
  Real del = (other.el - el);
  Real dgam0 = (other.gam0 - gam0);
  Real sum = dnu*dnu + deps*deps
           + da*da + dax*dax + day*day + daz*daz
           + dkappa*dkappa
           + del*del + dgam0*dgam0;
  return sum;
}


//////
// METHOD:  ScatterParams :: GSATO()
//
// PURPOSE: Computes (4.52) from Sato and Fehler using von Karman
// autocorellation function (PSATO).  The implementation here mimics
// the implementation in PSPhonon.
//
// Comment block from PSPhonon: (Note: Parameters el, nu, gam0, eps,
// and a are member variables and do not need to be passed, and psi
// and zeta are packed in toa.)
//
//! GSATO computes (4.52) from Sato and Fehler (exponential autocor)
//!   Inputs:  psi  =  spherical coor. angle (radians)
//!            zeta =  sph. coor angle from x3 (radians)
//!            el   =  S-wave wavenumber (=om/beta0)
//!            nu   =  density vs. velocity pert. scaling (see 4.48)
//!            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
//!            eps  =  RMS velocity perturbation
//!            a    =  correlation distance
//!   Returns: gpp,gps,gsp,gss  =  from eqn. (4.52)
//!            spol =  S-to-S scattered S polarization (radians)
//!                 =  0 for pure psi direction
//!
void ScatterParams::GSATO(S2::S2Point toa,
                          Real & gpp, Real & gps, 
                          Real & gsp, Real & gss, 
                          Real & spol) const {
  const Real pi4  = 4. * Geometry::Pi;
  const Real el4  = pow(el,4);
  const Real gam2 = pow(gam0,2);
  const Real psi = toa.Theta();

  Real xpp, xps, xsp, xss_psi, xss_zeta;
  Real arg;

  XSATO(toa, xpp, xps, xsp, xss_psi, xss_zeta); // Compute X-values

  Real xpp2 = xpp * xpp;                      // Square them all
  Real xps2 = xps * xps;
  Real xsp2 = xsp * xsp;
  Real xss_psi2  = xss_psi * xss_psi;
  Real xss_zeta2 = xss_zeta * xss_zeta;

  // GPP:
  arg = (2.*el/gam0)*sin(psi/2.);               // For PSATO()
  gpp = (el4/pi4) * xpp2 * PSATO(arg);
  if (gpp < 1.e-30) gpp=0.; // PSPhonon sets this lower-limit
                            // value. Not sure why. But I retain
                            // it here.
  // GPS and GSP share the same PSATO() argument (see scatparams.hpp's
  // GSATO() comment / decision 8 in the design note): computed once here
  // instead of twice, since the two original calls already had to agree
  // bit-for-bit (same 'arg', same const method, no mutable state) -- this
  // just removes the redundant second call, it does not change what
  // either value is.
  arg = (el/gam0) * sqrt(1. + gam2 - 2.*gam0*cos(psi));
  const Real psato_gps_gsp = PSATO(arg);

  // GPS:
  gps = (1./gam0) * (el4/pi4) * xps2 * psato_gps_gsp;
  if (gps < 1.e-30) gps=0.;

  // GSP:
  gsp = gam0 * (el4/pi4) * xsp2 * psato_gps_gsp;
  if (gsp < 1.e-30) gsp=0.;

  // GSS:
  arg = 2.*el*sin(psi/2.);
  gss = (el4/pi4) * (xss_psi2 + xss_zeta2) * PSATO(arg);
  if (gss < 1.e-30) gss=0.;
  spol=atan2(xss_zeta,xss_psi); // Allows the whole range of
                                // polarization angles (-Pi to Pi).
  return;

}


//////
// METHOD:  ScatterParams :: XSATO()
//
// PURPOSE: Computes (4.50) from Sato and Fehler
//
// COMMENT BLOCK FROM PSPhonon: (Note: Parameters nu and gam0 are
// implemented here as member variables and do not need to be passed
// in; and note also that psi and zeta are both packed into
// "take-off-angle" toa as toa.phi() and toa.theta(), respectively.)
//
//! XSATO computes (4.50) from Sato and Fehler
//!   Inputs:  psi  =  spherical coor. angle (radians)
//!            zeta =  sph. coor angle from x3 (radians)
//!            nu   =  density vs. velocity pert. scaling (see 4.48)
//!            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
//!   Returns: xpp, xps, xsp, xss_psi, xss_zeta  =  from eqn. (4.50)
//!
void ScatterParams::XSATO(S2::S2Point toa,
                          Real & xpp, Real & xps, 
                          Real & xsp, Real & xss_psi, 
                          Real & xss_zeta) const {

  const Real gam2  = gam0 * gam0;
  const Real cpsi  = cos(toa.Theta());     // psi <--> toa.Theta()
  const Real c2psi = cos(2.*toa.Theta());
  const Real spsi  = sin(toa.Theta());
  const Real czeta = cos(toa.Phi());       // zeta <--> toa.Phi()
  const Real szeta = sin(toa.Phi());
  const Real spsi2 = spsi*spsi;

  xpp = (1./gam2) * (nu * (-1. + cpsi + (2./gam2) * spsi2) 
                     - 2. + (4./gam2) * spsi2);

  xps = -spsi*(nu*(1.-(2./gam0)*cpsi)-(4./gam0)*cpsi);

  xsp = (1./gam2) * spsi*czeta*(nu*(1.-(2./gam0)*cpsi) 
                                - (4./gam0)*cpsi);

  xss_psi = czeta*(nu*(cpsi-c2psi)-2.*c2psi);

  xss_zeta = szeta*(nu*(cpsi-1.)+2.*cpsi);

  return;

}


//////
// METHOD:  ScatterParams :: PSATO()
//
// PURPOSE: Compute Sato and Fehler's P(m) function using von Karman
// autocorrelation.  (Note that when k=0.5, von Karman is equivalent
// to exponential, which is what is used in PSPhonon.)
//
// COMMENT BLOCK from PSPhonon: (Note that eps and a are member
// variables of the class, and thus do not need to be passed in.)
//
//! function EXPSATO computes (2.10) from Sato and Fehler
//!    Inputs:  eps  =  RMS velocity perturbation
//!             a    =  correlation distance
//!             m    =  wavenumber
//!    Returns: P(m) =  PSDF (Power Spectral Density Function)
//!
Real ScatterParams::PSATO(Real m) const {

  const Real denom = pow((1.+a*a*m*m),(kappa+1.5));

  return mPSatoNumer/denom;
}


//////
// METHOD:  ScatterParams :: RecomputePSatoNumer()
//
// PURPOSE: (Re)computes mPSatoNumer and mPSatoAnisoNumer -- PSATO()'s and
// PSATO_Aniso()'s constant numerators -- from the current eps, a, ax, ay,
// az and kappa. Called by every constructor (there is no public setter
// for any of those six members, so construction is the only time this
// needs to run).
//
// Each expression below is written exactly once, here, textually
// identical to what PSATO()/PSATO_Aniso() used to compute locally on
// every call, in the same left-to-right operation order (never
// re-associated: e.g. caching tgamma(kappa+1.5)/tgamma(kappa) as one
// number first would regroup (A*B)/C into A*(B/C), which can round
// differently). Writing it once, in one function called from every
// constructor, rather than copying the expression into each constructor
// body, also avoids the class of bug documented in parallel.hpp's
// IMPLEMENTATION NOTE: the same floating-point expression, written at
// two (or more) source locations, can be compiled to differently-
// rounded code (an FMA-contraction discrepancy was observed there
// between two call sites of one lambda).
//
void ScatterParams::RecomputePSatoNumer() {
  const Real pi32  = pow(Geometry::Pi,1.5); // pi^(3/2)

  mPSatoNumer = (8.*pi32*eps*eps*a*a*a)
              * tgamma(kappa+1.5)
              / tgamma(kappa);

  mPSatoAnisoNumer = (8.*pi32*eps*eps*ax*ay*az)
                    * tgamma(kappa+1.5)
                    / tgamma(kappa);
}


//////
// METHOD:  ScatterParams :: AnisoMajorantFactor()
//
// PURPOSE: Returns M = ax*ay*az/a^3 (a == min(ax,ay,az)), the factor by
// which the anisotropic spectrum PSATO_Aniso() can exceed the isotropic
// PSATO() evaluated at the same |m| (see AnisoRatio()). This is also the
// factor by which the anisotropic-medium majorant mean free path is
// shortened relative to the a=min isotropic mean free path used as the
// delta-tracking proposal distribution. M==1 when the spec is isotropic.
//
Real ScatterParams::AnisoMajorantFactor() const {
  return (ax*ay*az)/(a*a*a);
}


//////
// METHOD:  ScatterParams :: PSATO_Aniso()
//
// PURPOSE: Computes the anisotropic von Karman power spectral density
// function, Sato (2019) Eq. 3, generalizing PSATO() to three independent
// scale lengths. 'm' is the exchange wavevector expressed in the fabric
// frame (ie, its x,y,z components are measured along the ax,ay,az axes,
// respectively -- NOT necessarily the simulation's internal XYZ axes).
//
//   P(m) = 8 pi^(3/2) eps^2 ax ay az Gamma(kappa+3/2)
//          / [ Gamma(kappa) (1 + ax^2 mx^2 + ay^2 my^2 + az^2 mz^2)^(kappa+3/2) ]
//
// When ax==ay==az==a, this reduces exactly to PSATO(m.Mag()).
//
Real ScatterParams::PSATO_Aniso(const R3::XYZ & m) const {
  const Real mx = m.x();
  const Real my = m.y();
  const Real mz = m.z();
  const Real denom = pow((1. + ax*ax*mx*mx + ay*ay*my*my + az*az*mz*mz),
                         (kappa+1.5));

  return mPSatoAnisoNumer/denom;
}


//////
// METHOD:  ScatterParams :: AnisoRatio()
//
// PURPOSE: Returns PSATO_Aniso(m)/PSATO(|m|), the ratio of the
// anisotropic spectrum to the isotropic (a=min(ax,ay,az)) spectrum
// evaluated at the same wavevector. This ratio is always positive and
// never exceeds AnisoMajorantFactor() (equality is reached whenever the
// denominator term along the shortest fabric axis dominates, e.g. m
// purely along that axis, or trivially at m=0).
//
// Computed as AnisoMajorantFactor()*AnisoAcceptance(m) rather than
// PSATO_Aniso(m)/PSATO(m.Mag()): algebraically identical (the eps^2 and
// Gamma-function prefactors in PSATO_Aniso and PSATO cancel exactly,
// leaving the ax*ay*az/a^3 majorant factor times the acceptance-ratio
// bracket), but this form is well-defined at eps=0 -- PSATO_Aniso/PSATO
// is 0/0 there -- and it skips the two separate pow() calls (one per
// denominator) that PSATO_Aniso(m)/PSATO(m.Mag()) would otherwise need.
//
Real ScatterParams::AnisoRatio(const R3::XYZ & m) const {
  return AnisoMajorantFactor() * AnisoAcceptance(m);
}


//////
// METHOD:  ScatterParams :: AnisoAcceptance()
//
// PURPOSE: Returns AnisoRatio(m)/AnisoMajorantFactor(), ie,
//
//   [ (1 + a^2 |m|^2) / (1 + ax^2 mx^2 + ay^2 my^2 + az^2 mz^2) ]^(kappa+3/2)
//
// This is the delta-tracking accept probability for a proposal (direction
// and conversion type) drawn from the isotropic a=min(ax,ay,az) tables:
// always in (0,1], reaching 1 only when the proposal direction's fabric
// components make the isotropic and anisotropic spectra agree (always
// true when the spec is isotropic).
//
Real ScatterParams::AnisoAcceptance(const R3::XYZ & m) const {
  const Real mx = m.x();
  const Real my = m.y();
  const Real mz = m.z();
  const Real numer = 1. + a*a*m.MagSquared();
  const Real denom = 1. + ax*ax*mx*mx + ay*ay*my*my + az*az*mz*mz;

  return pow(numer/denom, kappa+1.5);
}


//////
// METHOD:  ScatterParams :: AxiRatio()
//
// PURPOSE: AnisoRatio(m), for an axisymmetric spec (ax==ay), computed
// directly from m2 = |m|^2 and mu2 = (m.Dot(up))^2 -- the local-vertical
// 'up' component of 'm', squared -- with no fabric-frame transform. Since
// ax==ay, the denominator's horizontal contribution ax^2*mx^2+ay^2*my^2
// depends only on the horizontal magnitude mx^2+my^2 = m2-mu2 (clamped at
// 0 to absorb rounding), whichever way that horizontal part is split
// between the fabric x and y axes; the vertical contribution is
// az^2*mu2. This is algebraically identical to AnisoRatio(m) for any
// fabric-frame 'm' with the same |m| and |m.up|, including one actually
// expressed via Scatterer::ToFabric(). See docs/superpowers/specs/
// 2026-10-17-axisymmetric-sampler-performance-design.md, decision 4.
//
// Not meaningful when ax!=ay (the horizontal split would then matter).
//
Real ScatterParams::AxiRatio(Real m2, Real mu2) const {
  const Real mh2 = std::max(m2 - mu2, (Real)0.0);
  const Real numer = 1. + a*a*m2;
  const Real denom = 1. + ax*ax*mh2 + az*az*mu2;

  return AnisoMajorantFactor() * pow(numer/denom, kappa+1.5);
}


//////
// METHOD:   ScatterParams :: AckSpit()
//
//   Raise a complaint if class used improperly.
//
void ScatterParams::AckSpit() {

  throw(std::logic_error(
    "Can't construct ScatterParams before frequency is known."));

}
