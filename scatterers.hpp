// scatterers.hpp
//
// This file develops a set of classes that represent seismic
// scattering events as "sources" of re-radiated elastodynamic
// phonons.
//
// The classes developed here derive from the base class developed in
// sources.hpp.  Including this file in your code module will
// automatically include sources.hpp as well.
//
//
#ifndef SCATTERERS_H_
#define SCATTERERS_H_
//
#include "sources.hpp"
#include "scatparams.hpp"
#include "phonons.hpp"      /* Needed for the Phonon member of
                               Scatterer::AnisoAttempt. */

//////
// CLASSES: Forward Declarations
//

//////
// TYPEDEFS:
//

//////
// CLASSES: Definitions
//

//////
// CLASS:  Scatterer -FROM- PhononSource
//
// PURPOSE:
//
//   This class represents "Scatterers," as treated theoretically in the
//   manner described by Sato and Fehler.  The Scatterer object derives
//   from, and thus is a kind of, the PhononSource object.  Ie, the
//   Scatterer object is a "source" of scattered phonons.
//
// MEMORY FOOTPRINT:
//
//   This class, like any class derived from PhononSource, could be
//   considered somewhat "heavy".  It maintains a large array of
//   probabilities correllating to each possible takeoff-angle, and the
//   Scatterer class in particular also maintains a similar array of
//   polarization values.  Because of this, care and consideration should
//   be given to reducing the total number of Scatterer objects that need
//   to be allocated.  In a model with many volume cells, giving each one
//   its own scatterer might not be a good idea.
//
// ALLOCATION PROCEDURE:
//
//   To ameliorate this situation, this class proveds a static method
//   called GetScattererMatchingParams(), the purpose of which is to
//   maintain a list of allocated Scatterers, and to return a pointer to
//   an existing Scatterer if one that sufficiently matches the requested
//   parameters already exists. Otherwise, it constructs and records a new
//   scatterer and returns a pointer to that one.  The model-building code
//   should use this member rather than the 'new' keyword to link volume
//   cells to Scatterers.
//
// CACHEING PROCEDURE:
//
//   To further ameliorate the memory demands of the scattering class,
//   the class will maintain the ability to flush the Probability
//   arrays if the net memory demands become too large.  Each
//   scattering object will keep a use-count so that more frequently
//   used Probability arrays can be preferentially kept compared to
//   less frequently used arrays.  In essence, the most frequently
//   used arrays will be kept in a cache, whereas the less frequently
//   used arrays will be re-generated as needed.
//
class Scatterer : public PhononSource {
protected:

  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Enums and Constants  (Scatterer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  enum in_types_e {IN_P, IN_S, NUM_INTYPES};
                                // IN_P and IN_S agree with raytype
                                // enum and can be used interchangably
                                // where P/S raytypes needed.

public:

  enum sampler_e {SAMPLER_AUTO, SAMPLER_GENERAL, SAMPLER_AXISYMMETRIC};
                                // Chooses which sampling algorithm a
                                // Scatterer uses for anisotropic
                                // ScatterParams (ax,ay,az not all
                                // equal). AUTO defers to the spec:
                                // isotropic specs use today's exact
                                // path; axisymmetric specs use
                                // SAMPLER_AXISYMMETRIC (falls back to
                                // SAMPLER_GENERAL until that sampler
                                // is implemented); all other specs use
                                // SAMPLER_GENERAL. An explicit, non-
                                // AUTO choice always forces the
                                // anisotropic (general) sampling path,
                                // even for an isotropic spec.

protected:

  enum out_types_e {GPP, GPS, GSP, GSS, NUM_OUTTYPES};
                                // Four basic conversions denote the
                                // four scattering distributions
                                // (shapes) used. (Note that Sato
                                // distinguishes between G_S->SH and
                                // G_S->SV, but we keep just a net
                                // distribution for S->S and
                                // distinguish polarization via a
                                // separate m_spol array.)


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Class-Static Infrastructure  (Scatterer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::


  static Scatterer  * cm_ll_first;  // Points to first element of the
                                    // linked list of Scatterer objects.

  static Scatterer   * cm_ll_last;  // Points to last element of the
                                    // linked list of Scatterer objects.

  static Count cm_ll_count;         // Keeps tally of how many elements 
                                    // are on the linked list.

  static bool cm_MFPOverride_b;         // True to overide mean-free paths
  static Real cm_MFPOverrides[RAY_NBT]; // with explicit values. (Applies
                                        // to all Scatterer objects.)

  static bool cm_NoDeflect_b;       // True to suppress deflection (phonon
                                    // direction, raytype, polarization
                                    // will not be changed by scattering
                                    // event). (Applies to all objects.)

                // NOTE: if MFP and deflection are BOTH overridden, then
                // ScatterParams have no remaining effect, and no more than
                // one object need be created. GetScattererMatchingParams()
                // should detect this and squash par to limit list to one
                // scatterer.

public:
  ;
  static Scatterer * 
  GetScattererMatchingParams(ScatterParams par);
                      // Retrieve or create a Scatterer object
                      // matching to within tolerances the provided
                      // parameters.

  static Count GetCollectionSize() {return cm_ll_count;}
                      // How many Scatterers are in the collection of
                      // Scatterer objects known to the class.

  static void OverrideMFP(Real mfp_P, Real mfp_S) {
    cm_MFPOverride_b = true;
    cm_MFPOverrides[RAY_P] = mfp_P;
    cm_MFPOverrides[RAY_S] = mfp_S;
  }

  static void SetNoDeflect() {
    cm_NoDeflect_b = true;
  }


protected:

  // :::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Member Variables  (Scatterer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::

  Scatterer        * mpllPrev;  // Linked-list pointers
  Scatterer        * mpllNext;  //

  Real  mMeanFreeP[RAY_NBT];    // Mean Free Path before scattering.
                                // Indexed by raytype. 

  Real    mDipoles[RAY_NBT];    // Dipole moments of the distributions.
  Real  mQuadpoles[RAY_NBT];    // Quadrupole moments.
                                //
                                // These values are computed because
                                // they make an interesting way to
                                // quantitatively summarize the
                                // scattering shapes, but they are not
                                // needed or used in the actual
                                // simulation. They are just computed
                                // so that they can be reported to the
                                // user if desired.

  ScatterParams     mParams;    // Parameters used to generate the
                                // probability distributions.
                                // Retained incase the distribution
                                // needs to be recomputed after a
                                // cache flush. (The ability to dump
                                // the distribution for memory
                                // conservation is a planned future
                                // feature.)

  std::vector<Real>  m_spol;    // Polarization angles for S->S
                                // conversions

  bool  mAnisotropic;          // True unless this Scatterer uses
                                // today's exact isotropic sampling
                                // path (GetRandomScatteredRelative-
                                // Phonon()). Set at construction from
                                // 'par' and the requested sampler_e
                                // choice; see IsAnisotropic().

  Real mReportMeanFreeP[RAY_NBT]; // Reporting-only mean free paths for
                                  // an anisotropic Scatterer: legacy
                                  // "MFP (P)"/"MFP (S)" values averaged
                                  // over incidence along the three
                                  // fabric axes (see ComputeReport-
                                  // Stats()). Distinct from mMeanFreeP,
                                  // which (for an anisotropic
                                  // Scatterer) holds the majorant MFP
                                  // used to draw tentative collision
                                  // distances, not a value meant for
                                  // display. Unset (and unused) for an
                                  // isotropic Scatterer.

  Real mReportDipoles[RAY_NBT];   // Reporting-only "DM (P)"/"DM (S)"
                                  // values for an anisotropic
                                  // Scatterer, the mean of the dipoles
                                  // of the same axis states that feed
                                  // mReportMeanFreeP. Unset (and
                                  // unused) for an isotropic Scatterer.

  std::vector<Real> mReportMFPs;  // Directional mean free paths, in
                                  // dump order, for an anisotropic
                                  // Scatterer; see GetReportMFPs().
                                  // Empty for an isotropic Scatterer.


public:

  // :::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Scatterer Class) :::
  // :::::::::::::::::::::::::::::::::::::::

  Scatterer(ScatterParams par, sampler_e choice = SAMPLER_AUTO);


protected:
;
  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Helper Functions  (Scatterer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

  void PopulateProbDists(ScatterParams par);
  void PopulateWholeProbs();
  void ComputeMFPs();
  void ComputeDipoles();
  void ComputeQuads();
  void ComputeReportStats();   // Populates mReportMeanFreeP,
                               // mReportDipoles, and mReportMFPs.
                               // Called by the constructor whenever
                               // IsAnisotropic() is true.

  static const S2::TesselSphere & ReportDirs();
                               // The fixed degree-6 direction set
                               // ComputeReportStats() integrates over,
                               // independent of a run's --toa-degree.
                               // Built lazily on first use.


public:
;
  // :::::::::::::::::::::::::::::::::::::::::::::
  // ::: Manipulate Methods  (Scatterer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::

  void SetMeanFreePathsPS(Real, Real);
  //            Though the MFPs are computed automatically by the
  //            constructor, this method can be used to change them
  //            after the fact.  This is mainly for diagnostic and
  //            development/debugging purposes.


  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Compute Methods  (Scatterer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  Real GetRandomPathLength(raytype intype);
  Phonon GetRandomScatteredRelativePhonon(raytype intype);

  bool IsAnisotropic() const {return mAnisotropic;}
                      // True unless this Scatterer is using today's
                      // exact isotropic sampling path.

  Real GetMajorantMeanFreePath(raytype rt) const {return mMeanFreeP[rt];}
                      // The mean free path used to draw tentative
                      // (proposal) collision distances. For an
                      // anisotropic (general-sampler) Scatterer this
                      // is the isotropic (a=a_min) mean free path
                      // divided by mParams.AnisoMajorantFactor().

  struct AnisoAttempt {
    bool    IsReal;  // True if this attempt produced a real (accepted)
                     // scattering event.
    Real    InPol;   // The incoming polarization-rotation angle passed
                     // in (echoed back for convenience).
    Phonon  Rel;     // The "relative" phonon describing the proposed
                     // (if IsReal) deflection/conversion, in the same
                     // sense as GetRandomScatteredRelativePhonon()'s
                     // return value. Meaningless when !IsReal.
    int     Conv;    // The out_types_e conversion drawn on the final
                     // (accepted, or -- without MFP override -- only)
                     // attempt.
    Index   Toa;     // The take-off-angle index drawn on that attempt.
    int     Draws;   // Number of proposal draws made by this call (1
                     // unless Scatterer::OverrideMFP() is in effect,
                     // in which case draws repeat until accepted).

    AnisoAttempt() :
      IsReal(false), InPol(0),
      Rel(S2::ThetaPhi(0,0), RAY_P),
      Conv(-1), Toa(0), Draws(0)
    {}
  };

  AnisoAttempt AttemptAnisotropicScatter(raytype intype,
                                         const R3::XYZ & loc,
                                         const S2::ThetaPhi & dir,
                                         Real pol);
                      // Exact (thinning/delta-tracking) sampler for
                      // anisotropic ScatterParams. 'loc' is used only
                      // to look up the local (East,North,Up) fabric
                      // frame via ECS; 'dir' and 'pol' describe the
                      // incoming phonon's direction and polarization-
                      // rotation angle exactly as Phonon::Transform()
                      // would construct them. Must not be called when
                      // !IsAnisotropic().

  Real IntegratedRate(raytype intype, const R3::XYZ & dir,
                      const R3::XYZ & pol, const S2::S2Set & dirs) const;
                      // Exact total scattering rate (1/mean-free-path)
                      // for an incoming ray of type 'intype' travelling
                      // in direction 'dir' with polarization direction
                      // 'pol' (both already expressed in fabric
                      // coordinates), computed by averaging, over the
                      // relative-direction grid 'dirs', the isotropic
                      // (a=a_min) GSATO() values times AnisoRatio() at
                      // the corresponding exchange wavevector, summed
                      // over the two conversions available to
                      // 'intype'. Used to validate AttemptAnisotropic-
                      // Scatter()'s realized event rate. A thin wrapper
                      // around IntegrateRateAndDipole() that discards
                      // the dipole.

  void IntegrateRateAndDipole(raytype intype, const R3::XYZ & dir,
                              const R3::XYZ & pol, const S2::S2Set & dirs,
                              Real & rate, Real & dipole) const;
                      // Computes the same 'rate' as IntegratedRate(),
                      // plus the dipole moment of the corresponding
                      // integrated scattering pattern, Sum(cos(psi)*G)
                      // / Sum(G), where psi is the deflection angle
                      // (angle between 'dirs[k]' and the forward
                      // direction, i.e. 'dirs[k].z()' in the local
                      // take-off-angle frame) and G is each direction's
                      // G-value*AnisoRatio() weight -- the same
                      // convention ComputeDipoles() uses for the
                      // isotropic case, applied here to the anisotropic
                      // integral.

  const std::vector<Real> & GetReportMFPs() const {return mReportMFPs;}
                      // Directional mean free paths, in dump order:
                      // for an axisymmetric spec, the 8 incidence-angle
                      // values (P0,P45,P90,S0,S45SV,S45SH,S90SV,S90SH);
                      // for any other anisotropic spec, the 9 axis
                      // values (P_E,P_N,P_U,S_E/N,S_E/U,S_N/E,S_N/U,
                      // S_U/E,S_U/N). Empty unless IsAnisotropic().

  static R3::XYZ ToFabric(const R3::XYZ & loc, const R3::XYZ & v);
                      // Expresses vector 'v' in the local fabric
                      // (East,North,Up) frame at 'loc', i.e. returns
                      // (v.Dot(ECS.GetEast(loc)), v.Dot(ECS.GetNorth(loc)),
                      // v.Dot(ECS.GetUp(loc))).

  static Real RotationForS1(const S2::ThetaPhi & dir, const R3::XYZ & s1);
                      // Returns the rotation angle 'rot' such that
                      // R3::OrthoAxes(dir.Theta(),dir.Phi(),rot).S1()
                      // is parallel to the projection of 's1' normal
                      // to 'dir'; returns 0 if that projection is
                      // ~zero (e.g. 's1' parallel to 'dir').


  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Testing Methods  (Scatterer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  void test_random_rayset(raytype intype, // Output to stdout a
                          int count);     // randomly generated set
                                          // of rays; For testing
                                          // scatter patterns.

  static void PrintAllScatteringStats();  // Walks the linked-list of
                                          // all scatterers and
                                          // outputs their stats.

  void PrintStats();    // Dump stats for particular scatterer.


}; // class Scatterer
////

///
#endif //#ifndef SCATTERERS_H_
//
