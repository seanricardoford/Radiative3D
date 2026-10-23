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
#include <memory>           /* std::unique_ptr */
#include "sources.hpp"
#include "scatparams.hpp"
#include "aniso_axi.hpp"    /* AxiProposal: SAMPLER_AXISYMMETRIC tables */
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
                                // path; axisymmetric specs (ax==ay!=az)
                                // use SAMPLER_AXISYMMETRIC; all other
                                // specs use SAMPLER_GENERAL. An
                                // explicit, non-AUTO choice always
                                // forces the named sampling path, even
                                // for an isotropic spec; forcing
                                // SAMPLER_AXISYMMETRIC on a spec with
                                // ax!=ay throws (see the constructor).

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

  static const raytype kOutTypes[4];   // Maps conversion types (GPP,
                                       // GPS, GSP, GSS) to output
                                       // raytypes (RAY_P or RAY_S).
                                       // Shared by
                                       // GetRandomScatteredRelative-
                                       // Phonon(), AttemptAnisotropic-
                                       // Scatter(), and
                                       // AttemptAxisymmetricScatter()
                                       // (previously three identical
                                       // local copies).


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

  static unsigned cm_BuildThreads;  // Construction-time worker count:
                                    // how many std::thread workers
                                    // (ParallelFor(), parallel.hpp) a
                                    // subsequently-constructed
                                    // Scatterer's AxiProposal and
                                    // ComputeReportStats() use to build
                                    // their tables/rows/states. Affects
                                    // only how they are built, never
                                    // what they contain -- see
                                    // SetBuildThreads().

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

  static void SetBuildThreads(unsigned n) {
    cm_BuildThreads = (n == 0) ? 1 : n;
  }
                      // Construction-time worker count used by every
                      // Scatterer constructed from this call onward:
                      // up to 'n' std::thread workers (ParallelFor(),
                      // parallel.hpp) build its AxiProposal's bin
                      // tables (if any) and its ComputeReportStats()
                      // rows/states, with results bit-identical for any
                      // value -- see the AxiProposal constructor and
                      // ComputeReportStats(). 0 is treated as 1 (there
                      // is no "automatic" sentinel at this layer; a
                      // caller such as Model::Model() resolves
                      // "automatic" to a concrete worker count before
                      // calling this).

  static unsigned GetBuildThreads() {return cm_BuildThreads;}
                      // Default 1 (serial construction).


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

  bool  mUseAxi;               // True when this (necessarily
                                // mAnisotropic) Scatterer uses the
                                // axisymmetric sampler (mpAxi) rather
                                // than the general sampler. Set at
                                // construction: SAMPLER_AXISYMMETRIC
                                // forces it; SAMPLER_AUTO sets it iff
                                // par.IsAxisymmetric(); SAMPLER_GENERAL
                                // never sets it.

  std::unique_ptr<AxiProposal> mpAxi;  // Axisymmetric-sampler proposal
                                       // tables; built in the
                                       // constructor iff mUseAxi, null
                                       // otherwise.

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

  // Scatterers are always heap-allocated via GetScattererMatchingParams()
  // and managed through the intrusive mpllPrev/mpllNext linked list
  // (never copied or moved); an anisotropic, axisymmetric-sampler
  // instance's mpAxi (AxiProposal) holds a const reference directly
  // into this object's mPDists member, which a move would silently
  // invalidate (the vector's contents move with the source object's
  // mPDists, leaving mpAxi's reference pointing at the moved-from,
  // now-empty vector). Delete copy and move explicitly so that
  // invariant can't be broken by a future caller.
  Scatterer(const Scatterer &) = delete;
  Scatterer & operator=(const Scatterer &) = delete;
  Scatterer(Scatterer &&) = delete;
  Scatterer & operator=(Scatterer &&) = delete;


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

  struct GRow {
    Real G[NUM_OUTTYPES]; // gpp,gps,gsp,gss (ScatterParams::GSATO()'s
                          // values at this row's direction).
    Real CosPsi;          // That direction's z-component: cos(psi), the
                          // deflection angle in the take-off-angle
                          // frame (see IntegrateRows()).
  };                      // One per report direction; see
                          // BuildReportRows().

  void BuildReportRows(const S2::S2Set & dirs, unsigned threads,
                       std::vector<GRow> & rows) const;
                      // Fills 'rows' (resized to dirs.size()) with one
                      // GRow per direction, via ParallelFor() (parallel.
                      // hpp) over direction indices on up to 'threads'
                      // workers. Each row is a pure function of
                      // dirs[k] and mParams (GSATO() touches no mutable
                      // state), so 'rows' is bit-identical to a serial
                      // build regardless of 'threads' -- see docs/
                      // superpowers/specs/2026-10-17-axisymmetric-
                      // sampler-performance-design.md, decision 2.
                      // Shared by ComputeReportStats() (dirs ==
                      // ReportDirs(), the rows reused across every
                      // state) and IntegrateRateAndDipole() (its own
                      // 'dirs'), so GSATO() is never evaluated more
                      // than once per (Scatterer,direction) pair.

  void IntegrateRows(raytype intype, const R3::XYZ & dir,
                     const R3::XYZ & pol, const S2::S2Set & dirs,
                     const std::vector<GRow> & rows,
                     Real & rate, Real & dipole) const;
                      // The row-based integrator shared by
                      // IntegrateRateAndDipole() and
                      // ComputeReportStats(): the same rate/dipole sum
                      // IntegrateRateAndDipole()'s header comment
                      // describes, but reading each direction's GSATO()
                      // values and cos(psi) from the precomputed 'rows'
                      // (BuildReportRows()) instead of calling GSATO()
                      // itself. The summation stays serial, in
                      // k-ascending order -- identical term-for-term to
                      // today's single combined loop -- so a call's
                      // result never depends on how (or on how many
                      // threads) 'rows' was built.

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
                      // anisotropic Scatterer this is the isotropic
                      // (a=a_min) mean free path divided by the
                      // sampler's majorant factor: mParams.
                      // AnisoMajorantFactor() for the general sampler,
                      // or mpAxi->MajorantRate(rt)'s reciprocal for the
                      // axisymmetric sampler (see the constructor).

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

  static Real CanonicalPPolarization(const S2::ThetaPhi & dir,
                                     const R3::XYZ & up);
                      // Returns RotationForS1(dir, (up.d)d - up), where
                      // d=dir.XYZ(): the polarization-rotation angle
                      // that puts the local vertical 'up' into the
                      // (S1,S3) plane with a non-positive S1 component,
                      // so the axisymmetric sampler's P tables (which
                      // have no chi dimension) can assume chi==0. Used
                      // in place of a real polarization angle for
                      // incoming P phonons, which have none.

protected:

  AnisoAttempt AttemptAxisymmetricScatter(raytype intype,
                                          const R3::XYZ & loc,
                                          const S2::ThetaPhi & dir,
                                          Real pol);
                      // The SAMPLER_AXISYMMETRIC implementation of
                      // AttemptAnisotropicScatter(); see that method's
                      // header comment. Called from
                      // AttemptAnisotropicScatter() when mUseAxi.

  Phonon BuildRelativePhonon(out_types_e conv, Index toa_index,
                            raytype intype) const;
                      // Builds the "relative" Phonon for an accepted
                      // anisotropic-scatter event: looks up the output
                      // raytype via kOutTypes[conv], the polarization
                      // (m_spol[toa_index] for GSS, else 0), and -- if
                      // Scatterer::SetNoDeflect() is in effect -- a
                      // non-deflected dummy Phonon instead. Shared by
                      // AttemptAnisotropicScatter() and
                      // AttemptAxisymmetricScatter() (previously two
                      // identical finishing blocks).

public:


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
