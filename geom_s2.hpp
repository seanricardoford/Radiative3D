// geom_s2.hpp
//
// The geom library defines namespaces containing classes and other
// constructs related to the geometry of spaces in various dimensions.
// This file defines:
//
// o Namespace: S2:
//
//    For geometry on the surface of a sphere.  More info is in the
//    comments just inside the namespace statement.
//
//
#ifndef GEOM_S2_H_
#define GEOM_S2_H_
//
#include "geom_base.hpp"
#include <vector>
#include <cmath>


//////
// *** FORWARD DECLARATIONS:
//

namespace R3{       /* Definitions in geom_r3.hpp */
  class XYZ;
};


//////
// *** CLASSES:
//

//__________________________________________________________________________
//**************************************************************************
// NAMESPACE: S2
// PURPOSE:
//
//   To provide a set of classes, types, and constants for
//   representing geometry on the surface of a sphere (so-called S2
//   space).
//
//   Also provides classes for discretizing the space on the surface
//   of the sphere, such as the TesselSphere class, which uses a
//   tessellation algorithm to provide a set of points more-or-less
//   uniformely distributed on the surface. (Turns out it's not
//   actually uniform, but might be close enough.)
//
//__________________________________________________________________________
//**************************************************************************

namespace S2 {

//////
// CONSTANTS:
//

  enum TESS_Algorithm {TESS_TETRA,  // Tetrahedron Starting Pattern 
                       TESS_OCTA,   // Octahedron
                       TESS_ICO};   // Icosahedron
                       // Used by the TesselSphere class to determine
                       // the initial arrangement of triangles before
                       // subdivision begins.

//////
// CLASSES:
//


//CLASS:  S2::Node
//ENCAPS: Unit-Sphere Node: a point on the surface of the unit sphere
//Capabilities:
//  * Self-normalizes on the assignment operator.
//  * Does not normalize on arrithmatic operators. (allows for easy
//    averaging of points, i.e., a=b+c+d puts the three-way average of
//    b,c,d into a)
//Examples:
//  x = a + b + c
//    Sets x equal to the NORMALIZED vector average of a & b & c.
//  a + b + c
//    This represents the UN-normalized vector sum of a + b + c
//  somefunction(a + b)
//    Here, the function will receive an UN-normalized node object.
//    (But might be difference in rec by ref or rec by val.)
//Normalized vs. Unnormalized State:
//  In general, named instances of Node will always be in a
//  normalized state, since their value can only be set by
//  initialization or assignment.  Unnamed temporaries, though, such
//  as the results of addition operators, CAN be in an unnormalized
//  state.  There is but one exception to this that I know of: It is
//  possible to initialize a Node from an unnamed temporary like so:
//      Node node(Node(1,0,0) + Node(0,1,0));
//  In this case, the expected call to the copy-constructor is
//  suppressed and the new object simply binds to the unnamed
//  temporary result of the addition, and node will end up in an
//  unnormalized state.  (For this reason, I didn't bother overriding
//  the default copy constructor, as it doesn't really buy me
//  anything.  But anyway, aside from that oddball kind of thing,
//  named objects will always be normalized.
//Zero State:
//  It is also possible for the node to be in a "zero" state, with
//  x=y=z=0, and iscale=1.  This happens if either 1) you initialize a
//  Node without arguments, or 2) you take a two-way average of
//  antipodes (or any set of n points whose vector sum is zero).  At
//  present, this state does not raise an error condition.  Although
//  perhaps it should.
//
class Node {
protected:

  // :::::::::::::::::::::::::::::::::
  // ::: Member Data  (Node Class) :::
  // :::::::::::::::::::::::::::::::::

  Real x_;
  Real y_;
  Real z_;
  int iscale; //Normalized when iscale==1


public:

  // ::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Node Class) :::
  // ::::::::::::::::::::::::::::::::::

  Node() :
    x_(0),y_(0),z_(0),iscale(0) {
  }
  Node(Real xin, Real yin, Real zin) :
    x_(xin),y_(yin),z_(zin),iscale(-1) {
    normalize();
  }
  ~Node(){
  };


  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (Node Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  Real x() const {return x_;}
  Real y() const {return y_;}
  Real z() const {return z_;}

  void printvalue() const;


  // ::::::::::::::::::::::::::::::::::::::::
  // ::: Operator Overloads  (Node Class) :::
  // ::::::::::::::::::::::::::::::::::::::::

  const Node operator + (const Node &) const;
  Node & operator = (const Node &);

private:

  // :::::::::::::::::::::::::::::::::::::
  // ::: Private Methods  (Node Class) :::
  // :::::::::::::::::::::::::::::::::::::

  void normalize();

};


//////
//CLASS:  S2::ThetaPhi
//ENCAPS: A single Theta,Phi tuple.
//Capabilities:
//  * Initializable expressly from Theta, Phi parameters
//  * Initializable implicitly from a Tessel::Node parameter
//  * Initializable implicitly from a Tessel::LatLon parameter
//Details:
//  * When initializing from parameters or an object that holds x,y,z
//    values, uses a RHS coordinate system where +z is North pole, +x
//    goes through the equator/prime-meridian, and +y points at 90-deg
//    East longitude.
//  * When initializing from a (Lat,Lon) object, assume phi==0
//    represents prime-meridian, and phi increases with eastward
//    longitude.
//
class ThetaPhi {
protected:

  // :::::::::::::::::::::::::::::::::::::
  // ::: Member Data  (ThetaPhi Class) :::
  // :::::::::::::::::::::::::::::::::::::

  Real _theta;
  Real _phi;

public:

  // ::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (ThetaPhi Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  ThetaPhi(Real theta, Real phi) :
    _theta(theta), _phi(phi) {
  }
  ThetaPhi(const Node & n) {
    _theta = acos(n.z()); // Assumes n is normalized
    _phi = atan2(n.y(),n.x());
  }
  ThetaPhi() :
    _theta(0), _phi(0) {
  }
  ~ThetaPhi() {
  }


  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (ThetaPhi Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  Real Lat() const {return (90.0 - Geometry::RtoD * _theta);}
  Real Lon() const {return (Geometry::RtoD * _phi);}
  Real Theta() const {return _theta;}
  Real Phi() const {return _phi;}
  Real x() const {return sin(_theta)*cos(_phi);}
  Real y() const {return sin(_theta)*sin(_phi);}
  Real z() const {return cos(_theta);}

  R3::XYZ XYZ() const;
            // Returns the ThetaPhi object converted to an XYZ
            // object. (Unit vector pointing in same direction as
            // ThetaPhi object.)

  ThetaPhi ThetaHat() const;
            // Returns a unit-vector in the theta^hat (direction of
            // increasing co-latitude) direction relative to the
            // direction of the ThetaPhi object.

  ThetaPhi PhiHat() const {
    //                Like ThetaHat(), but returns phi^hat
    //                unit vector.
    return ThetaPhi(Geometry::Pi90, 
                    (_phi < Geometry::Pi270) 
                        ? _phi + Geometry::Pi90
                        : _phi - Geometry::Pi270);
  }


  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (ThetaPhi Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  void SetTheta(Real theta) {_theta = theta;}
  void SetPhi(Real phi) {_phi = phi;}
  void Set(Real theta, Real phi) {
    _theta = theta;
    _phi = phi;
  }

};


//////
// CLASS:  S2::ThetaPhi_Set -FROM- std::vector<ThetaPhi>
//
// PURPOSE: To encapsulate a VECTOR of ThetaPhi objects, with some
//          added features.
//
class ThetaPhi_Set : public std::vector<ThetaPhi> {
public:
  int GetBestIndexFromPoint(ThetaPhi point);
};


//////
//CLASS:  S2::Triangle
//
//  Encapsulates Great-Circle Triangles on the unit sphere with the
//  following capabilities:
//
//  * Capable of subdividing itself into smaller triangles.
//  * Can return a list of the centerpoints of all subtriangles
//    represented as unit vectors.
//  * Functions recursively
//
class Triangle {
protected:

  // :::::::::::::::::::::::::::::::::::::
  // ::: Member Data  (Triangle Class) :::
  // :::::::::::::::::::::::::::::::::::::

  Node n1,n2,n3;  // Three nodes make up a triangle
  Node q1,q2,q3;  // And three nodes quadrisect a triangle

  int degree;     //Degree of this triangle (When degree==0, we do
                  //not subdivide further)

  bool solidAngles;  // True if this triangle (and, by propagation,
                     // every triangle produced by subdivide()) should
                     // compute and keep a per-leaf solid angle
                     // alongside each leaf's centerpoint. Set once at
                     // construction (default false; see constructor)
                     // and never changed afterwards. Off by default
                     // so existing callers (and their output) are
                     // completely unaffected.

  ThetaPhi_Set _CenterPoints;  // List of center-points within the
                               // fully subdivided triangle

  std::vector<Real> _SolidAngles;  // Parallel to _CenterPoints: the
                               // solid angle (steradians) of the leaf
                               // triangle that produced the
                               // correspondingly-indexed centerpoint.
                               // Left empty unless 'solidAngles' is
                               // true.

public:

  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Public Member Data  (Triangle Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  const ThetaPhi_Set & Centers; // RO access to _CenterPoints
                // NOTE:/WARNING: *** Assigning read-only references
                // to class members comes with certain dangers: this
                // may invalidate the ability to write an assignment
                // operator for this class (which at present I don't
                // think I need, but...)  or even to use the default
                // (is there one?)  assignment operator.  Should look
                // into this.  Access function may be better idea.

public:

  // ::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Triangle Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  Triangle(Node, Node, Node, int degree, bool solidAngles = false);

  ~Triangle();


  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (Triangle Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  const std::vector<Real> & SolidAngles() const {return _SolidAngles;}
                // One entry per leaf triangle, in the same order as
                // Centers (index-aligned with it). Empty unless this
                // Triangle was built with solidAngles=true.

private:

  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Private Methods  (Triangle Class) :::
  // :::::::::::::::::::::::::::::::::::::::::

  void subdivide();
  void reserve_lists();
  void populate_lists();

};


//////
// CLASS:  S2::TesselSphere
//
// ENCAPSULATES the properties of a tessellated sphere, ie, a sphere
// represented by a finite number of discrete surface elements.
// Provides a way (eventually several ways) of getting a set of points
// suitable for sampling functions defined on the surface of the unit
// sphere, such as the moment-tensor radiation patterns for earthquake
// events, or Sato and Fehler's scattering functions.  The icosahedron
// based quadrisection algorithm produces points which are
// approximately, though not exactly, uniformly spread across the
// surface of the sphere.  (Surface-area-per-point variability is
// approximately xxx% (TBD))
//
class TesselSphere : public ThetaPhi_Set {
public:
  TesselSphere(TESS_Algorithm algo, int degree, bool withSolidAngles = false);
  ~TesselSphere();

  bool HasSolidAngles() const {return _HasSolidAngles;}
                // True iff constructed with withSolidAngles=true, in
                // which case SolidAngles() below is populated.

  const std::vector<Real> & SolidAngles() const {return _SolidAngles;}
                // Omega_k (steradians), index-aligned with *this (ie,
                // with the inherited ThetaPhi_Set of take-off
                // directions): SolidAngles()[k] is the solid angle of
                // the leaf triangle whose centerpoint is (*this)[k].
                // Empty unless HasSolidAngles() is true.

private:
  void init_icosahedron(int degree);
  void grab_centerpoints_from_tri(const Triangle & tri);

  bool _HasSolidAngles;             // Set once at construction time;
                                    // see withSolidAngles above.
  std::vector<Real> _SolidAngles;   // Backing store for SolidAngles();
                                    // left empty unless
                                    // _HasSolidAngles is true.

};

};// end namespace S2
///
//

#endif //#ifndef GEOM_S2_H_
