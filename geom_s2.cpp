// geom_s2.cpp
//
#include "geom_s2.hpp"
#include "geom_r3.hpp"
#include <iostream>
#include <cmath>
#include <cerrno>     /* Error checking in math functions */
using namespace std;
using namespace S2;

// IN THIS FILE:
// Class Implementations for:
//
//   o  Classes in the S2:: namespace
//
// Search on "&&&&" to jump to the next class implementation block.
//


//////////////////////////////////////////////////////////////////////////
//************************************************************************
//* %%%%                                                             #%%#
//* #%%#    IMPLEMENTATIONS for classes in the S2:: namespace        #%%#
//* #%%#                                                             #%%#
//*         Including:
//*
//*         o  class TesselSphere
//*         o  class ThetaPhi_Set
//*         o  class Triangle
//*         o  class Node
//*
//************************************************************************


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  TesselSphere                                        ****
// ****                                                              ****
//

TesselSphere::TesselSphere(TESS_Algorithm algo, int degree,
                           bool withSolidAngles) :
  _HasSolidAngles(withSolidAngles) {
  switch (algo) {
    case TESS_TETRA :
      cerr << "Not Implemented; algo = " << algo << endl;
      break;
    case TESS_OCTA :
      cerr << "Not Implemented; algo = " << algo << endl;
      break;
    case TESS_ICO :
      init_icosahedron(degree);
      break;
    default :
      cerr << "Not Implemented; algo = " << algo << endl;
    }
}

TesselSphere::~TesselSphere(){
}

void TesselSphere::init_icosahedron(int degree) {
  // Twelve nodes make up the vertices of an icosahedron. The nodes
  // can be aligned so that each node extends in only two of the three
  // cartesian basis directions, which forms the basis of the labeling
  // convention we use here. For details see LabNotebook entry at
  // https://rainbow.phys.uconn.edu ...
  // ... /geophysics/wiki/index.php/Radiative3D_//_2013-03_MARCH

  const Real phi = (1.+sqrt(5.0))/2.0;  // Golden ratio when other
                                        // sidelength is 1.0

  Node NF = Node(   1,    0,  phi );
  Node NB = Node(  -1,    0,  phi );
  Node SF = Node(   1,    0, -phi );
  Node SB = Node(  -1,    0, -phi );
  Node FL = Node( phi,   -1,    0 );
  Node FR = Node( phi,    1,    0 );
  Node BL = Node(-phi,   -1,    0 );
  Node BR = Node(-phi,    1,    0 );
  Node RN = Node(   0,  phi,    1 );
  Node RS = Node(   0,  phi,   -1 );
  Node LN = Node(   0, -phi,    1 );
  Node LS = Node(   0, -phi,   -1 );

  // Twenty faces make up the icosahedron.  Here we create a 2-dimensional
  // array to represent the faces as node * tripples.

  //There are altogether 12 faces that touch the "edges" of the three
  //golden rectangles, and then 8 more that span the corners, but not
  //the edges.  We list the edge ones first, in order of N,S,F,B,L,R,
  //then the corner ones, starting with the four in the northern
  //hemisphere, and then the four in the southern.

  const int num_faces = 20;
  Node * fa [num_faces] [3] = {   // Face Array:
    {&NF, &NB, &LN},    // The two that touch the North edge
    {&NF, &NB, &RN},
    {&SF, &SB, &LS},    // The two that touch the South edge
    {&SF, &SB, &RS},
    {&FL, &FR, &NF},    // The two that touch the Front edge
    {&FL, &FR, &SF},
    {&BR, &BL, &NB},    // The two that touch the Back edge
    {&BR, &BL, &SB},
    {&LN, &LS, &FL},    // The two that touch the Left edge
    {&LN, &LS, &BL},
    {&RN, &RS, &FR},    // The two that touch the Right edge
    {&RN, &RS, &BR},
    {&NF, &LN, &FL},    // The four that touch corners in North hemisphere
    {&NF, &RN, &FR},
    {&NB, &LN, &BL},
    {&NB, &RN, &BR},
    {&SF, &LS, &FL},    // The four that touch corners in South hemisphere
    {&SF, &RS, &FR},
    {&SB, &LS, &BL},
    {&SB, &RS, &BR}

  };

  // Now create a temporary fully-subdivided triangle for each face,
  // and harvest the centerpoints arrays:
  Triangle * tri = 0;
  int ctr;

  if (_HasSolidAngles) {
    // _SolidAngles grows face-by-face below via grab_centerpoints_
    // from_tri()'s insert() (num_faces calls, each appending one
    // face's pow(4,degree) leaf values); without a reservation, the
    // vector's capacity-doubling growth strategy leaves old and new
    // buffers briefly live together at the final reallocation, right
    // when the vector is largest -- a transient peak-RSS cost on top
    // of the vector's own steady-state size. Reserving the exact
    // final size up front (same total this loop will insert) avoids
    // that reallocation entirely. No-op when solid angles aren't
    // requested, so equal mode is unaffected.
    _SolidAngles.reserve( num_faces * pow(4,degree) );
  }

  for (ctr=0; ctr<num_faces; ctr++) {
    tri = new Triangle(*fa[ctr][0], *fa[ctr][1], *fa[ctr][2], degree,
                       _HasSolidAngles);
    grab_centerpoints_from_tri(*tri);
    delete tri;
  }
  tri = 0;

}

void TesselSphere::grab_centerpoints_from_tri(const Triangle & tri) {
  ThetaPhi_Set::const_iterator src_begin =  tri.Centers.begin();
  ThetaPhi_Set::const_iterator src_end   =  tri.Centers.end();
  insert(end(), src_begin, src_end);

  // Carry the face's per-leaf solid angles along in lock-step with its
  // centerpoints, so _SolidAngles stays index-aligned with *this. When
  // _HasSolidAngles is false, tri.SolidAngles() is empty and this is a
  // no-op, which is what keeps SolidAngles() empty in that case.
  _SolidAngles.insert(_SolidAngles.end(),
                      tri.SolidAngles().begin(), tri.SolidAngles().end());
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  ThetaPhi                                            ****
// ****                                                              ****
//

//////
// METHOD:   ThetaPhi :: XYZ()
//
//   Returns an R3::XYZ unit-vector pointing in the direction of the
//   ThetaPhi object.
//
R3::XYZ S2::ThetaPhi::XYZ() const {
  return R3::XYZ(sin(_theta)*cos(_phi),
                 sin(_theta)*sin(_phi),
                 cos(_theta));
}


//////
// METHOD:  ThetaPhi :: ThetaHat()                   [S2:: Namespace]
//
//   Returns a unit-vector in the theta^hat (direction of increasing
//   co-latitude) direction relative to the direction of the ThetaPhi
//   object.
//
S2::ThetaPhi S2::ThetaPhi::ThetaHat() const {

  Real rth;       // Return-values: theta and phi
  Real rph;       //

  if (_theta < Geometry::Pi90) {    // Simple case

    rth = Geometry::Pi90 + _theta;
    rph = _phi;

  }                                 // Else need to flip phi in order to
  else {                            // keep our values in the favored
                                    // domain
    rth = Geometry::Pi270 - _theta;
    rph = (_phi < Geometry::Pi180) ? _phi + Geometry::Pi180
      : _phi - Geometry::Pi180;

  }

  return ThetaPhi(rth, rph);

}                                        



//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  ThetaPhi_Set                                        ****
// ****                                                              ****
//

int ThetaPhi_Set::GetBestIndexFromPoint(ThetaPhi point) {
  // Search all points in the set and return the best match to the
  // argument point.  Not exactly a fast process, but works when you
  // need it.  Calculation based on cartesian distances.
  int best_index = 0;
  Real least_dist = 2.0; // (Unit-sphere antipode distance, ie, max
                           // possible distance.)
  R3::XYZ ptXYZ = point.XYZ();
  for (Index k=0; k<size(); k++) {
    Real dist = ptXYZ.DistFrom(  (*this)[k].XYZ()  );
    if (dist < least_dist) {
      least_dist = dist;
      best_index = k;
    }
  }
  return best_index;
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  Triangle                                            ****
// ****                                                              ****
//

Triangle::Triangle(Node n1, Node n2, Node n3, int degree, bool solidAngles)
  : n1(n1),n2(n2),n3(n3),degree(degree),solidAngles(solidAngles),
    Centers(_CenterPoints) {
  //
  reserve_lists();
  //
  if (degree > 0) {
    subdivide();
  } else {
    populate_lists();
  }
}

Triangle::~Triangle(){
}

void Triangle::reserve_lists() {
  //Reserve space for the data we want this Triangle to return:
  _CenterPoints.reserve( pow(4,degree) );
  if (solidAngles) {
    // Same reallocation-avoidance reasoning as TesselSphere::
    // init_icosahedron()'s _SolidAngles reservation: subdivide()
    // below grows _SolidAngles by repeated insert() calls (one per
    // recursive subtriangle), and without a reservation that would
    // otherwise double the backing buffer's capacity (old and new
    // buffers briefly live together) on the way to its final size.
    // No-op when this triangle wasn't asked to compute solid angles,
    // so equal mode is unaffected.
    _SolidAngles.reserve( pow(4,degree) );
  }
}

void Triangle::populate_lists() {
  //Called when degree==0 to fill the data lists with the data from
  //THIS triangle.  Higher degree triangle will populate their lists
  //by concatenating the lower degree lists.
  Node center;
  center = n1 + n2 + n3; //Normalization happens on assignment
  _CenterPoints.push_back(center);

  if (solidAngles) {
    // Solid angle (steradians) subtended at the sphere's center by
    // this leaf's own three (unit-vector) vertices n1, n2, n3, via Van
    // Oosterom & Strackee's atan2 formula for the solid angle of a
    // triangle on the unit sphere given its three vertex vectors:
    //
    //   Omega = 2*atan2( |n1.(n2 x n3)|, 1 + n1.n2 + n2.n3 + n3.n1 )
    //
    // n1.(n2 x n3) is the scalar triple product, ie, (twice) the
    // signed volume of the parallelepiped spanned by the three unit
    // vectors; its absolute value makes Omega come out positive
    // regardless of the n1,n2,n3 vertex winding order. This must be
    // evaluated on THIS triangle's own three vertices -- not on the
    // quadrisection nodes q1,q2,q3, which populate_lists() never
    // touches and which don't even exist yet when degree==0 (subdivide()
    // computes them only when degree>0). A throwaway prototype
    // confirmed this reproduces the expected whole-sphere and
    // per-triangle statistics; see task-1-brief.md.
    const R3::XYZ v1(n1.x(), n1.y(), n1.z());
    const R3::XYZ v2(n2.x(), n2.y(), n2.z());
    const R3::XYZ v3(n3.x(), n3.y(), n3.z());
    const Real tripleProduct = v1.Dot(v2.Cross(v3));
    const Real cosSum        = 1.0 + v1.Dot(v2) + v2.Dot(v3) + v3.Dot(v1);
    _SolidAngles.push_back(2.0 * atan2(fabs(tripleProduct), cosSum));
  }
}

void Triangle::subdivide() {
  Triangle * tri;

  //Compute quadrisect nodes:      n2                     //////
  //                              /  \                        //
  //                            q1 -- q2                      //
  //                           /  \  /  \                     //
  //                         n1 -- q3 -- n3                   //
  q1 = n1 + n2;
  q2 = n2 + n3;
  q3 = n3 + n1;

  //Get CenterLists from each subtriangle (and, in lock-step, their
  //SolidAngles lists -- 'solidAngles' is passed down unchanged, so
  //every descendant triangle computes (or skips) the same way this one
  //was asked to, and each insert below keeps _SolidAngles index-aligned
  //with _CenterPoints exactly as the latter is built):
  //FIRST: (lower left subtriangle)
  tri = new Triangle(n1,q1,q3,degree-1,solidAngles);
  _CenterPoints.insert(_CenterPoints.end(),
			tri->Centers.begin(), tri->Centers.end());
  _SolidAngles.insert(_SolidAngles.end(),
                      tri->SolidAngles().begin(), tri->SolidAngles().end());
  delete tri;
  //SECOND: (top triangle)
  tri = new Triangle(q1,n2,q2,degree-1,solidAngles);
  _CenterPoints.insert(_CenterPoints.end(),
			tri->Centers.begin(), tri->Centers.end());
  _SolidAngles.insert(_SolidAngles.end(),
                      tri->SolidAngles().begin(), tri->SolidAngles().end());
  delete tri;
  //THIRD: (bottom right triangle)
  tri = new Triangle(q3,q2,n3,degree-1,solidAngles);
  _CenterPoints.insert(_CenterPoints.end(),
			tri->Centers.begin(), tri->Centers.end());
  _SolidAngles.insert(_SolidAngles.end(),
                      tri->SolidAngles().begin(), tri->SolidAngles().end());
  delete tri;
  //FOURTH: (middle triangle)
  tri = new Triangle(q2,q3,q1,degree-1,solidAngles);
  _CenterPoints.insert(_CenterPoints.end(),
			tri->Centers.begin(), tri->Centers.end());
  _SolidAngles.insert(_SolidAngles.end(),
                      tri->SolidAngles().begin(), tri->SolidAngles().end());
  delete tri;
  tri = 0;

  //NOTE: Some online references on STL vector<> template classes:
  //http://www.codeguru.com/cpp/cpp/cpp_mfc/stl/article.php/c4027
  //                  /C-Tutorial-A-Beginners-Guide-to-stdvector-Part-1.htm
  //http://www.cplusplus.com/reference/vector/vector/
  //http://stackoverflow.com/questions/201718/concatenating-two-stl-vectors

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  Node                                                ****
// ****                                                              ****
//

//Binary Addition: operator +
//What it does:
//  Simple member-wise addition of elements and scale counter. Returns
//  constant reference to a USNode temp variable, which survives
//  function-end as a nameless temporary in the outside scope. Result
//  of addition is NOT normalized, but this is by design, since the
//  nameless temporary might be added to yet another USNode, and
//  normalizing would ruin the equal weighting of what is now a
//  three-way addition.  Normalization will occur on assignment.  Thus
//  only nameless temporaries will ever be in a non-normalized
//  state. (At least if I've though this through right.)
const Node Node::operator + (const Node & other) const {
  //cout << "  :n:  Adding two USNodes...\n";
  Node temp;
  temp.x_ = x_ + other.x_;
  temp.y_ = y_ + other.y_;
  temp.z_ = z_ + other.z_;
  temp.iscale = iscale + other.iscale;
  return (temp);
}

//Assignment: operator =
//What it does:
//  Simple member-wise assignment, followed by normalization.  This
//  enforces the notion that USNode objects are unit-vectors, or
//  points on the surface of the unit sphere.
Node & Node::operator = (const Node & rhs) {
  // No pointers; so no need to check for &rhs == this
  x_ = rhs.x_;
  y_ = rhs.y_;
  z_ = rhs.z_;
  iscale = rhs.iscale;
  normalize();
  return(*this);
}

//Member Function: normalize()
//Purpose:
//  Normalize the vector so that it is a unit vector. 
void Node::normalize() {
  //cout << "  :n:  Taking iscale " << iscale << " down to 1.\n";
  if ((x_==0)&&(y_==0)&&(z_==0)) {
    iscale=0;   // This is a zero-vector; NOT on the unit sphere.
    return;
  } //ELSE Continue:
  Real norm = sqrt(x_*x_ + y_*y_ + z_*z_);
  x_ /= norm;
  y_ /= norm;
  z_ /= norm;
  iscale = 1;
}

// FUNCTION: printvalue()
// Prints the value of the node on stdout.
void Node::printvalue() const {
  cout << "  :n:  Node value: \n\tx = " << x_
       << "\ty = " << y_
       << "\tz = " << z_ << endl
       << "\tiscale = " << iscale << endl;
}
