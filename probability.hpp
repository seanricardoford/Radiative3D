// probability.hpp
//
// Develops a class designed to contain probability "distributions,"
// or indexed sets over which each index has some given probability of
// being selected.  Provides methods for selecting random indices
// using a random number generator.
//
#ifndef PROBABILITY_H_
#define PROBABILITY_H_
//
#include <cstdint>
#include <random>
#include <vector>
#include "typedefs.hpp"  // CHECK!!!

//////
// *** FORWARD DECLARATIONS:
//

//////
// *** CLASSES:
//

//////
// CLASS: RandomEngine
//
//   Deterministic random stream used by stochastic simulation code.
//
class RandomEngine {
public:
  explicit RandomEngine(std::uint64_t seed) : mEngine(seed) {}

  std::uint64_t NextUInt64() { return mEngine(); }

  Real Uniform01() {
    return std::uniform_real_distribution<Real>(0.0, 1.0)(mEngine);
  }

  static std::uint64_t SeedForStream(std::uint64_t seed,
                                     std::uint64_t stream) {
    std::uint64_t value = seed + 0x9e3779b97f4a7c15ULL * (stream + 1);
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
  }

  static RandomEngine &Default() {
    static thread_local RandomEngine engine(
        SeedForStream(0x72616433645f6469ULL, 0));
    return engine;
  }

private:
  std::mt19937_64 mEngine;
};

//////
// CLASS:   ProbDist
//
//   Encapsulates a probability distribution, or an array of
//   non-negative numbers used to establish the relative probability
//   of a set of discrete, indexed elements.  Provides also a set of
//   methods useful for randomly selecting an index from the set.
//
//   STORAGE: Numerical probabilities are stored one of two different
//   ways, with a boolean member encoding the storage method.  The
//   class will of its own accord switch between the storage methods
//   depending on which representation suits a desired action.  The
//   two representations are (1) differential and (2) accumulated, aka
//   integrated.  If the representation is differential, then the
//   elements of the array contain the relative probility of occurence
//   of the element.  If the representation is integrated, then the
//   elements of the array contain the cumulative sum of relative
//   probability as the array index increases.  The differential
//   representation is typically used for population the distribution,
//   and the integrated representation is used for randomly selecting
//   an index (since the random chooser requires this form).  Usually,
//   the ProbDist object begins life in differential form and then
//   switches to integral form when the first method that requires
//   that form is called.
//
//   NORMALIZATION: The numerical values stored in the arrays are NOT
//   normalized, meaning (equivalently) that the sum of the
//   differential values, and the final integrated value, are not
//   necessarily unitary. (Ie, might not sum to 1).  This is intended
//   design, and the implication is that when code-users set
//   individual probability elements, they are setting RELATIVE
//   probability weights, and not necessarily absolute (normalized)
//   probabilities. This means the distribution will have some
//   "magnitude," and that that magnitude will be accessible to the
//   user as it may presumably be of some use to the user.  The
//   chooser methods, which select an element index at random, will
//   scale appropriately and function as expected regardless of this
//   magnitude.
//
//   TERMINOLOGY: "[Relative] Probability Weight" will refer to the
//   un-normalized probability of a particular single element.
//   "Differential Probability" will refer to the normalized
//   probability of a particular single element.  The Probability
//   Weights (PW) and the Differential Probabilities (DP) stand in the
//   relation PW = Magnitude * DP.  "Cumulative Probability" will
//   refer to the sum of the Probability Weight of a particular
//   element and all elements with lower indices than the particular
//   element.  Cumulative Probabilities are NOT normalized.  There is
//   at present no adopted terminology for normalized cumulative
//   probability as there is at present no use case for it.
//
class ProbDist {
protected:

  // :::::::::::::::::::::::::::::::::::::
  // ::: Member Data  (ProbDist Class) :::
  // :::::::::::::::::::::::::::::::::::::

  bool mbIntegrated;    // false: Representation is un-normalized
                        //        differential (I.e, as Relative
                        //        Probability Weights).
                        // true:  Representation is un-normalized
                        //        integrated (I.e., as Cumulative
                        //        Probabilities).

  bool mbAllocated;     // False until array is resized to mSize.
  int  mSize;           // Target size of the array. 
                        //
                        // We use late allocation, by which I mean the
                        // size of the array is a required constructor
                        // parameter, but the array (based on
                        // std::vector) isn't resized to make room
                        // until the first element is set, post-
                        // construction.  This simplifies allocation
                        // of vectors of ProbDist objects because
                        // potentially very LARGE objects don't need
                        // to be passed to the vector::resize(n,val)
                        // method.
                        //

  std::vector<Real> mDist;  // Contains the distribution, or array
                            // of probability values.


public:

  // ::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (ProbDist Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  ProbDist(int size) :
    mbIntegrated (false),
    mbAllocated  (false),
    mSize        (size)
  {}


  // :::::::::::::::::::::::::::::::::::::
  // ::: Set Methods  (ProbDist Class) :::
  // :::::::::::::::::::::::::::::::::::::

  void SetRelativeProb(Index idx, Real prob) {
    if (!mbAllocated) {         // Set a relative probability weight,
      mDist.resize(mSize);      // with appropriate checks.
      mbAllocated=true;}        //
    if (!mbIntegrated) {
      if (prob<0.0)   //
        {AckSpit();}  // Complain if negative (invalid probability)
      mDist[idx] = prob;
    } else {
      AckSpit();      // Complain if we've already integrated
    }                 //
  }


  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Get-Info Methods  (ProbDist Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  Real GetMagnitude();      // Return the sum-total (non-normalized)
                            // probability of the distribution.

  Real GetRelativeProb(Index);  // Return a (non-normalized) relative
                                // probability weight.  This is the
                                // same value as what the user would
                                // have set for this index.

  Real GetDiffProb(Index);  // Return a (normalized) differential
                            // probability for a given index (equal to
                            // the relative prob divided by
                            // magnitude).


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Generate Result Methods  (ProbDist Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  Index GetRandomIndex(RandomEngine &); // Use the supplied random stream.
  Index GetRandomIndex();               // Diagnostic convenience overload.

protected:

  // :::::::::::::::::::::::::::::::::::
  // ::: Mutations  (ProbDist Class) :::
  // :::::::::::::::::::::::::::::::::::

  void Integrate();     // Switch internal representation from
                        // differential to cumulative.


  // ::::::::::::::::::::::::::::::::::::::::
  // ::: Error Handling  (ProbDist Class) :::
  // ::::::::::::::::::::::::::::::::::::::::

  void AckSpit(); // Register a complaint if class used improperly

};

///
#endif //#ifndef PROBABILITY_H_
//
