// parallel.hpp
//
// ParallelFor(): a minimal, deterministic-coverage parallel-for over
// std::thread, built for Radiative3D's construction-time table
// building -- AxiProposal's P/S bin tables (aniso_axi.cpp) and
// Scatterer::ComputeReportStats()'s direction rows and report states
// (scatterers.cpp). See docs/superpowers/specs/2026-10-17-axisymmetric-
// sampler-performance-design.md, decision 2.
//
// Every caller in this codebase writes each index's result to a
// disjoint location (a table entry, an array slot, owned by that index
// alone) and combines the results afterward serially, in a fixed index
// order -- that combination, not ParallelFor itself, is what makes the
// final result independent of thread count (see aniso_axi.cpp's
// AxiProposal constructor and scatterers.cpp's ComputeReportStats()/
// IntegrateRateAndDipole()). ParallelFor's own contract is narrower: it
// guarantees only that body(i) runs exactly once for every i in
// [begin,end) (when nothing throws), never which thread runs a given
// index or in what order indices complete.
//
// IMPLEMENTATION NOTE -- body(i) has exactly one call site on purpose:
//
//   An earlier version of this file special-cased 'threads <= 1' with
//   its own plain "for (i=begin;i<end;i++) body(i);" loop, separate
//   from the claim-and-call loop ('worker', below) used for the
//   multi-thread case. That turned out to be an actual bit-identity
//   bug, not just a style choice: with that second call site, a
//   caller's body -- e.g. aniso_axi.cpp's BuildBinTable()/BoundRatio()
//   -- compiled to measurably different code at the two sites (GCC
//   -O3's default fused-multiply-add contraction applied at one call
//   site but not the other), so the *exact same* floating-point
//   expression rounded differently depending on whether 'threads' was
//   1 or more -- a real, reproduced discrepancy (confirmed by printing
//   intermediate values with -ffp-contract=off, which made both sides
//   agree, and by checking which of the two values matched a frozen
//   reference build: the single-call-site version below is the one
//   that does). Giving body(i) a single call site -- 'worker' is
//   called directly, inline, for the no-thread case, instead of being
//   bypassed -- removes the discrepancy entirely: see
//   test_axi_tables_thread_invariant and test_report_stats_thread_
//   invariant (tests/anisotropy_tests.cpp), which exist specifically
//   to catch a regression here. Do not reintroduce a separate inline
//   loop as a "simplification".
//
#ifndef PARALLEL_H_
#define PARALLEL_H_
//
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>
#include "typedefs.hpp"    /* Index */

//////
// TEMPLATE FUNCTION:  ParallelFor()
//
// PURPOSE: Calls body(i) exactly once for every i in the half-open
//          range [begin,end), using up to 'threads' std::thread
//          workers. 'threads' <= 1, or a range of at most one index,
//          runs the claim-and-call loop directly on the calling
//          thread instead of spawning anything -- no std::thread is
//          created in either case.
//
//          Indices are claimed from a single shared
//          std::atomic<Index> counter, so which thread ends up running
//          a given index, and the order in which indices complete, is
//          scheduling-dependent -- callers must not depend on either.
//          What callers *can* rely on: every index in [begin,end) runs
//          exactly once (absent a thrown exception), and body() is
//          never called with an index outside that range, regardless
//          of 'threads' -- and, in particular, that body(i)'s own
//          generated code (and so its floating-point results) does not
//          depend on 'threads' either (see the IMPLEMENTATION NOTE
//          above).
//
//          If body() throws for one or more indices: every spawned
//          thread is still joined (never detached, never
//          std::terminate()); indices not yet claimed at the time a
//          failure is recorded are left unclaimed ("skipped"); and,
//          once every thread has joined, the exception object from the
//          numerically lowest index that threw is rethrown to the
//          caller -- deterministic regardless of which thread hit it
//          first or the order in which indices happened to finish. Any
//          exception from a higher-numbered index is discarded.
//
template <class Body>
void ParallelFor(Index begin, Index end, unsigned threads, Body body) {

  if (end <= begin) {return;}
  const Index n = end - begin;

  const bool runInline = (threads <= 1 || n <= 1);
  const unsigned nThreads = runInline ? 1 :
    ((threads < (unsigned) n) ? threads : (unsigned) n);

  std::atomic<Index> cursor(begin);
  std::atomic<Index> failAt(end);     // 'end' doubles as the "nothing
                                      // has failed yet" sentinel: it is
                                      // never itself a valid claimed
                                      // index.
  std::exception_ptr failException;
  std::mutex excMtx;

  // The sole call site for body(i) -- see the IMPLEMENTATION NOTE
  // above for why 'runInline' below calls this function directly
  // rather than using a separate loop.
  auto worker = [&]() {
    while (true) {
      const Index i = cursor.fetch_add(1);
      if (i >= end) {return;}
      if (i >= failAt.load(std::memory_order_relaxed)) {continue;}
                      // A lower index has already failed: this index
                      // stays claimed (the counter keeps advancing
                      // towards 'end'), but running its body is
                      // skipped -- the outcome is already decided.
      try {
        body(i);
      } catch (...) {
        std::lock_guard<std::mutex> lock(excMtx);
        if (i < failAt.load(std::memory_order_relaxed)) {
          failAt.store(i, std::memory_order_relaxed);
          failException = std::current_exception();
        }
      }
    }
  };

  if (runInline) {
    worker();
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nThreads);
    for (unsigned t = 0; t < nThreads; t++) {
      pool.emplace_back(worker);
    }
    for (std::vector<std::thread>::iterator it = pool.begin();
        it != pool.end(); ++it) {
      it->join();
    }
  }

  if (failException) {std::rethrow_exception(failException);}

}//

///
#endif //#ifndef PARALLEL_H_
//
