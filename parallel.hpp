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
// IMPLEMENTATION NOTE -- indices are claimed in contiguous chunks, from
// a 64-bit cursor:
//
//   An earlier version of 'worker' claimed a single index per
//   std::atomic::fetch_add(), so a large-'n' call -- AxiProposal's
//   per-representative loops at 81,920 representatives, or a
//   degree-9 take-off-angle set's 5,242,880 directions -- touched the
//   one shared counter once per index. Profiling (docs/superpowers/
//   specs/2026-10-17-axisymmetric-sampler-performance-design.md's
//   follow-up, task-5b-brief.md) found that contention, not the
//   callers' own arithmetic, dominating wall time at those sizes.
//   Each worker below instead claims a contiguous run of 'grain'
//   indices per fetch_add ('grain = max(1, n/(threads*8))', sized so
//   each worker gets about 8 chunks -- enough that a worker finishing
//   its share early still finds smaller remaining work to steal),
//   cutting the number of atomic operations from 'n' down to roughly
//   'threads*8'. The claimed-exactly-once and lowest-failing-index-
//   rethrow guarantees are unchanged: a chunk boundary is just where
//   one claim's for-loop ends and the next begins, 'failAt' is still
//   checked before every single body(i) call (not just once per
//   chunk), and every chunk -- including one claimed after a failure
//   was recorded -- still counts as "claimed" whether or not its
//   indices' bodies actually run, so the set of indices handed to
//   *some* worker is still exactly [begin,end), partitioned without
//   gaps or overlap, regardless of 'threads' or 'grain'.
//
//   The cursor and its failure sentinel use a fixed 64-bit integer
//   type ('Cursor', below) regardless of Index's own width (Index is
//   'unsigned int' -- typedefs.hpp): claiming 'grain' indices at a
//   time, rather than one at a time, means the last few fetch_add()
//   calls before the range is exhausted can legitimately push the
//   cursor past 'end' by up to 'grain' per racing thread. A cursor no
//   wider than Index could, in principle, wrap around through zero
//   after enough such overshoot and start claiming (invalid,
//   already-finished) low indices again; a 64-bit one has no
//   realistic path to wrapping for any 'end' that itself fits in an
//   Index.
//
#ifndef PARALLEL_H_
#define PARALLEL_H_
//
#include <algorithm>
#include <atomic>
#include <cstdint>
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
//          Indices are claimed in contiguous chunks from a single
//          shared 64-bit atomic cursor (see the IMPLEMENTATION NOTE
//          above), so which thread ends up running a given index, and
//          the order in which indices -- or chunks -- complete, is
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

  // 'Cursor': a fixed 64-bit type regardless of Index's own width -- see
  // the IMPLEMENTATION NOTE above on why chunked claiming needs this.
  typedef std::uint64_t Cursor;
  const Cursor cBegin = (Cursor) begin;
  const Cursor cEnd   = (Cursor) end;
  const Cursor grain  = std::max<Cursor>(
    1, (Cursor) n / ((Cursor) nThreads * 8));

  std::atomic<Cursor> cursor(cBegin);
  std::atomic<Cursor> failAt(cEnd);   // 'cEnd' doubles as the "nothing
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
      const Cursor lo = cursor.fetch_add(grain);
      if (lo >= cEnd) {return;}
      const Cursor hi = (cEnd - lo > grain) ? lo + grain : cEnd;
      for (Cursor i = lo; i < hi; i++) {
        if (i >= failAt.load(std::memory_order_relaxed)) {break;}
                      // A lower index has already failed: the rest of
                      // this chunk stays claimed (no other worker will
                      // ever see these indices), but running their
                      // bodies is skipped -- the outcome is already
                      // decided.
        try {
          body((Index) i);
        } catch (...) {
          std::lock_guard<std::mutex> lock(excMtx);
          if (i < failAt.load(std::memory_order_relaxed)) {
            failAt.store(i, std::memory_order_relaxed);
            failException = std::current_exception();
          }
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
