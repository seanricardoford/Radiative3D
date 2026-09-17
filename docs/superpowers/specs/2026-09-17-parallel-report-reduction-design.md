# Parallel Report Reduction Design

## Status

Approved for implementation on the `parallel-performance` branch; implemented,
verified, and merged into local `master` at `8b9713d`. The feature branch was
deleted after the fast-forward merge. See `PLANS.md` and
`docs/DEVELOPMENT.md` for the current state and follow-up work.

## Goal

Make multi-worker Radiative3D simulations scale better by removing shared
reporting and seismometer accumulation from the normal propagation hot path,
while preserving the existing command-line interface, seeded per-phonon random
streams, and serial defaults.

## Problem

The current worker implementation parallelizes phonon propagation but routes
every event through the global `DataReporter`. Each callback acquires one
mutex, even when its textual report is disabled. Collection callbacks hold the
same mutex while testing all seismometers and mutating shared time bins. The
worker scheduler and phonon diagnostic identifiers also perform frequent
shared atomic operations. The result is enough synchronization and cache-line
traffic to make the four-worker 10M-phonon benchmark slower than the serial
reference.

## Chosen architecture

Each simulation worker receives a private reporting context. The context owns
all mutable simulation-side reporting state:

- lost, timeout, and invalid counters;
- invalidity diagnostic flags; and
- one private copy of every seismometer's time-bin accumulation.

The existing `DataReporter` retains immutable run configuration and seismometer
geometry. Event callbacks receive the worker context explicitly. They perform
geometry tests against the immutable configuration and write only to the
worker's private bins. After all workers join, `DataReporter` reduces the
contexts into the canonical seismometer traces and aggregate counters before
post-simulation output is generated.

Normal simulations with textual micro-reports disabled will not acquire a
report mutex in event callbacks. If a user explicitly enables micro-reports,
the existing serialized output stream remains protected; this preserves output
integrity while documenting that detailed event reporting is intentionally a
serialization point.

The worker scheduler will allocate phonon indices in chunks rather than
performing one shared atomic decrement before and after every phonon. The
stable mapping from phonon index to random stream remains unchanged. The
global phonon diagnostic counter will be retained initially for compatibility;
its cost will be measured separately before any identifier semantics change.

## Alternatives considered

### Sharded locks

Separate locks for report categories or seismometers would reduce contention
but would still serialize mutable time-bin updates and retain lock overhead on
the common reports-disabled path. It is a lower-risk fallback, not the primary
design.

### Asynchronous report queue

A queue would decouple detailed textual reports from propagation, but it would
not solve shared seismometer-bin accumulation and would add queue lifecycle and
backpressure complexity. It may be considered later for explicitly enabled
micro-report workloads.

## Numerical and behavioral requirements

- `--workers=1` remains the serial reference and keeps existing defaults.
- The same base seed and phonon index produce the same per-phonon random stream
  regardless of worker count.
- Integer counters and invalidity flags must agree between serial and parallel
  runs for the same seeded workload.
- Seismometer energy bins may differ by normal floating-point reduction order
  at the last few bits when multiple workers are used; this is accepted.
- Report line ordering remains unspecified for multiple workers.
- No MPI dependency or C++ standard migration is introduced.
- Worker-local state must be initialized before worker threads start and reduced
  only after all worker threads have joined.

## Testing and acceptance

The implementation is accepted when:

1. focused unit tests cover worker-local counter/bin accumulation and reduction;
2. an end-to-end seeded simulation compares stable summary counters between one
   and multiple workers;
3. `make -j2` and `make test` pass on the branch;
4. the exact 10M-phonon Lop Nor benchmark is rerun for one, two, and four
   workers using the same seed and configuration;
5. the benchmark report records initialization/simulation timing consistently,
   identifies the machine and compiler, and distinguishes numerical agreement
   from wall-clock speedup; and
6. `docs/DEVELOPMENT.md`, `docs/MANUAL.md` if user-visible semantics change,
   and `PLANS.md` describe the new reduction boundary and limitations.
