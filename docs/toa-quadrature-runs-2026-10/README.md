# TOA Quadrature Runs, 2026-10

Captured output for Task 5 of `.superpowers/sdd/2026-10-25-toa-area-quadrature/`
(scenario-level report and degree-9 build check for `--toa-quadrature`).
Task 6 quotes these files; every number in that documentation task must
trace back to something here.

This capture is a **K=12 recapture** (controller-directed, final-review
fix wave): it supersedes the earlier K=4 capture with `TQ_REPLICATES=12`
instead of 4, and a seed scheme that gives every configuration its own
disjoint seed block (`tests/compare_toa_quadrature.sh`'s header comment,
"Seeds") so that, unlike the K=4 capture, no two configurations -- and no
two replicates or modes within one configuration -- ever draw from the
same seed. "Fix round 1" below documents an earlier, separate recapture
(solid-angle-vector reservation and the per-station appendix, taken at
K=4); its own narrative and figures are untouched by this K=12 change.

## Provenance

- **Command:** `tests/compare_toa_quadrature.sh > docs/toa-quadrature-runs-2026-10/compare-toa-quadrature.out 2> docs/toa-quadrature-runs-2026-10/compare-toa-quadrature.err`
  (default environment: `TQ_DEGREE=7 TQ_REPLICATES=12 TQ_THREADS=auto TQ_SEED_BASE=596 TQ_LOPNOR_PHONONS=200000 TQ_CRUSTPINCH_PHONONS=50000 TQ_CONFIGS="lopnor:iso lopnor:axi4 crustpinch:iso crustpinch:axi4" TQ_BUILD_CHECK=1`)
- **Date:** 2026-10-26, 12:50:25 to 12:52:27 -0700 (Pacific Daylight Time;
  start/end of the capture, from `date +%s.%N` immediately before/after)
- **Host:** md-5558269, Darwin 24.6.0 arm64, Apple M2 Ultra, 24 logical CPUs
  (`TQ_THREADS=auto` resolved to 24)
- **Git revision:** working tree at commit `fa33661` plus this final-review
  fix wave's changes to `tests/compare_toa_quadrature.sh` (the
  `TQ_REPLICATES` default and the per-configuration seed-block scheme;
  the `make` inside the script built exactly this tree); committed as
  part of this fix wave's commit(s) immediately after capture -- `git
  log` on branch `quadrature` shows them as `fa33661`'s child(ren) at
  capture time. (`git rev-parse --short HEAD` at the moment of capture
  still read `fa33661`, since the fix wave was not yet committed; citing
  that alone would omit the script change this fix wave made, so it's
  spelled out here instead.)
- **Wall time:** 121.2s (whole-script elapsed time, `date +%s.%N` before/after
  the invocation above; includes a `make` of the current build, 96 energy/MFP
  runs, and 12 degree-9 build-check runs)

## Files

- `compare-toa-quadrature.out` -- the script's stdout: the Markdown energy
  table (summary + per-station detail appendix), MFP table, and degree-9
  build-check table, for the default configuration above.
- `compare-toa-quadrature.err` -- the script's stderr: per-run progress
  messages (config, replicate, seed) for all 108 runs.
- `toa_quadrature_solid_angles.txt` -- copied from `build/` after a fresh
  `tests/test_anisotropy_unit.sh` run (Task 1's `tessel_solid_angles`):
  per-degree solid-angle sum/min/max/rms statistics for the icosahedral
  tessellation, degrees 0-8.
- `toa_quadrature_report.txt` -- copied from `build/` after the same fresh
  run (Task 2's `area_report_directional_mfps`): equal- and area-mode
  mean free paths for the 8 `IntegratedRate`/dipole report states
  (P0/P45/P90/S0/S45SV/S45SH/S90SV/S90SH).
- `toa_quadrature_pframe.txt` -- copied from `build/` after the same fresh
  run (Task 2's `area_p_frame_residual`): equal- and area-mode mean free
  paths for the crust and mantle anisotropic layers at 15-degree P-frame
  rotation steps, plus each layer's exact value.
- `toa_quadrature_iso.txt` -- copied from `build/` after the same fresh run
  (Task 3's `area_iso_mfp_exact`): relative error against the
  exact isotropic MFP, equal vs. area, at degrees 5-7, for two
  (a, kappa) pairs.

All four `toa_quadrature_*.txt` files were regenerated immediately before
this capture by running `tests/test_anisotropy_unit.sh` fresh (73/73
`PASS`, 0 `FAIL` -- 66 checks in the main batch process plus 7
`isolated_*` checks each in their own process; `build/anisotropy_tests
--list` also reports 73). Their content is byte-identical to the
previously-committed copies (confirmed by `md5`): this K=12 recapture
changes the scenario script's replicate count and seed scheme only,
never a computed solid-angle, report-MFP, P-frame-MFP, or
isotropic-exactness value.

## Fix round 1 (controller-directed)

Two changes, both summarized below:

1. **`geom_s2.cpp`**: `Triangle::reserve_lists()` and
   `TesselSphere::init_icosahedron()` now reserve `_SolidAngles`'s capacity
   up front (`pow(4,degree)` and `20*pow(4,degree)` respectively, both
   gated on solid angles being requested at all) instead of letting it grow
   by repeated `insert()`/`push_back()` calls with no reservation -- which
   left old and new backing buffers briefly live together at every
   capacity-doubling reallocation, right when the vector was largest. Pure
   capacity change; no computed value moves (confirmed above).
2. **`tests/compare_toa_quadrature.sh`**: the energy table's summary row
   alone couldn't be interpreted (per-station noise and a genuine
   scenario-wide shift look the same in one median number), so a
   per-station detail appendix and a scenario-wide "all stations" row were
   added; see the script's header comment and the "Energy table:
   per-station detail" section of `compare-toa-quadrature.out`.

**Deterministic storage figure, vs. the noisy measured delta.** The data
the reserve in item 1 is sized for is exact and computable at compile
time, independent of any run: at degree 9 there are 20*4^9 = 5,242,880
directions, and `_SolidAngles` holds one `Real` (`double`, 8 bytes) per
direction, so the data itself is 5,242,880 x 8 bytes = 41,943,040 bytes =
**41.9 MB**. That is now allocated exactly once, by the exact reserve.
Before this fix, `TesselSphere::_SolidAngles`'s unreserved growth (20
`insert()` calls, doubling capacity each time it was exceeded) landed on
a *final* capacity of 8,388,608 entries, not 5,242,880 -- 67.1 MB
steady-state, a 25.2 MB overshoot above the 41.9 MB actually needed --
plus a transient spike at the last reallocation (old 4,194,304-entry /
33.6 MB buffer and new 8,388,608-entry / 67.1 MB buffer briefly coexisting,
for a momentary ~100.7 MB). (Verified directly: a standalone
`std::vector<double>` reproducing that exact 20-insert sequence with this
machine's compiler/standard library confirms this capacity progression --
262144 -> 524288 -> 1048576 -> 2097152 -> 4194304 -> 8388608 -- landing
on size=5,242,880, capacity=8,388,608.) So the fix's theoretical benefit
spans a wide range -- anywhere from the 25.2 MB
steady-state overshoot alone up to the full ~58.7 MB difference between
the old transient peak and the new exact allocation -- depending on
whether "maximum resident set size" happens to catch that momentary
spike.
That range is itself comparable to, or larger than, the build-check
table's actual measured RSS deltas (`compare-toa-quadrature.out`,
"Degree-9 build check" section), which is why those deltas are noisy: at
only 3 replicates per cell, run-to-run allocator noise on a several-GB
process is large enough to swamp a signal in this range. Taking the
pre-fix and post-fix build-check captures together, the measured delta
went the *wrong* way for `iso` after the fix (+68.4 ->
+78.0 MB -- a fix that can only reduce allocation overhead cannot
legitimately increase the steady-state delta) while `axi10` dropped to
+28.6 MB, below even the 41.9 MB exact-data figure. **Task 6 should cite
the 41.9 MB exact-data figure (and, if useful, the 67.1 MB old-capacity /
25.2 MB-overshoot figures above) for the storage cost itself, and should
not treat the build-check table's specific per-variant MB deltas as more
precise than this noise floor allows.**

