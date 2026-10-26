# TOA Quadrature Runs, 2026-10

Captured output for Task 5 of `.superpowers/sdd/2026-10-25-toa-area-quadrature/`
(scenario-level report and degree-9 build check for `--toa-quadrature`).
Task 6 quotes these files; every number in that documentation task must
trace back to something here.

This capture is the **fix-round-1** recapture: it supersedes an earlier
capture taken at commit `2eb8eea`, before `geom_s2.cpp` reserved the
solid-angle vectors (see "Fix round 1" below) and before
`tests/compare_toa_quadrature.sh` grew its per-station energy-detail
appendix table.

## Provenance

- **Command:** `tests/compare_toa_quadrature.sh > docs/toa-quadrature-runs-2026-10/compare-toa-quadrature.out 2> docs/toa-quadrature-runs-2026-10/compare-toa-quadrature.err`
  (default environment: `TQ_DEGREE=7 TQ_REPLICATES=4 TQ_THREADS=auto TQ_SEED_BASE=596 TQ_LOPNOR_PHONONS=200000 TQ_CRUSTPINCH_PHONONS=50000 TQ_CONFIGS="lopnor:iso lopnor:axi4 crustpinch:iso crustpinch:axi4" TQ_BUILD_CHECK=1`)
- **Date:** 2026-10-26, 11:01:57 to 11:03:14 -0700 (Pacific Daylight Time;
  start/end of the capture, from `date +%s.%N` immediately before/after)
- **Host:** md-5558269, Darwin 24.6.0 arm64, Apple M2 Ultra, 24 logical CPUs
  (`TQ_THREADS=auto` resolved to 24)
- **Git revision:** working tree at commit `2eb8eea` plus this fix round's
  changes to `geom_s2.cpp` and `tests/compare_toa_quadrature.sh` (the
  `make` inside the script built exactly this tree); committed immediately
  after capture as the fix commit this README's own change is part of --
  `git log` on branch `quadrature` shows it as `2eb8eea`'s only child at
  capture time. (`git rev-parse --short HEAD` at the moment of capture
  still read `2eb8eea`, since the fix was not yet committed; citing that
  alone would omit the two files this fix round changed, so it's spelled
  out here instead.)
- **Wall time:** 76.2s (whole-script elapsed time, `date +%s.%N` before/after
  the invocation above; includes a `make` of the current build, 32 energy/MFP
  runs, and 12 degree-9 build-check runs)

## Files

- `compare-toa-quadrature.out` -- the script's stdout: the Markdown energy
  table (summary + per-station detail appendix), MFP table, and degree-9
  build-check table, for the default configuration above.
- `compare-toa-quadrature.err` -- the script's stderr: per-run progress
  messages (config, replicate, seed) for all 44 runs.
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
copying by running `tests/test_anisotropy_unit.sh` fresh (73/73 `PASS`, 0
`FAIL` -- 66 checks in the main batch process plus 7 `isolated_*` checks
each in their own process; `build/anisotropy_tests --list` also reports 73).
Their content is byte-identical to the previous (`2eb8eea`) capture's copies
(confirmed by `md5`): this fix round changes memory behavior and reporting
detail only, never a computed solid-angle, report-MFP, P-frame-MFP, or
isotropic-exactness value.

## Fix round 1 (controller-directed)

Two changes, both documented in full in
`.superpowers/sdd/2026-10-25-toa-area-quadrature/task-5-report.md`'s
"Fix round 1" section:

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

