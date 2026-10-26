# TOA Quadrature Runs, 2026-10

Captured output for Task 5 of `.superpowers/sdd/2026-10-25-toa-area-quadrature/`
(scenario-level report and degree-9 build check for `--toa-quadrature`).
Task 6 quotes these files; every number in that documentation task must
trace back to something here.

## Provenance

- **Command:** `tests/compare_toa_quadrature.sh > docs/toa-quadrature-runs-2026-10/compare-toa-quadrature.out 2> docs/toa-quadrature-runs-2026-10/compare-toa-quadrature.err`
  (default environment: `TQ_DEGREE=7 TQ_REPLICATES=4 TQ_THREADS=auto TQ_SEED_BASE=596 TQ_LOPNOR_PHONONS=200000 TQ_CRUSTPINCH_PHONONS=50000 TQ_CONFIGS="lopnor:iso lopnor:axi4 crustpinch:iso crustpinch:axi4" TQ_BUILD_CHECK=1`)
- **Date:** 2026-10-26, 02:18:26 to 02:19:42 -0700 (Pacific Daylight Time;
  start/end of the capture, from `date +%s.%N` immediately before/after)
- **Host:** md-5558269, Darwin 24.6.0 arm64, Apple M2 Ultra, 24 logical CPUs
  (`TQ_THREADS=auto` resolved to 24)
- **Git revision:** `2bdaa41` (`git rev-parse --short HEAD`; branch `quadrature`,
  tip of Task 4 at capture time, before this task's own commit)
- **Wall time:** 76.1s (whole-script elapsed time, `date +%s.%N` before/after
  the invocation above; includes a `make` of the current build, 32 energy/MFP
  runs, and 12 degree-9 build-check runs)

## Files

- `compare-toa-quadrature.out` -- the script's stdout: the Markdown energy
  table, MFP table, and degree-9 build-check table, for the default
  configuration above.
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
each in their own process; `build/anisotropy_tests --list` also reports 73),
so they reflect revision `2bdaa41`, the same revision as the scenario
capture above.
