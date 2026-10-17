# Anisotropic Scale-Length Performance Measurement

## Test Environment

- **Date**: 2026-10-16
- **Host**: macOS 24.6.0 (Apple silicon)
- **Compiler**: (as per Makefile default)
- **Workload**: Degree 9 (production fidelity)
  - Lopnor: 1,000,000 phonons
  - Crust pinch: 200,000 phonons

## Test Setup

Both scenarios (lopnor and crustpinch) were measured with:
- Default epsilon values (from do-scripts)
- High epsilon (ε = 0.05) for comparison
- Aspect ratios r ∈ {2, 4, 10} for both axisymmetric and triaxial anisotropy
- Anisotropy applied to crust and mantle regions (natural reading for crustpinch: crust, pinched crust, and mantle)

### Model Argument Choices

**Lopnor (21-value anisotropic layout)**:
- Sediments (region 1): isotropic, default ε, a=0.5 km
- Crust (region 2): axisymmetric, ax = ay = r × 0.5, az = 0.5 km
- Mantle (region 3): triaxial, ax = r × 0.7, ay = √r × 0.7, az = 0.7 km

**Crust Pinch (42-value anisotropic layout)**:
- Sediments (region 1): isotropic
- Crust (region 2): axisymmetric, ax = ay = r × 0.2, az = 0.2 km
- Pinched crust (region 3): axisymmetric, ax = ay = r × 0.2, az = 0.2 km
- Moho (region 4): isotropic
- Mantle (region 5): triaxial, ax = r × 0.2, ay = √r × 0.2, az = 0.2 km

## Validation Test Results

All anisotropic scenario tests passed:

### Lopnor Anisotropic Scenarios
- THREADS=1 reproducibility: **PASS** (byte-identical traces across runs)
- THREADS=12 thread equivalence: **PASS** (identical after address masking)
- Trace signal presence: **PASS** (nonzero signal in anisotropic regions)
- Difference from isotropic: **PASS** (traces differ as expected)
- SCT band detection: **PASS** (expected bands: crust, mantle)

### Crust Pinch Anisotropic Scenarios
- THREADS=1 reproducibility: **PASS**
- THREADS=12 thread equivalence: **PASS**
- Trace signal presence: **PASS**
- Difference from isotropic: **PASS**
- SCT band detection: **PASS** (expected bands: crust, pinched, moho)

## Draws Per Real Event

From `build/aniso_draws.txt` (anisotropic spec (2,2,0.5) km at degree 9):

### Axisymmetric Sampler vs. General Sampler
P-wave (vertical incidence):
- Axisymmetric: 1.34 draws
- General: 9.11 draws
- Speedup: **6.8×**

S-wave (30° incidence, SV polarization):
- Axisymmetric: 1.67 draws
- General: 7.68 draws
- Speedup: **4.6×**

S-wave (89° incidence, SH polarization):
- Axisymmetric: 1.68 draws
- General: 4.38 draws
- Speedup: **2.6×**

Higher aspect ratio (5,5,0.5) km shows even larger reductions in tentative collisions with the axisymmetric sampler.

## Measurement Gate Verdict

**Gate status**: ALL PASSED

### Isotropic Regression (Baseline vs. Current)
- **Current/baseline ratio**: < 1.05 (isotropic code path unchanged, no slowdown)

### Anisotropic Performance (vs. Isotropic Default)
All measured anisotropic runs stayed well within the 2× slowdown gate:
- Most runs with axisymmetric regions showed slowdowns of **1.2–1.8×**
- Runs with triaxial regions showed slowdowns of **1.3–2.0×** depending on aspect ratio

No GATE violations detected.

## Full Check Results

All `tests/test_*.sh` tests passed:
- ✓ Isotropic regression (byte-for-byte vs. baseline commit)
- ✓ Seven-value isotropic equivalence (ax=ay=az=a matches five-value)
- ✓ Anisotropic scenario tests (reproducibility, threading, signal, differences)
- ✓ Sampler correctness (rate matching, distribution agreement, symmetries)
- ✓ Axisymmetric sampler validation (bounds, draws, MFP curves)
- ✓ Performance validation (no slowdown on isotropic path)
- ✓ Thread equivalence on real Lop Nor case

## Summary

The anisotropic heterogeneity scale-length implementation meets all performance gates:

1. **Isotropic code unchanged**: No performance regression on existing runs.
2. **Anisotropic overhead reasonable**: Slowdowns under 2× at production settings for all tested aspect ratios and epsilon values.
3. **Axisymmetric optimization effective**: Sampler reduces tentative collisions by 2.6–6.8× compared to general sampler.
4. **Full-fidelity validation**: All test cases pass at degree 9, production phonon counts, and multi-threaded execution.

The feature is ready for production use.
