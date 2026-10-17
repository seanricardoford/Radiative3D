# Anisotropic Scale-Length Performance Measurement

Measured by running `tests/test_anisotropy_unit.sh axi_draws` and
`tests/benchmark_aniso.sh` against the `anisotropy` branch build. All numbers
below come from those runs; nothing here is estimated.

| Date | Host | Compiler | Workload |
|---|---|---|---|
| 2026-10-16 | Apple M2 Ultra, 24 threads | Apple clang version 17.0.0 (clang-1700.6.4.2) | Degree 9; lopnor 1,000,000 phonons; crustpinch 200,000 phonons |

## Benchmark table

`tests/benchmark_aniso.sh` times the isotropic run and the axisymmetric
(`ax=ay=r·a, az=a`) and triaxial (`ax=r·a, ay=√r·a, az=a`) runs at aspect
ratios r ∈ {2, 4, 10}, for both scenarios' do-script ε (0.01) and ε = 0.05.
Anisotropy is applied to the crust and mantle regions for lopnor, and to the
crust, pinched crust, and mantle regions for crustpinch (sedi and moho stay
isotropic in both). Output, verbatim:

| Scenario | eps | Mode | r | Seconds | Slowdown vs isotropic | Gate |
|---|---|---|---|---:|---:|---|
| lopnor | 0.01 | iso | - | 23.01 | 1.00x | ok |
| lopnor | 0.01 | axi | 2 | 109.68 | 4.767x | GATE |
| lopnor | 0.01 | axi | 4 | 109.96 | 4.778x | GATE |
| lopnor | 0.01 | axi | 10 | 111.24 | 4.834x | GATE |
| lopnor | 0.01 | tri | 2 | 27.71 | 1.204x | ok |
| lopnor | 0.01 | tri | 4 | 28.06 | 1.219x | ok |
| lopnor | 0.01 | tri | 10 | 29.41 | 1.278x | ok |
| lopnor | 0.05 | iso | - | 29.09 | 1.00x | ok |
| lopnor | 0.05 | axi | 2 | 122.82 | 4.223x | GATE |
| lopnor | 0.05 | axi | 4 | 134.36 | 4.620x | GATE |
| lopnor | 0.05 | axi | 10 | 183.60 | 6.313x | GATE |
| lopnor | 0.05 | tri | 2 | 38.68 | 1.330x | ok |
| lopnor | 0.05 | tri | 4 | 50.46 | 1.735x | ok |
| lopnor | 0.05 | tri | 10 | 91.30 | 3.139x | GATE |
| crustpinch | 0.01 | iso | - | 8.21 | 1.00x | ok |
| crustpinch | 0.01 | axi | 2 | 29.83 | 3.634x | GATE |
| crustpinch | 0.01 | axi | 4 | 29.81 | 3.632x | GATE |
| crustpinch | 0.01 | axi | 10 | 29.84 | 3.636x | GATE |
| crustpinch | 0.01 | tri | 2 | 9.35 | 1.140x | ok |
| crustpinch | 0.01 | tri | 4 | 9.36 | 1.140x | ok |
| crustpinch | 0.01 | tri | 10 | 9.44 | 1.150x | ok |
| crustpinch | 0.05 | iso | - | 8.26 | 1.00x | ok |
| crustpinch | 0.05 | axi | 2 | 29.91 | 3.619x | GATE |
| crustpinch | 0.05 | axi | 4 | 30.05 | 3.636x | GATE |
| crustpinch | 0.05 | axi | 10 | 30.75 | 3.722x | GATE |
| crustpinch | 0.05 | tri | 2 | 9.49 | 1.149x | ok |
| crustpinch | 0.05 | tri | 4 | 9.77 | 1.182x | ok |
| crustpinch | 0.05 | tri | 10 | 10.61 | 1.284x | ok |

## Isotropic no-slowdown check

Each scenario's default isotropic case was timed three times against the
pre-anisotropy baseline executable and three times against the current
build; the values below are the median current/baseline ratio. Both are
below the 1.05 threshold, so the check passed.

- lopnor: 0.9972
- crustpinch: 0.9962

## Draws per real event

From `build/aniso_draws.txt`, produced by `TEST(axi_draws_per_event)` in
`tests/anisotropy_tests.cpp`. That test runs at take-off-angle degree 6
(`UseTOADegree(6)`), using the Lop Nor crust `ScatterParams` with two
anisotropic specs, `(ax, ay, az)` = (2, 2, 0.5) km and (5, 5, 0.5) km. Each
value is the mean number of tentative-collision draws per accepted
(real) event, averaged over 3000 trials, for the axisymmetric (AXI) and
general (GENERAL) samplers:

### Spec (2, 2, 0.5)

| Ray | theta | chi | AXI draws | GENERAL draws |
|---|---:|---:|---:|---:|
| P | 0 | - | 1.35167 | 9.22567 |
| P | 30 | - | 1.341 | 8.663 |
| P | 60 | - | 1.31333 | 7.42433 |
| P | 89 | - | 1.31067 | 6.638 |
| S | 30 | 0 | 1.70467 | 7.67033 |
| S | 30 | 45 | 1.68033 | 7.59733 |
| S | 30 | 90 | 1.713 | 8.24867 |
| S | 89 | 0 | 1.59033 | 4.95533 |
| S | 89 | 45 | 1.651 | 4.54567 |
| S | 89 | 90 | 1.644 | 4.41467 |

### Spec (5, 5, 0.5)

| Ray | theta | chi | AXI draws | GENERAL draws |
|---|---:|---:|---:|---:|
| P | 0 | - | 2.27167 | 53.7167 |
| P | 30 | - | 2.27133 | 48.794 |
| P | 60 | - | 2.05233 | 33.7237 |
| P | 89 | - | 1.99933 | 22.1213 |
| S | 30 | 0 | 4.41067 | 42.0337 |
| S | 30 | 45 | 4.479 | 44.4553 |
| S | 30 | 90 | 4.59733 | 45.337 |
| S | 89 | 0 | 3.163 | 14.0497 |
| S | 89 | 45 | 3.019 | 13.619 |
| S | 89 | 90 | 3.15733 | 14.1717 |

### Comparison against the spec's prototype table

The design spec ("Axisymmetric sampler") gives a prototype measurement at
fine degree 7, Lop Nor-like crust, kappa 0.5, covering P at 0-89 degrees
incidence and S at 30 and 89 degrees with chi = 0, 45, 90 degrees:

| `ah, az` (km) | General sampler (spec) | Axisymmetric sampler (spec) |
|---|---:|---:|
| 2, 0.5 | 4.5-9.1 | 1.5-1.9 |
| 5, 0.5 | 13-53 | 2.7-5.1 |

The measured ranges above (table rows for theta in {30, 89}, every chi) are,
for comparison, AXI 1.31-1.71 / GENERAL 4.41-9.23 at (2, 0.5), and AXI
1.99-4.60 / GENERAL 13.6-53.7 at (5, 0.5). Both samplers' measured draw
counts agree with the spec's prototype table to within the expected
run-to-run variation (different direction sets, random seeds, and a
degree-6 rather than degree-7 fine tessellation here): AXI 1.3-1.7 against
the spec's 1.5-1.9, and 2.0-4.6 against 2.7-5.1; GENERAL 4.4-9.2 against
4.5-9.1, and 13.6-53.7 against 13-53.

## Measurement gate verdict

**FAILED.** Every axisymmetric row in the benchmark table exceeds the 2×
slowdown gate (lopnor axi: 4.223x–6.313x; crustpinch axi: 3.619x–3.722x), as
does the lopnor ε = 0.05 triaxial r = 10 row (3.139x).

Per the spec's measurement gate, the follow-up for a `GATE` row depends on
which sampler is involved. Neither is implemented here:

- Axisymmetric rows: revisit the axisymmetric sampler's cell and bin
  resolution within its memory budget. The diagnosis below finds the
  axisymmetric gate failures are actually dominated by per-scatterer
  `AxiProposal` construction cost, not sampling/bin-resolution cost, so
  the more direct remedy is parallelizing or restructuring that
  construction across scatterers (the spec's own suggested remedy, "if
  that is material for a model with many axisymmetric layers, compute
  them in parallel across scatterers"); that is a recommended follow-up,
  not implemented in this fix pass.
- The lopnor ε = 0.05 triaxial r = 10 row: consider extending
  incidence-tailored proposals to triaxial P waves, designed separately.

## Diagnosis: fixed construction cost

The axisymmetric times in the benchmark table are nearly flat in r at
ε = 0.01: 109.68 s, 109.96 s, 111.24 s for lopnor (r = 2, 4, 10), and
29.83 s, 29.81 s, 29.84 s for crustpinch. If the slowdown were coming from
sampling cost, it should grow with aspect ratio the way the triaxial rows
do; instead it is essentially constant. That pattern suggests a fixed
per-scatterer construction cost (e.g. building the `AxiProposal` tables
and/or report statistics) rather than sampling cost.

To check this, each scenario was run with `--num-phonons=0` at degree 9
(`TOA_DEGREE=9 tests/run_aniso_case.sh lopnor 24 OUTDIR 0 MODEL_ARGS`),
isolating setup/construction time from any phonon-sampling time, for the
lopnor isotropic default, the lopnor axi r = 2 ε = 0.01 layout, and the
lopnor tri r = 2 ε = 0.01 layout (all from `tests/benchmark_aniso.sh`'s own
model-argument construction):

| Scenario | Zero-phonon seconds | Full-phonon seconds (benchmark table) |
|---|---:|---:|
| iso | 21.947 | 23.01 |
| axi r=2 | 108.125 | 109.68 |
| tri r=2 | 26.467 | 27.71 |

In each case the zero-phonon run accounts for nearly all of the full run's
time (iso: 21.947 of 23.01 s; axi r=2: 108.125 of 109.68 s; tri r=2: 26.467
of 27.71 s) — roughly 95-98%. This is consistent with the flat-in-r pattern
above and suggests that reducing construction cost may matter more than
reducing bin resolution for closing the axisymmetric gate.

### Measured per-scatterer construction cost

To measure the per-scatterer cost directly rather than infer it, the same
zero-phonon lopnor axi r=2 ε=0.01 run (degree 9) was repeated and its
scatterer dump inspected: it builds 20 axisymmetric (`AXI`) scatterers out
of 21 total. Timed against the zero-phonon isotropic lopnor run on the same
host:

| Run | User+sys seconds (zero phonons, degree 9) |
|---|---:|
| iso | 21.70 |
| axi r=2 | 106.11 |

That is (106.11 - 21.70) / 20 ≈ **4.22 seconds per axisymmetric scatterer**
for `AxiProposal` construction plus report-statistics integration combined,
at degree 9. The spec's own estimate was "roughly 1 s per axisymmetric
scatterer"; the measured cost here is about 4x that (consistent with the
review's independent measurement of about 3.8 s/scatterer for `AxiProposal`
alone plus about 0.3 s for report statistics, measured at degree 6 -- this
construction cost does not depend on `--toa-degree`, since the bound-table
build always uses a fixed degree-6 (or lower) set of representative
directions regardless of the run's own take-off-angle degree). This
confirms, with a direct per-scatterer measurement rather than an inference
from the flat-in-r pattern, that per-scatterer construction cost -- not
bin/cell resolution or sampling cost -- is what the axisymmetric gate
failures above come from.

### Status: axisymmetric construction-cost remedy not implemented

The spec names the remedy for this cost directly: "if that is material for
a model with many axisymmetric layers, compute them in parallel across
scatterers." That parallelization (or an equivalent restructuring of the
bound-table build, e.g. hoisting the per-representative invariants noted in
the final review) is a recommended follow-up for closing the axisymmetric
measurement gate. It is not implemented in this fix pass, which is scoped
to correctness and documentation findings only; the axisymmetric sampler's
*results* are unaffected and remain exact (see the sampler-correctness
checks in `tests/anisotropy_tests.cpp`) -- only its one-time, per-scatterer
model-build cost is unoptimized.

