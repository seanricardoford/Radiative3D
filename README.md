# Radiative3D
Radiative transport in 3D Earth models

Radiative3D is a 3D radiative transport software tool being developed by Christopher Sanborn and the Solid Earth Geophysics Research Group at the University of Connecticut. Radiative3D can be used to produce synthetic waveforms, travel-time curves, or volumetric visualizations of energy propagation through three-dimensional Earth models. Radiative3D uses ray tracing to simulate propagation dynamics in large-scale structure, and uses a stochastic multiple scattering process to simulate the effects of statistically-described small-scale structure. Radiative3D simulates realistic source events described by moment tensor elements, allowing it to be used to simulate a variety of focal mechanisms, including explosions, double-couple earthquakes, CLVD's, etc.

#### Publications:

* [Modelling Lg blockage, GJI, 2018](https://academic.oup.com/gji/article/214/2/1426/5000173)
* [Simulations with Radiative3D, 2017](https://opencommons.uconn.edu/dissertations/1460/)
* [Combined effects of deterministic and statistical structure, GJI, 2017](https://academic.oup.com/gji/article/210/2/1143/3833065)

## What's New:

##### Janaury 2020:

* Added support for model architecture composed of concentric spherical shells, in which the elastic velocity profiles vary quadraticaly with the radial coordinate (_v<sub>P,S</sub> = a r^2 + c_).  The grid defines the velocities at the top and bottom of each spherical layer. This allows for whole-Earth models, efficiently defined in terms of a reasonably small number of layers.

* Raspbian Buster on the Raspberry Pi 4 is now my primary development and testing environment.  Radiative3D is fast, efficient, and runs on inexpensive hardware!  Of course, it still runs great on Intel and AMD-based workstations as well.

## Build Process

Radiative3D builds with GCC on MacOS (OS X), Linux, and Raspbian.  (And perhaps also Windows.)

```
$ git clone https://github.com/christophersanborn/Radiative3D.git
$ cd Radiative3D
$ make
$ make test
# Optional plotting regression (requires Octave and gnuplot)
$ make test-plotting
```

Results in a binary named `main`.

Run with:

```
$ ./main [args]
```

User Manual here: [Radiative3D Manual Page](docs/MANUAL.md)

There are also supporting scripts (e.g. `do-crustpinch.sh`) to help with managing command line options and organizing the various output files and post-processing of data. The "do-scripts" are written in Bash and may depend on the installation of additional command line tools. Their plotting stages invoke Octave in non-GUI mode with the gnuplot toolkit. (See the [User Manual](docs/MANUAL.md) for details.)

Capability-focused Lop Nor recipes are also included:

```bash
./do-lopnor-big.sh big-demo
./do-lopnor-big-parallel.sh parallel-demo
./do-lopnor-anisotropic.sh anisotropic-demo
./do-lopnor-anisotropic-equal.sh lopnor-anisotropic-equal
```

The big recipe is a serial, fully visualized `10M`-phonon baseline. The
parallel recipe uses the same `10M` workload, model, and fixed seed with eight
explicit shared-memory workers, making it suitable for a speed comparison.
The anisotropic recipe follows the standard `do-lopnor.sh` waveform workflow
with global horizontal and vertical scattering correlation lengths of `1.25`
and `0.625`, respectively, in the model's distance unit.
The equal-length recipe uses `1.25` for both directions as an isotropic-limit
runtime comparison.

## Current capabilities

- Shared-memory simulation workers are selected with `--workers=N` or
  `--threads=N`; use `--workers=1` for a serial reference run.
- Reproducible per-phonon random streams are selected with `--seed=VALUE`.
- Ellipsoidal anisotropic scattering is configured with paired
  `--scatter-horizontal=L_H` and `--scatter-vertical=L_V` options.
- The focused regression suite is available through `make test`.
- Plotting compatibility can be checked with `make test-plotting` when GNU
  Octave and gnuplot are installed.

## Development handoff

The current development state, architecture, reproducibility rules, testing
workflow, and follow-up roadmap are documented in
[`docs/DEVELOPMENT.md`](docs/DEVELOPMENT.md), [`AGENTS.md`](AGENTS.md), and
[`PLANS.md`](PLANS.md).
