### NAME
Radiative3D - A program that uses radiative transport theory to model seismic energy propagation in 3D Earth models.

### SYNOPSIS
When run directly, Radiative3D is invoked as follows:

    ./main [OPTIONS]...

Often, however, Radiative3D runs are scripted, and a few example "do-scripts" are included as samples.  The scripts: (1) set up model and simulation parameters and assemble them as command line options, (2) start the simulation in a timestamped output directory, (3) record program output and other information into log files, and (4) run visualization scripts after run completion to produces figures or videos. One such example do-script is <tt>do-crustpinch.sh</tt>.  It is recommended to copy the script and modify to taste.  It is run in the following way:

    ./do-crustpinch.sh [RUN_ID]

where _`RUN_ID`_ is an optional alphanumeric identifier that will become part of the name of the output directory. 

### DESCRIPTION
Radiative3D models energy propagation through both deterministic and statistical structure which are specified and characterized on a gridded dataset referred to as the "Earth model".  Arguments to Radiative3D control various aspects of the operation, including specification of the Earth model, specification of the broader "physics model", and the selection and control of various inputs and outputs.

#### Arguments

* <tt>-F, --frequency=_FREQ_</tt>

  Sets the frequency of all generated phonons.  Radiative3D operates at a single frequency.  To generate broadband results, do multiple data runs and combine the results.

* <tt>-N, --num-phonons=_NUM_PHONONS_</tt>

  Number of phonons to emit from source before exiting.  Can take a multiplier suffix of `K`, `M`, or `B` for thousand, million, and billion, respectively.  E.g., use `--num-phonons=10M` for ten million phonons.

* <tt>-T, --timetolive=_time_</tt>

  Sets the time to live for each phonon, or in other words the amount of sim-time to model.

* <tt>--threads=_N_</tt>

  Sets the number of phonon-simulation worker threads. If omitted, Radiative3D
  selects a worker count automatically from the available hardware and the
  requested phonon count. Use <tt>--threads=1</tt> for a single-thread
  debugging or reference run. Explicit values must be positive.

* <tt>--seed=_N_</tt>

  Sets the non-negative random seed used to derive independent phonon random
  streams. If omitted, a seed is generated automatically and printed in the
  simulation parameters so the run can be repeated. Event-report ordering is
  not guaranteed when multiple workers are active.

  Phonon <tt>i</tt>'s random stream is seeded from <tt>seed XOR i</tt>. Under
  seeds <tt>s</tt> and <tt>s XOR d</tt>, phonon <tt>i</tt>'s stream equals
  phonon <tt>i XOR d</tt>'s stream under the other seed, so when <tt>d</tt>
  is smaller than the phonon count the two runs share all, or all but a few,
  of their streams (all of them when the phonon count is a multiple of a
  power of two larger than <tt>d</tt>) -- seeds 596 and 597 (<tt>d = 1</tt>),
  for example. For independent replicate runs, use seeds that differ above
  bit 32, e.g. <tt>596 + k&middot;2^32</tt>.

### Multithreaded performance benchmark

The simulation-only Lop Nor benchmark compares a serial reference run, an
explicit 12-worker run, and automatic worker selection. It does not invoke
Octave or generate figures:

    LOPNOR_BENCH_PHONONS=210000000 \
    LOPNOR_BENCH_OUTPUT_DIR=/path/to/benchmark-output \
    tests/benchmark_lopnor_threads.sh

The benchmark uses the same fixed-seed scenario for all three modes. Record
the host and compiler with the results because performance depends on the
machine:

| Date | Host | Compiler | Workload |
|---|---|---|---|
| 2026-10-16 | macOS host, 24 logical CPUs | Apple Clang/G++ C++11 build | 210,000,000 phonons |

| Mode | Workers | Phonons | Seconds | Phonons/s | Speedup |
|---|---:|---:|---:|---:|---:|
| Serial | 1 | 210,000,000 | 221.917 | 946,300 | 1.00x |
| Explicit | 12 | 210,000,000 | 27.623 | 7,602,360 | 8.03x |
| Automatic | 15 | 210,000,000 | 28.122 | 7,467,463 | 7.89x |

The benchmark was rerun after removing unnecessary report-mutex acquisition
when reporting is disabled. The original implementation serialized these
frequent report calls even with <tt>--reports=ALL_OFF</tt>; the corrected
results show the expected strong scaling up to the useful core range of this
host.
Automatic selection uses the 16 performance logical CPUs reported by macOS and
reserves one; the explicit 12-worker run was slightly faster in this sample.

* <tt>-A, --toa-degree=_TOA_DEGREE_</tt>

  Sets the take-off angle degree, which determines the number of discrete [[take-off angle]]s that PhononSource objects (e.g. event sources and scattering sources) have to choose from.  Increasing this number increases the angular fineness of these objects.  "Degree," in this context, is an integer that represents the number of times the angular space is subdivided, and in the default (and at present only) method for discretizing take-off angles, the algorithm starts with 20 angles and subdivides by four with each iteration (degree), such that N_TOA = 20 * 4 ^ (degree).

* <tt>--toa-quadrature=_MODE_</tt>

  Selects how the take-off-angle set is weighted when building scatterer and event-source probability tables, mean free paths, and dipole moments. Two modes are recognized:

  * <tt>--toa-quadrature=equal</tt> _(default)_

    Every take-off direction carries the same weight, exactly as in all earlier versions.

  * <tt>--toa-quadrature=area</tt>

    Each direction is weighted by the solid angle of the icosahedral triangle it represents, <tt>w = Ω·N/4π</tt> (Ω the triangle's solid angle, N the number of directions in the set being weighted). This mode changes every output that depends on the take-off-angle set, including isotropic regions, and prints three additional lines identifying it: <tt>TOA Quadrature: area</tt> in the parameter echo, a line reading <tt>(Area-weighted take-off-angle quadrature.)</tt> under the take-off-angle count, and <tt>#  TOA quadrature: area</tt> in the scatterer-dump header (see _Scatterer dump_, below). See _Validity limits_, below, for its effect on reported and simulated values.

  <tt>--event-test</tt> ignores this option: it always displays per-direction pattern values under <tt>equal</tt> weighting, since a weighted value would distort that display.

##### _Grid args:_

* <tt>--grid-file=_FILENAME_</tt> _(not yet implemented)_

  This specifies that the model grid should be read from file _`FILENAME`_ and used for the construction of the Earth model.  The grid file format is described in: _[[Grid formats in Radiative3D]]_. (Note: at the current time, reading of grid files is not yet supported, and users should hard-code models. See the `--grid-compiled` option.)

* <tt>--grid-compiled[=_INDEX_]</tt>

  Specifies that the model grid to be used is compiled-in and defined in the function <tt>Grid::ConstructGridManual()</tt> located in the user-modifiable source file <tt>user.cpp</tt>.  Defining a model in source code allows full access to the [[grid building API]] and can facilitate construction of quite complex models, in comparison to models defined via [[Grid formats in Radiative3D|grid files]] (c.f. <tt>--grid-file</tt> option).  It can also enable aspects of the model to be parameterized by a set of model arguments provided by the <tt>--model-args</tt> option in order to facilitate scripting of structural variations in a model.  The optional integer-valued _`INDEX`_ parameter, if included, will be passed to <tt>Grid::ConstructGridManual()</tt> and can be used as a selection index if the user wishes to define multiple user-coded models and choose between them at runtime without need of a recompile.

* <tt>--dump-grid</tt>

  If specified, a parsable plaintext dump of the grid will be written to <tt>stdout</tt> prior to simulation start.

##### _Model args:_

* <tt>--model-args=_ARG_LIST_</tt> <br>
  <tt>--model-compiled-args=_ARG_LIST_</tt> _(deprecated form)_

  Some model grids allow parameterized values, enabled structural variations of the model to be easily scripted. This option allows a comma-separated list of double-precision values to be passed in on the command line and made available as parameters to these grids.

* <tt>--range=_RANGE_</tt> <br>
  <tt>--cylinder-range=_RANGE_</tt> _(deprecated form)_

  When using a LAYERED model defined by a plumb-line grid, the model cells are horizontal layers separated by planar tilted interfaces.  The otherwise infinite lateral extent of these layers is truncated by imposing a cylindrical bounding radius centered on a vertical through the origin.  This option sets the bounding radius, effectively defining the distance range that is modeled.  This option only applies to LAYERED models and is ignored for TETRA models based on tetrahedral model cells.

* <tt>--flatten</tt>  *** _DEPRECATED_ ***

  Apply an [[Earth-flattening transformation]] to depth coordinates and velocities when constructing the model.  This is intended to enhance ray-turning in layered Earth models as a proxy for the spherical curvature that is lacking from those models.  See [[Aki and Richards 1980 - Quantitative Seismology|Aki and Richards]] for more info on EFT's. (Note: the velocity transform is performed only for layered model types, but the depth transform occurs for both layered or tetrahedral based models, though it should be noted it is generally unnecessary for tetra models, as the model mesh can be constructed to follow real Earth curvature.)

  Deprecation note:  the <tt>--flatten</tt> option is deprecated, as coordinate mapping choices, including flattening, are now attributes of the Grid object, and as such are specified in the grid definition file, or via Grid class methods if hard-coding a grid.  If a mapping choice is specified in this way, then the <tt>--flatten</tt> option will have no effect and will be silently ignored.  This option may be removed in a future release.

##### _Sim args:_

* <tt>--overridemfp=_P_,_S_</tt>

  If specified, use the provided mean free path values for scattering model-wide, despite what the calculated mean free paths would be from the otherwise specified heterogeneity parameters.  Scattering shapes (deflection profiles) are unaffected by this argument.

  For an anisotropic (axisymmetric or triaxial) region, phonons still draw path lengths directly from the override mean free path, with no null collisions; at each event the sampler redraws a proposal (conversion type and take-off angle) until one is accepted, so the deflection pattern still follows the region's direction- and polarization-dependent anisotropic rates, not a uniform distribution. Combining <tt>--overridemfp</tt> with an anisotropic region whose <tt>eps = 0</tt> fails at model-construction time with an explicit error (every scattering probability in such a region is identically zero, so the redraw-until-accepted loop above would otherwise either divide zero by zero or spin forever); use a positive <tt>eps</tt>, or an isotropic region, instead.

* <tt>--nodeflect</tt>

  If specified, scattering deflections are "squashed," meaning that while scattering events still occur, they do not result in any deflection or other modifications (raytype or polarization changes) to the phonon's ray trajectory.  This is commonly used, sometimes in concert with <tt>--overridemfp</tt> for producing videos (which plot scattering events), to easily visualize "clean" (i.e. no scattering) wavefront propagation. Or in other words, this option changes scattering events into mere "checkpoint" events used to illustrate evolving wavefronts.

  For an anisotropic region, delta tracking (null-collision thinning) still runs as normal, so the rate of (non-deflecting) scattering events still reflects the region's direction- and polarization-dependent anisotropic rates; only the deflection itself is suppressed. With both <tt>--overridemfp</tt> and <tt>--nodeflect</tt> given together, an anisotropic region collapses to the same single dummy isotropic scatterer it would without anisotropy, as before.

##### _Event args:_

* <tt>-L, --source-loc=_X_,_Y_,_Z_</tt>

  Specifies the location of the event source in cartesian X,Y,Z coordinates.  In the future, other coordinate systems may be supported (such as Lat, Lon, and depth).

* <tt>-E, --source=_SOURCE-ARGS_</tt>

  Determines the event source type via the provided source arguments.  The format of the source args is a keyword followed by a comma-sepated (no whitespace) list of arguments that are required by the given keyword.  The keywords are as follows:

  * <tt>--source=EQ</tt>

    Specifies a "basic earthquake" source.  The source type will be a simple double-couple with principle axes in the X-Y plane.  (Ie, energy radiates predominantly in the X-Y directions, with little energy radiated virtically from the source.)

  * <tt>--source=EXPL</tt>

    Specifies an "explosion" source, which radiates equal energy in all directions (isotropic energy release).

  * <tt>--source=USGS,_rr_,_tt_,_pp_,_rt_,_rp_,_tp_</tt>

    Specifies the full moment-tensor using the USGS standard of specifying elements in a reference with axes pointing UP, SOUTH, and EAST (NOTE: NEED TO FIGURE OUT WHAT THE DIRECTIONS ACTUALLY ARE...)  This method of specifying depends on a specified Earth-center and polar-axis being defined—these parameters are set by the <tt>--earth-frame</tt> argument.

  * <tt>--source=SDR,_strike_,_dip_,_rake_,_[isofrac|isoangle_,_[moment]]_</tt>

    Specifies a double-couple [[moment tensor]] with an optional isotropic component added in.  The double-couple orientation is given in terms of _`strike`_, _`dip`_, and _`rake`_ angles, given in degrees.  The optional fifth argument, _`moment`_, sets the magnitude of the moment tensor, which otherwise defaults to 1.0.  This parameter has no bearing on the actual simulation, but does affect the numbers displayed when the moment tensor is output to the console (the parameter has purely aesthetic utility).  The optional fourth argument, _`isofrac`_, specifies the isotropic component of the moment tensor as a fraction of the total squared-magnitude of the moment tensor. Valid values are in the range [-1.0 to 1.0], with 0.0 being the default.  Positive values indicate an explosive component, whereas negative values indicate an implosive component. Zero means pure double-couple, or no isotropic component.  This argument can alternatively be specified as an angle, rather than a fraction, owing to the fact that [[isotropic]] and [[deviatoric]] moments form orthogonal subspaces of the total moment tensor vector space.  If this manner of specification is desired, the value _`isoangle`_ should be given in the split range of [-90, -1) or (1, +90].  (The range between -1 and 1 is inaccessible because it will be interpreted as _`isofrac`_ rather than _`isoangle`_, but this is unlikely to be problematic in practice.)

##### _Reporting / data output args:_

* <tt>--reports[=_REPORTS_FLAGS_]</tt>

  Enables real-time event reporting.  These are reports of various simulation events detailing the trajectories of individual phonons, such as generation at source, reflection off interfaces, scattering events, collection by seismometers, etc., written to an output stream as they are simulated.  The aggregate of these micro-reports can be used for, among other things, producing videos of energy propagation within the Earth model.  Multiple keywords can be provided as a comma-separated list in order to specify which event types are desired. Event reports are written to <tt>stdout</tt> unless directed to a file by the <tt>--report-file</tt> option.  The keywords are:

  * <tt>--reports</tt>, <tt>--reports=ALL_ON</tt>

    Turns ON all real-time reports. (If called without any keyword, <tt>ALL_ON</tt> is assumed.)
    
  * <tt>--reports=ALL_OFF</tt>

    Turns OFF all real-time reports. _(DEFAULT)_

  * <tt>--reports=GEN,SCT,REF,CEL,COL,LST,TMO</tt>

    This shows all event types explicitly turned on.  They are, respectively: Generate, scatter, reflection, cell-to-cell transfer, collection (interaction with a surface where a seismometer may reside), phonon lost (departing boundary of model), and phonon timeout (phonon lifetime exceeds window of interest as determined by the <tt>--timetolive</tt> option).

  * <tt>--reports=SCATTERS</tt>

    This is a shortcut equivalent to <tt>--reports=GEN,SCT,REF</tt> to report generate, scatter, and reflection events.  Useful for visualizing bulk energy transport, e.g. in video form.

* <tt>--report-file=_filename.dat_</tt>

  Specifies the name of the file that real-time event reports will be written to.  Without this option, reports go to <tt>stdout</tt>.  If option is set, reports go only to the file and will not be printed to <tt>stdout</tt>.

* <tt>--output-dir=_path_</tt>

  Specifies a directory name in which in which all files created by Radiative3D will be created.  If this option is not specified, then all files are relative to the current working directory.

* <tt>--mparams-outfile=_filename.octv_</tt>

  If this is specified, various model and mission parameters will be written to _filename.octv_ in a format readable by GNU Octave for easy processing.

* <tt>--summary=_SUMMARY_FLAGS_</tt>

  Enables or disables various post-simulation summary reports.  These include such things as the trace data from the various seismometers that were requested. Multiple keywords can be provided as a comma-separated list.  The keywords are:

  * <tt>--summary=ON</tt>

    Turns ON all summary reports

  * <tt>--summary=OFF</tt>

    Turns OFF all summary reports

##### _Seismometer args:_

* <tt>--binsize=_number_</tt>

  When binning energy for seismogram output, this is the width of each bin in time units (typically seconds).  (Mutually exclusive with <tt>--binspercycle</tt>.)

* <tt>--binspercycle=_number_</tt>

  An alternate way to specify bin size.  This lets you specify the number of time bins per cycle at the selected frequency.  Can be fractional.  E.g. 1.5 time bins per cycle means three bins will cover two cycles worth of time.  (Mutually exclusive with <tt>--binsize</tt>.)

* <tt>--seismometer=_SEIS_ARGS_</tt>

  Adds a seismometer according to _SEIS_ARGS_.  This flag can be passed multiple times.

* <tt>--seis-array=_AZ,R0,gap,gather,#seismometers_</tt>

  Adds an array of seismometers according to _SEIS_ARRAY_ARGS_.  This flag can be passed multiple times.

* <tt>--seis-p2p=_X1,Y1,Z1,X2,Y2,Z2,R0,gather,num_seismometers_</tt><br>
  <tt>--seis-p2p=_X1,Y1,Z1,X2,Y2,Z2,R0,gather1,gather2,num_seismometers_</tt><br>
  <tt>--seis-p2pw=_X1,Y1,Z1,X2,Y2,Z2,R0,gather,num_seismometers_</tt><br>
  <tt>--seis-p2pw=_X1,Y1,Z1,X2,Y2,Z2,R0,gather1,gather2,num_seismometers_</tt><br>

  Adds a "point-to-point" array along the line (X1,Y1,Z1) to (X2,Y2,Z2).  If two gather radii are specified, then the gather radius increases linearly from _gather1_ at the beginning point to _gather2_ at the end point.  If the option is specified as <tt>--seis-p2pw</tt>, then the gather radii are specified in elastic wavelengths, which are a function of elastic velocity at the seismometer's emplacement location at the specified frequency (and will result in differing gather radii for P and S wave modes).  Otherwise, the gather radii are specified in the model's chosen length unit, typically kilometers.  (Note however that, at present, the <tt>--seis-p2p</tt> option, if specified with only a single gather radius, still assumes this radius is in wavelengths.  This is for backwords compatibility with old scripts, but this usage is deprecated and will be changed in a future version.)

  This option can be specified multiple times.

#### Earth models

Earth models are constructed as a mesh of nodal locations at which known or hypothesized material properties are specified. Based on these nodal locations, Radiative3D divides space into a set of model cells inside which material properties are interpolated based on the values specified at the nodes.  The material properties specified include the elastic properties of _Vp_, _Vs_, _density_, and _Q_, and four additional properties defining the [[heterogeneity spectrum]] of the material: _nu_, _epsilon_, _a_, and _kappa_.

##### _Heterogeneity parameters and modes:_

Each region's heterogeneity spectrum (a von Karman autocorrelation function, following Sato and Fehler) is given by:

- **nu** -- the coupling between the fractional velocity and density perturbations (Sato and Fehler's Eq. 4.48).
- **eps** -- the RMS fractional velocity perturbation of the heterogeneity.
- **a** (or, anisotropically, **ax**, **ay**, **az**) -- the correlation (scale) length(s) of the heterogeneity, in the model's chosen length unit.
- **kappa** -- the von Karman parameter (Hurst exponent) controlling the spectrum's high-wavenumber falloff.
- **Q** -- the intrinsic (thermal) attenuation quality factor for the region.

Scale lengths may be given either as a single isotropic value <tt>a</tt>, or as three separate values <tt>ax, ay, az</tt> along the region's **fabric frame**, an orthogonal triad of local East (<tt>x</tt>), local North (<tt>y</tt>), and local Up (<tt>z</tt>) directions at the scattering point. In flat and flattened coordinate mappings the fabric frame coincides everywhere with the model's internal <tt>+x, +y, +z</tt> axes. In curved mappings, the fabric frame instead follows the local vertical, so that horizontally-stretched heterogeneity remains horizontal (rather than tilting with the internal axes) as it is carried across the curvature of the Earth.

Radiative3D auto-detects the heterogeneity mode of a region from its scale lengths, comparing <tt>ax</tt>, <tt>ay</tt>, and <tt>az</tt> exactly:

- **isotropic** -- <tt>ax == ay == az</tt>. Scattering uses the same code path and random-draw sequence as before anisotropic scale lengths were introduced.
- **axisymmetric** -- <tt>ax == ay != az</tt> (i.e., symmetric about the vertical axis, but distinct from horizontal). Scattering is sampled from the axisymmetric, incidence-tailored sampler.
- **general (triaxial)** -- any other combination of distinct <tt>ax, ay, az</tt>. Scattering is sampled from the general anisotropic sampler.

All three scale lengths must be strictly positive. (This is a behavior change for the five-value layout: a region specified with <tt>a &lt;= 0</tt>, previously accepted, now fails construction with the "must be positive" message above. Scripts relying on the old behavior should use a small <tt>eps</tt> instead of a non-positive <tt>a</tt> -- an <tt>eps = 0</tt> region scatters at zero rate, same as before, without needing a workaround scale length.)

##### _Seven-value (anisotropic) model-args layouts:_

Only three of the included example models accept the seven-value layout: Lop Nor Baseline, Lop Nor Moho, and Crust Pinch. Every other compiled-in model (Upthrust, Halfspace, SphereEarth, ToySphere, Lop Nor Moho2, and ScatParamsStudy) remains isotropic-only and still takes only the five-value layout, though hand-coded C++ models built directly against the <tt>Elastic::HetSpec</tt>/<tt>HSnexyzk</tt> API can still construct anisotropic specs. For the three models that do accept it, the seven-value layout gives <tt>nu, eps, ax, ay, az, kappa, Q</tt> per region, in place of the usual five, <tt>nu, eps, a, kappa, Q</tt>, to specify per-axis heterogeneity scale lengths. Every region of a given model invocation uses the same layout width (all five-value, or all seven-value); the total argument count selects both the model variant and the layout width:

| Model | Regions | Five-value count | Seven-value count |
|-------|--------:|------------------:|--------------------:|
| Lop Nor Baseline | 3 | 15 | 21 |
| Lop Nor Moho | 4 | 20 | 28 |
| Crust Pinch | 5 | 25 | 35 |
| Crust Pinch + layer thicknesses (<tt>SediThick, CrustThick, MohoThick</tt>) | 5 | 28 | 38 |
| Crust Pinch + layer thicknesses + pinch fractions (<tt>PinchFrac, DepthFrac, SediFrac, MohoFrac</tt>) | 5 | 32 | 42 |

A seven-value layout with <tt>ax = ay = az = a</tt> for every region is isotropic and produces output byte-identical to the corresponding five-value layout. Grid dumps (<tt>--dump-grid</tt>) for a model containing any non-isotropic region add two extra columns, <tt>ay</tt> and <tt>az</tt>, following <tt>kappa</tt> (every row, including <tt>***</tt> missing-data rows, gets both columns, so the column count stays consistent for downstream tools such as Octave's <tt>load</tt>).

##### _Scatterer dump:_

On every run, Radiative3D prints a scatterer dump (a table of the distinct <tt>Scatterer</tt> objects built for the model, headed <tt>BEGIN SCATTERER DUMP</tt>) before simulation begins. Each row's first 11 columns (<tt>nu, eps, a, kappa, el, gam0, address, MFP (P), MFP (S), DM (P), DM (S)</tt>) are unchanged from before anisotropic scale lengths were introduced; the <tt>a</tt> column holds <tt>ax</tt> (identical to <tt>a</tt> for an isotropic region).

If a model contains any anisotropic region, the column-header line gains the text <tt>ay az mode directional-MFPs</tt>, and an anisotropic row (one built from a non-isotropic <tt>ScatterParams</tt>) appends, after its existing 11 columns:

- <tt>ay</tt>, <tt>az</tt> -- the region's remaining two scale lengths.
- A mode label, <tt>AXI</tt> for an axisymmetric spec (<tt>ax == ay != az</tt>) or <tt>ANI</tt> for any other (triaxial) anisotropic spec.
- The row's directional mean free paths, computed by direct numerical integration of the anisotropic scattering coefficients over a fixed degree-6 direction set (81,920 directions, independent of the run's <tt>--toa-degree</tt>) at model build:
  - For an <tt>AXI</tt> row, 8 values: P at incidence angles of 0, 45, and 90 degrees from vertical (the <tt>az</tt> axis), then S at 0 degrees (horizontally polarized -- SV and SH coincide there by the medium's symmetry about vertical), and SV and SH separately at 45 and at 90 degrees.
  - For an <tt>ANI</tt> row, 9 values: P travelling along East, North, and Up, then S travelling along each of those three axes, polarized along each of the other two (e.g. S travelling East polarized Up is SV, polarized North is SH).

For an anisotropic row, the leading <tt>MFP (P)</tt>, <tt>MFP (S)</tt>, <tt>DM (P)</tt>, and <tt>DM (S)</tt> columns are also computed by this same direct integration, independent of mode: <tt>MFP (P)</tt> is 3 divided by the sum of the three P-axis rates above (East, North, Up), <tt>MFP (S)</tt> is 6 divided by the sum of the six S-axis rates, and the <tt>DM</tt> columns are the plain mean of the corresponding dipoles. These values describe the underlying physical scattering pattern; they are distinct from the majorant mean free path the general sampler uses internally to draw tentative collision distances. A model with only isotropic regions prints exactly today's 11-column rows and header, with no appended text.

Every value in this dump follows the run's selected take-off-angle quadrature (<tt>--toa-quadrature</tt>), and in <tt>area</tt> mode the dump also gains a <tt>#  TOA quadrature: area</tt> header line, printed once immediately before the column-header line above.

##### _Validity limits:_

The Born/Bourret scattering approximation that Radiative3D's transport theory relies on requires <tt>k^2 * a_max^2 * eps^2 &lt;&lt; 1</tt>, where <tt>k</tt> is the background wavenumber and <tt>a_max</tt> is the largest of a region's scale lengths. Model parameters should be chosen so that this quantity stays well below one. Additionally, Margerin (2006) cautions that for extremely flat scatterers (<tt>az -&gt; 0</tt> with the horizontal scale lengths remaining large), the medium behaves like fine layering, a regime in which radiative transport theory itself breaks down; such parameter choices should be avoided.

By default (<tt>--toa-quadrature=equal</tt>), reported and simulated anisotropic rates (the scatterer-dump directional mean free paths, and the rates the simulation actually realizes) both use an equal-weight average over the run's icosahedral take-off-angle set -- the same convention <tt>ComputeMFPs()</tt> already uses for isotropic regions. That convention is not bias-free there either: because the icosahedral triangles' solid angles are not all equal, the equal-weight average introduces a direction-independent offset of about 1% (e.g. +1.1% for Lop Nor crust S) rather than none, as detailed in <tt>docs/toa-quadrature-bias.md</tt>. The same non-uniformity also distorts which directions carry that average: under <tt>equal</tt> weighting, a node's effective share of sampled or emitted phonons is modulated by mean(Ω)/Ω_k relative to a uniform density on the sphere, a factor that ranges from 0.829 to 1.078 across the icosahedral node set at degrees 6-8 (<tt>docs/toa-quadrature-runs-2026-10/toa_quadrature_solid_angles.txt</tt>). The <tt>--toa-quadrature=area</tt> option (above) weights each direction by its own Ω_k instead, removing both effects: for an isotropic region with <tt>a</tt> = 0.5, <tt>kappa</tt> = 0.3, equal weighting's directional-MFP error against the exact integral (-0.04% P, +1.09% S at take-off-angle degree 5) falls to about -0.00005% (P) and +0.0005% (S) under area weighting (<tt>docs/toa-quadrature-runs-2026-10/toa_quadrature_iso.txt</tt>). For an anisotropic region the bias is no longer direction-independent: the icosahedral node density is non-uniform and does not tend to a uniform limit as the degree increases, so a direction- and polarization-dependent quadrature bias of up to several percent persists at every degree. For example, for a <tt>(2, 2, 0.5)</tt> km Lop Nor-like crust, the horizontal-incidence, horizontally-polarized (SH) S mean free path converges (confirmed at <tt>--toa-degree</tt> 6, 7, and 8) to about 398 km under this quadrature, against an exact numerical integral of about 375 km -- roughly 6%, which is a sizeable fraction of the anisotropy signal at that aspect ratio. Under <tt>--toa-quadrature=area</tt> this same state converges instead to that exact integral itself, about 374.7 km, within 2e-4 relative at take-off-angle degree 6 (<tt>docs/toa-quadrature-runs-2026-10/toa_quadrature_report.txt</tt>). This does not violate exactness: the sampler realizes exactly the rate the report prints, so simulated and reported rates always agree with each other; it is a property of the fixed take-off-angle discretization, not a defect in the sampler. For anisotropic P incidence specifically, there is no physical polarization-rotation angle to measure this bias against: the sampler must still pick some frame to build its take-off-angle grid in, and the equal-weight quadrature makes the realized rate depend on that otherwise-meaningless choice -- for a Lop Nor r=10-like crust and mantle, scanning the frame angle about the ray spans 1442-1518 km at 60 degrees from vertical and 742-826 km at 89 degrees, respectively, unchanged at <tt>--toa-degree</tt> 5, 7, and 8 (see <tt>docs/toa-quadrature-bias.md</tt>). Under <tt>--toa-quadrature=area</tt> this dependence does not disappear, but shrinks to a residual that falls roughly 4x with each additional take-off-angle degree (degrees 5-7, the design spec's feasibility-spike table; see <tt>docs/toa-quadrature-bias.md</tt>'s Frame orientation section): at degree 6 the measured spread across all 24 sampled frame angles is about 0.013% for the crust case and 0.028% for the mantle case, each within 2e-4 of the frame-independent exact rate (<tt>docs/toa-quadrature-runs-2026-10/toa_quadrature_pframe.txt</tt>). Both the axisymmetric and general samplers use the canonical (highest-rate) orientation for P incidence, so their results always agree with each other, including for an already-scattered P phonon carrying an arbitrary polarization-rotation angle from its previous scattering event. The <tt>--toa-quadrature=area</tt> option above implements an area-weighted take-off-angle quadrature, applied consistently to both the isotropic and anisotropic paths and to the event source; it does not remove the P-incidence frame-orientation dependence outright, but reduces it to the shrinking residual described above -- the canonical-orientation convention in the previous sentence keeps both samplers in agreement in either mode.

An axisymmetric region's proposal tables (built once per distinct axisymmetric <tt>Scatterer</tt>, not per phonon) are now built in parallel, using the run's thread count (<tt>--threads</tt>, or the automatic worker limit when it is omitted); construction draws no random numbers and produces identical tables at any thread count. This construction carries a fixed per-scatterer cost that does <i>not</i> shrink as <tt>--toa-degree</tt> grows -- the bound-table build always uses degree-6-or-lower representative directions regardless of the run's own degree, the same cap the directional-MFP integration above uses -- roughly 25-30 ms at 24 threads, roughly 0.35 s single-threaded; the general (triaxial) sampler's own isotropic-table build carries no such fixed term and does shrink with degree, so which sampler builds faster depends on <tt>--toa-degree</tt>, thread count, and how many axisymmetric regions a model has. At take-off-angle degree 9 with 24 threads this parallel construction no longer dominates model-build time: the full Lop Nor model (20 axisymmetric regions) now builds in about 4.5 seconds total, against about 3.9 seconds for the same medium built with the general sampler instead -- see <tt>docs/anisotropy-performance.md</tt> for the measurement and its build-time history. At this program's <i>default</i> <tt>--toa-degree</tt> of 7, the same fixed per-scatterer term is large enough, relative to a smaller build, that a run using the general sampler can instead finish faster overall: a scenario-level measurement at degree 7 found the axisymmetric sampler's total run time 1.24-1.51x slower than the general sampler's in half of the configurations tested -- short, construction-dominated runs, where this fixed build cost is a larger fraction of the total (<tt>docs/anisotropy-performance.md</tt>). This is a one-time cost at model build, not a per-phonon cost, and does not affect simulated run time; it is negligible next to axisymmetric <i>sampling</i>'s own cost on a scattering-heavy run, where the axisymmetric sampler runs up to about 2.1x faster than the general sampler on the same medium (same reference) -- so the practical guidance is that the axisymmetric sampler is faster overall on scattering-heavy runs (high <tt>eps</tt>, large aspect ratio, many phonons), and can cost up to a few tenths of a second more than the general sampler to build on a short or lightly-scattering run at this program's default degree.

Radiative3D currently supports two distinct modeling and meshing approaches: stacked layers of spatially uniform elastic properties, or a warped cartesian grid (WCG) of nodes in which material properties vary linearly based on a tetrahedral tessellation of the space using the node locations as the corners of the tetrahedral model cells.  The former approach is more akin to "1-D" models, and the latter approach fully utilizes the 3-D modeling capabilities of Radiative3D.

Currently, Radiative3D lacks file importers for reading Earth models into memory. To construct Earth models, the model structure must be hard coded using our model-building API.  Changes to the model structure thus require a recompile.  These "user-coded" models are coded in the source file called <tt>user.cpp</tt>, and inspection of that file will reveal several example models to choose from, and explain the model building API.

#### Dimensions and units

For user interaction (input and output), Radiative3D is generally agnostic about choice of physical units.  The only requirements are consistency of choice, and an understanding that the units you use for input will determine the units you get for output.  Specific details are as follows:

##### _Distance and time:_

Radiative3D is designed for Earth-scale models and temporal resolutions of ~1 Hz or greater.  As such there is a general presumption that length quantities will be in kilometers and time quantities will be in seconds, but there is nothing specifically constraining this.  What you use for input is what you will get for output.  E.g., if you use kilometers for your model dimensions, and want temporal output in seconds, then make sure you use kilometers/second as your velocity units.

##### _Frequency:_

For input/output, Radiative3D uses ordinary frequency measurements (i.e., cycles-per-time, not radians-per-time).  Frequency measures will share the time dimension with velocities.  I.e, if you are inputting velocities in length-per-second, then input frequency using cycles/second, or Hertz.

##### _Density:_

For most use cases, Radiative3D sets no restriction on choice of units for densities.  Where densities are used in computation, only ratios come into play.  Thus it is user choice what units to use.

The only exception to this is when Earth models are defined in terms of shear and bulk moduli (a feature not yet supported), instead of P and S velocities.  In this case, the choice of density units must be in agreement with the choice of moduli units to produce the desired length and time units.

##### _Attenuation Q:_

Radiative3D interprets Q values as the inverse of the per-radian fractional rate-of-energy-loss due to intrinsic (e.g. thermal) attenuation. For instance, for a wavefront in a uniform-Q medium, where the energy decays with time _t_ as:

<math>E(t) = E_0\exp\left[-\frac{\omega t}{Q}\right]</math>,

inverse Q is equivalent to:

<math>Q^{-1} \equiv \frac{-1}{E} \; \frac{\mathrm{d}}{\mathrm{d}(\omega t)} E(t)</math>.

Or, thinking of Q as an angular time constant, Q represents the number of radians _(ω*t)_ of wavefront evolution that result in a 1/_e_ energy attenuation of the wave.

Thinking in terms of ordinary cycles, instead of radians, Q is the number of cycles _(f*t)_ of wavefront evolution resulting in an exp[-2π] energy attenuation (about a 99.8% reduction) of the wave.

### EXAMPLES

#### Crust pinch: envelopes

An envelope and travel-time curve simulation in a crustal pinch model can be run via the example script called <tt>do-crustpinch.sh</tt>.  From the Radiative3D install directory, run:

    ./do-crustpinch.sh test01

The script will create an output directory at <tt>./data/YYYYMMDD-HHMMSS-test01-R3D/</tt> and will run Radiative3D inside that directory, collecting all output into various log files and data files.  Upon completion of the simulation, it will process the output and produce a series of figure files depicting envelope traces and travel time curves.  The output directory name is a combination of a time stamp and the run ID provided as the first argument to the script ("test01" in this case).

Note that to get good quality envelopes, tens of millions of phonons will be propagated, and you can expect the simulation to run for a few hours.  For this reason, it may be prudent to run the simulation inside a GNU screen session, especially if you are running it on a remote computer.

#### Crust pinch: video

A video simulation in a crustal pinch model can be run via the example script called <tt>do-crustpinch-vids.sh</tt>.  From the Radiative3D install directory, run:

    ./do-crustpinch-vids.sh test02

The script will create an output directory at <tt>./data/YYYYMMDD-HHMMSS-test02-R3VID/</tt> and will run Radiative3D inside that directory, collecting all output into various log files and data files.  Upon completion of the simulation, it will process the output and produce a series of still frames depicting phonon propagation as a function of time, and will concatenate those frames onto two videos, one depicting a viewpoint from above, and the other depicting a side-elevation viewpoint.

To get good quality images, a few tens of thousands of phonons will be propagated.  You can expect the simulation to run for about ten to twenty minutes, and the video generation to take another twenty minutes after that.  Videos are generally much less computationally intensive to produce than envelopes and travel-time curves are, because for videos every phonon is useful, whereas for envelopes only those phonons that interact with a given seismometer are useful in producing the envelope trace.

### INSTALLATION

#### Via Subversion

Source code for Radiative3D can be obtained via Subversion.  When a version 1.0 release becomes available, the following command will retrieve the source code from our repositories:

    svn co https://rainbow.phys.uconn.edu/svn/Radiative3D/tags/Radiative3D-v1.0.0/ [destdir]

This will download the v1.0.0 release and place the contents in destination directory <tt>destdir</tt>, or else it will create directory <tt>Radiative3D-v1.0.0</tt> if the optional argument <tt>[destdir]</tt> is omitted.

_NOTE: At the present time, we are not ready for a version 1.0 release.  If you want access to the code as it stands right now, please contact the authors, and we will provide you the repository address of the latest development version, along with a long list of "caveat emptor" warnings about using alpha versions of a not-yet-validated codebase.  We do welcome your use of the code, however, as well as any feedback you may wish to provide._

Radiative3D is written in C++ and compilation is automated with GNU Make.  Compiling is simply a matter of typing "<tt>make</tt>" on the command line.

#### Requirements

Radiative3D itself has no requirements other than a reasonably recent C++ compiler.  It is tested to compile and run on the following platforms: Fedora Linux, Ubuntu Linux, and Macintosh OS X.  While we strive for platform-independence in the development of Radiative3D, the accompanying scripts that automate Radiative3D simulation runs and produce visualizations from its output may have cross-platform issues.  Linux is the _recommended_ platform for these scripts. The following packages are required to run the automation and visualization scripts in unmodified form:

* __Bash__ interpreter.

  * This is standard on Linux and OS X systems, but not on Windows.  The scripts are tested under Linux but may use some non-standard commands under OS X.  Cygwin may provide a usable Bash interpreter for Windows.

* __GNU Octave__ — for figure generation

  * GNU Octave is an open-source implementation of the Matlab programming language, and is available for Linux, OS X, and Windows.  Octave is usually available through the package management utilities on modern Linux systems.  On OS X it can be installed via Homebrew, a third-party open source package management system.
  
* __Gnuplot__

  * This is a standard plotting utility on Unix systems and provides the back end to Octave's plotting system. It is available through the system package management tools or via Homebrew on Mac.

* __LaTeX__ (recommended)

  * The visualization scripts use <tt>pdfcrop</tt> to get rid of whitespace in pdf figure output. This tool is typically part of a LaTeX distribution. On Linux, you may need to install package "texlive-extra-utils" to get it. On OS X it should be included in the full installation of MacTex. 

* __epstool__ and __fig2dev__ (optional)

  * These aren't actually needed by the plotting scripts, but Octave will complain about their absence during plotting if they are not installed.  _epstool_ is its own package on the Linux systems we have tested, and _fig2dev_ is part of the _transfig_ package.  Both packages are available through Homebrew on OS X.
  
* __ffmpeg__ (optional; for video generation scripts)

  * ffmpeg is an open-source video encoding tool that we use to concatenate still frames into videos.  It is available for Linux, OS X, and Windows.

### BUGS

#### Limitations of algorithm

The code employs [[radiative transport]] for modeling the propagation of energy within three-dimensional solid bodies.  Radiative transport is an extension of [[ray theory]], which can be thought of as the limiting case of full-wave treatment for infinite frequencies.  Since we use it to model finite-frequency physics, ray theory, (and hence radiative transport), may be inappropriate and give erroneous results for certain Earth-model configurations.  In particular, radiative transport may have trouble with: steep velocity gradients or model features with sizes on the order of a wavelength or smaller, free-surface effects, sharp interface effects, and the propagation of free-surface or interface waves (diffracted phases), or other instances in which waves would propagate along a 2D manifold rather than a 3D space.

The code models phonons at a single frequency per run.  To produce broadband synthetics, multiple simulation runs must be conducted, and the results summed.  Dispersion is not modeled by the code, however it can be explicitly input by using custom Earth models for each desired frequency.

#### Validation status

Radiative3D is in need of validation.  Until it is validated, it should be assumed that results produced with Radiative3D not only may contain, but _do_ contain, errors and inaccuracies in multiple areas, and that extreme caution should be used in drawing any definitive conclusions from the results.  Validation is an important aspect of code development, and it would be foolish to assume that logical and/or implementation errors will not be found — and corrected — during this stage.

Validation would consist of an iterative process of code review and code testing, in which models of varying degrees of complexity are simulated and results compared against expectations.  In some cases, test models will be simple enough that expected results can be computed analytically. In other cases, the test models will be complex, but can be simulated by other methods, and expectations arrived at in that manner.

It should also be noted that validation, once complete, will not be an assertion that the results produced by Radiative3D are physically correct, but only that they can be reasonably expected to portray an accurate picture of the radiative transport method applied to an Earth model.  [[Radiative transport]] itself, and any limitations of the supplied Earth model, both impose constraints on the physical realism that can be expected from the resultant simulation output.

As regards the types of issues that may be revealed during validation, our greatest concern is with the apportionment of energy between the three components of motion for seismic envelope output. (E.g., between radial, transverse, and vertical components.)  There are known (but correctable) issues with the current handling of S-wave polarization for curved ray paths, which will affect channel apportionment.  Additionally, the authors would like to do a close review of the energy apportionment when rays reflect off the surface.  We have a greater degree of confidence in the timing of phase arrivals and in the summed (three-component) energy arriving in any given time window, though these must also be subject to the validation process.

#### Portability

Radiative3D should compile and run under multiple platforms (Linux, OS X, Windows...), however at present the scripts and visualization routines are not platform-independent. The automation scripts are in Bash, and use some commands that are not standard on OS X, which means Linux is the only 100% supported platform.  The graphics are produced using GNU Octave, which, while nominally platform-independent, it has been our empirical observation that graphical results differ from one platform to another, and even from one Linux distribution to another.  The graphics subsystem in Octave seems to be in a continual state of flux.

For this reason, it is our belief that reimplementing our graphics and automation scripts in a more consistently cross-platform language, such as Python, would be a very good idea.  For now, though, this is on our to-do list, but is not a top priority.  

### AUTHORS

Radiative3D is written by [[User:ChristopherSanborn|Christopher J. Sanborn]] and the [[Solid Earth Geophysics Research Group]] at UConn, including significant contributions by [[User:Swalsh|Steven Walsh]].  The code implements algorithms and takes inspiration from two major forerunner projects: [[Raytrace3d]] by William Menke and [[PSPhonon]] by Shearer and Earle.

### COPYRIGHT

The authors have yet to decide on an official copyright to use.  For now, if you have been given a copy of the code by one of the authors, you can assume you have permission to use it and publish results from it.  If you wish to share the code with others, please contact the authors first for permission.  The authors retain ALL RIGHTS to the code and disclaim ALL LIABILITIES for use of the code.
