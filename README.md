MATLAB source code and example simulation outputs for the manuscript:

**Optically programmed acoustic field synthesis at an electrically passive aperture**

The code simulates ultrasound generation and beamforming by a hexagonally packed fiber-bundle aperture using the k-Wave MATLAB toolbox. It supports on-axis focusing, beam steering, two-beam synthesis with concentric subapertures, and helical-wavefront generation. Acoustic excitation is represented by prescribed single-cycle tone bursts; the code does not simulate optical absorption or thermoelastic conversion explicitly.

Repository: <https://github.com/HuangYanwei-pixian/Supplementary_OPAS>

The default run is a self-contained simulation demonstration, whose main script is `simulation_opas_beamforming.m`. 
Experimental source data will be made public after publication.

## 1. System requirements

### Software and tested versions

| Component | Tested version or requirement |
|---|---|
| Operating system | Windows 11 Pro, version 25H2 |
| MATLAB | R2024a; also checked with R2024a Update 7 |
| k-Wave MATLAB toolbox | Version 1.4 |
| Parallel Computing Toolbox | Version 24.1 (R2024a), for the GPU execution path |
| GPU driver | NVIDIA driver 591.86 on the reference workstation; the script enables CUDA forward compatibility for the newer GPU |

Other MATLAB releases and operating systems have not been validated for this repository. MATLAB requires a separate license. k-Wave is installed separately; see the [k-Wave website](https://www.k-wave.org/) and [documentation](https://www.k-wave.org/documentation.php).

The script uses MATLAB implementations of `kspaceFirstOrder3D` and `kWaveArray`. No separate k-Wave C++/CUDA executable or code compilation is required.

**CPU execution:** select `model = 1`. GPU execution requires Parallel Computing Toolbox and a supported or appropriately configured NVIDIA GPU.

### Reference workstation

| Component | Configuration |
|---|---|
| CPU | 12th Gen Intel Core i7-12700K, 3.60 GHz |
| RAM | 32 GB |
| GPU | NVIDIA GeForce RTX 5060 |
| Operating system | Windows 11 Pro, version 25H2 |

A GPU is optional when the CPU path is used. The workstation above is a tested configuration, not a measured minimum requirement. Memory requirements increase with grid size, source frequency, spatial sampling density, and recorded time-series size.

## 2. Installation guide

1. Install and activate MATLAB R2024a. Install Parallel Computing Toolbox if using GPU mode.
2. Install k-Wave 1.4, either through MATLAB Add-On Explorer or by downloading and extracting the toolbox. If installed manually, add the directory containing `kWaveGrid.m` and `kspaceFirstOrder3D.m` to the MATLAB path:

   ```matlab
   addpath('C:\path\to\k-Wave');  % replace with the actual toolbox directory
   savepath;                     % optional: persist the path between sessions
   ```

3. Download this repository using **Code > Download ZIP**, or clone it with Git. On Windows, use a short destination path, such as `C:\work\Supplementary_OPAS`, because some example output filenames are long.
4. Open MATLAB and change its Current Folder to the repository's `simulation_code` directory:

   ```matlab
   cd('C:\work\Supplementary_OPAS\simulation_code');
   ```

5. Confirm that MATLAB can find the required files:

   ```matlab
   which kWaveGrid
   which kWaveArray
   which kspaceFirstOrder3D
   which simulation_opas_beamforming
   ```

   Each command should return a file path. The main script adds `array_functions` automatically.

**Typical setup time:** approximately 5–10 minutes on a desktop with MATLAB already installed. This is an estimate, not a timed installation benchmark; it covers obtaining k-Wave and this repository and configuring paths. MATLAB/Parallel Computing Toolbox installation, license activation, and driver installation are excluded and depend on the installation and network environment.

## 3. Demo

### Demo inputs

The demo creates its synthetic input data at runtime: fiber-element coordinates, channel-specific emission delays, tone-burst source signals, a homogeneous acoustic medium, and the computational grid. No external experimental dataset is required. Example simulated output data are included as MATLAB `.fig` files in [`simulation_code/sim_results`](simulation_code/sim_results/); these are saved outputs, not input datasets read by the solver.

The default configuration is:

| Setting | Value |
|---|---|
| Aperture | 61 elements, 125 µm center-to-center spacing |
| Active disc diameter | 105 µm |
| Source waveform | One-cycle tone burst at 13 MHz |
| Source amplitude parameter | 1 MPa |
| Programmed axial focal depth | 1 mm |
| Azimuth/elevation steering | 0° / 0° |
| Medium | Sound speed 1,480 m/s; density 1,000 kg/m³ |
| Nominal computational region | 2 × 2 × 4 mm³ |
| Spatial sampling | `ppw = 5` |
| Simulated time | `t_end = 3e-6` s |
| CFL number | `cfl = 0.2` |
| Execution | `model = 2` (GPU) |

### Run the demo

From the `simulation_code` directory, run:

```matlab
simulation_opas_beamforming
```

The script clears the workspace and closes existing figures at startup. Save any unrelated work first. To select CPU mode, edit `model` near the top of the script as described in Section 1.

### Expected output

MATLAB displays the aperture layout, a channel-delay map, the normalized depth–lateral pressure field, and axial and lateral pressure profiles. The default field should show an on-axis focal region near the programmed 1 mm depth. The actual peak position and beam width depend on the finite aperture, excitation, and numerical sampling.

Six `.fig` files are saved for the default run, with prefixes:

- `Time delays`
- `Ultrasound Field (Azi_norm)`
- `Profile (Axial)` and `Profile (Axial_norm)`
- `Profile (Azi)` and `Profile (Azi_norm)`

The files are written to `simulation_code/sim_results`. Re-running the same parameter combination overwrites files with the same names. To inspect an output, open it in MATLAB, for example:

```matlab
openfig(fullfile('sim_results', ...
    'Ultrasound Field (Azi_norm)_enum61_freq13MHz_f1000um_def0(azi)_0(ele)deg.fig'));
```

Pressure fields are based on `sensor.record = {'p_max'}`: they represent **peak positive pressure**, not peak-to-peak pressure. Normalized plots are divided by their respective maximum values. The specified source amplitude is a model input and is not a calibration of experimentally measured pressure.

### Expected run time

The reported reference calculation time on the workstation in Section 1 is **20.067 s**, as printed by k-Wave under `total computation time` for the reference demo. This is the k-Wave calculation time, not an end-to-end timing of the entire script. MATLAB startup, array/source construction, GPU initialization or first-use compatibility compilation, plotting, and figure saving add overhead. Allow extra time for the first GPU run. CPU execution and larger simulations may take substantially longer.

## 4. Parameter settings and instructions for use

Edit the assignments in the script, save it, and run it again. Setting variables in the Command Window before running the script does not override its parameters, because the script begins with `clear` and assigns its own values.

All distances are in **metres**, frequencies in **hertz**, pressures in **pascals**, and angles in **degrees**, except `opas.time_step`, which is expressed in **nanoseconds**.

### 4.1 Aperture and acoustic-source parameters

| Parameter | Default | Meaning and guidance |
|---|---:|---|
| `opas.d_num` | `9` | Number of fiber centers along the principal diameter; use a positive odd integer. `5` gives 19 elements and `9` gives 61. The total is `(3*d_num^2 + 1)/4`. |
| `opas.d_core` | `125e-6` | Fiber **outer diameter including cladding**, also used as the element spacing. Despite the variable name, this is not the 105 µm optical core diameter. |
| `opas.element_diameter` | Calculated | Set in `array_fiber_bundle_coordinates.m` to `0.84*opas.d_core`, giving 105 µm for the default spacing. Change this relation separately if a different active-diameter/spacing ratio is required. |
| `opas.source_f0` | `13e6` | Tone-burst frequency for single-beam and helical modes. This is an excitation-model parameter, not a measured local spectral centroid. |
| `opas.source_amp` | `1e6` | Pressure-amplitude scaling of the prescribed source waveform. |
| `opas.source_cycles` | `1` | Number of cycles in the tone burst. Larger values give a longer, generally narrower-band waveform and may require a longer simulation. |
| `c0` | `1480` | Homogeneous-medium sound speed. The script also assigns it to `opas.speed` and `medium.sound_speed`; keep these consistent. |
| `medium.density` | `1000` | Medium density in kg/m³. |

The supplied medium is homogeneous and lossless; absorption, heterogeneous tissue properties, and nonlinear material parameters are not configured by default.

### 4.2 Focusing, steering, and structured fields

| Parameter | Default | Meaning and guidance |
|---|---:|---|
| `opas.source_focus` | `1e-3` | Programmed axial focal depth. A value of `0` disables the spherical focusing contribution. Keep the target inside the simulated region. |
| `opas.source_deflect_azi` | `0` | Azimuthal steering angle in the x–z plane. For example, use `10` for +10° steering. |
| `opas.source_deflect_ele` | `0` | Elevational steering angle in the y–z plane. The ordinary output samples one fixed x–z plane, so an elevation-steered peak may leave that plane. |
| `opas.multi_beam_enable` | `0` | Set to `1` to use concentric inner and outer subapertures. Use the 61-element configuration for the supplied 19/42-element example. |
| `opas.helical_wavefront_enable` | `0` | Set to `1` for the helical-wavefront example. Do not enable it simultaneously with two-beam mode. |
| `opas.tpl_charge` | `-3` | Signed topological-charge parameter used to construct the helical delay pattern. |
| `show_lateral_depth` | `1e-3` | Depth used to extract a lateral profile or an x–y field. Update this independently when changing the focal depth. It is assigned later in the script. |
| `opas.time_delay_digitized` | `1` | Enables quantization of the displayed/exported hardware delay values; see the implementation note below about the delays currently used by the solver. |
| `opas.time_step` | `2` | Hardware-delay quantization step in ns. This is distinct from the k-Wave numerical time step. |

In the code's coordinate convention, a positive azimuth angle directs the beam toward positive x. Steering is constructed by adding a linear delay gradient to the spherical focusing delay; geometric focus predictions are approximate for a finite aperture and large steering angles.

For **two-beam mode**, set `opas.multi_beam_enable = 1` and keep `opas.helical_wavefront_enable = 0`. The supplied parameters are:

| Parameter | Inner subaperture | Outer subaperture |
|---|---:|---:|
| Spiral element IDs | 1–19 | 20–61 |
| Frequency | `source_f0_inner = 15e6` | `source_f0_outer = 7.5e6` |
| Axial focal depth | `source_focus_inner = 0.75e-3` | `source_focus_outer = 1.5e-3` |
| Azimuth angle | `source_deflect_azi_inner = -10` | `source_deflect_azi_outer = 10` |
| Elevation angle | `source_deflect_ele_inner = 0` | `source_deflect_ele_outer = 0` |

All parameter names in this table are fields of `opas`. The group-specific parameters replace the single-beam focus, steering, and frequency settings in this mode. The 15/7.5 MHz model frequencies are illustrative simulation settings; they are not the experimental nominal values or the local spectra of the combined field. The single `show_lateral_depth` setting extracts only one lateral slice at a time.

For a **helical-wavefront example**, set:

```matlab
opas.multi_beam_enable = 0;
opas.helical_wavefront_enable = 1;
opas.tpl_charge = -3;
```

This mode produces an x–y pressure map and a peak-arrival-time map at `show_lateral_depth`. Arrival times are computed from the largest absolute pressure sample; they are not complex spectral phase measurements.

### 4.3 Numerical simulation parameters

| Parameter | Default | Meaning and guidance |
|---|---:|---|
| `model` | `2` | `1`: MATLAB CPU solver using single precision. `2`: MATLAB GPU solver using `gpuArray-single`. See the CPU compatibility instruction in Section 1. |
| `grid_size_x`, `grid_size_y`, `grid_size_z` | `2e-3`, `2e-3`, `4e-3` | Nominal physical region. Include the aperture, all intended foci, and the field region of interest. Actual dimensions are adjusted to even grid-point counts. |
| `ppw` | `5` | Spatial points per wavelength at the frequency used to construct the grid. Higher values give finer sampling but increase memory and runtime rapidly. |
| `cfl` | `0.2` | Courant–Friedrichs–Lewy number used by `kgrid.makeTime` to choose the time step. Lower values increase temporal sampling and runtime. |
| `t_end` | `3e-6` | **Simulated physical duration**, not wall-clock computation time. Increase it for deeper targets, longer source pulses, or larger emission delays. |

The grid spacing is `dx = c0/(ppw*f_grid)`, where `f_grid` is `source_f0` in ordinary/helical mode and `max(source_f0_inner, source_f0_outer)` in two-beam mode. The `ppw` setting is a sampling parameter, not a guaranteed accuracy percentage. A single-cycle pulse is broadband: quantitative work should check whether the relevant high-frequency content is resolved and whether results converge when spatial and temporal sampling are refined.

Increasing `ppw` by a factor of two produces approximately eight times as many 3D grid points and about twice as many time steps at fixed CFL and duration, before FFT-size and implementation effects. Changing source frequency also changes the grid and computational cost. For an initial demonstration, retain `ppw = 5` and `cfl = 0.2`.

Choose `t_end` to include the maximum emission delay, propagation time to the most distant recorded location, the source-pulse duration, and a margin. Increasing `t_end` does not enlarge the spatial region. The solver uses an external perfectly matched layer (`PMLInside = false`) with automatically selected thickness.

For unattended runs, adding `'PlotSim', false` to `input_args` disables k-Wave's live animation while retaining the final visualization code.

### 4.4 Using a different configuration or your own inputs

For a different regular fiber bundle, change the aperture, source, medium, and beamforming settings above. No measured input file needs to be imported for these forward simulations.

An irregular aperture or measured source waveform requires code changes: replace the coordinate-generation step and/or the `toneBurst` source construction. Preserve one source-signal row per physical array element in the same order used to add elements to `kWaveArray`, and sample waveforms at `1/kgrid.dt`. The existing script does not provide an experimental-data import or analysis interface.

## 5. Code organization and computational workflow

```text
Supplementary_OPAS/
├── README.md
├── LICENSE
└── simulation_code/
    ├── simulation_opas_beamforming.m
    ├── array_functions/
    │   ├── array_fiber_bundle_coordinates.m
    │   ├── array_focusing_deflection.m
    │   ├── array_multi_beam.m
    │   ├── array_helical_wf.m
    │   └── array_time_digitize.m
    └── sim_results/
        └── *.fig
```

The main computation is:

1. Generate the hexagonal element coordinates and spiral element IDs.
2. Calculate spherical focusing and linear steering delays; optionally partition the aperture or add a helical delay profile.
3. Calculate the displayed/exported quantized hardware delays.
4. Build the k-Wave spatial/time grid and delayed tone-burst signals for each source element.
5. Map the finite disc elements to the numerical grid using `kWaveArray`.
6. Propagate pressure using `kspaceFirstOrder3D` on the selected CPU/GPU path.
7. Extract pressure fields/profiles and, in helical mode, peak-arrival times; display and save figures.

6. Current implementation notes and manuscript reproduction

The following behaviors apply to the documented source revision:

- The solver's source construction reads `time_delays_show`, which contains delays before hardware quantization. Changing `time_step` changes the displayed/exported hardware delays but does not currently change those solver inputs. Do not interpret these runs as a validation of 2 ns hardware quantization. The numerical signal offsets are also rounded to the k-Wave time grid.
- `addDiscElement` receives `[0,0,0]` as its axis point, so the finite discs point toward the computational-grid origin rather than all having parallel normals. This affects the simulated source geometry.
- The ordinary sensor uses index `Ny/2`; for the even-sized k-Wave grid this is one grid spacing away from y=0. Its profiles therefore represent that sampled plane.

The default demonstration uses 13 MHz, whereas Supplementary Fig. S2 uses 27 MHz. Match the source frequency, aperture geometry, and beamforming settings before making quantitative comparisons.

For the same physical beamforming configuration, numerically extracted quantities such as lateral and axial full width at half maximum (FWHM) can vary with the spatial and temporal resolution of the simulation, including `ppw` and the time step controlled by `cfl`, as well as the sampling and interpolation used to extract the profiles. These quantities should be assessed for convergence as the numerical resolution is refined. Supplementary Note 1 includes both analytical estimates (Table S1) and numerical simulations (Fig. S2); the analytical estimates additionally rely on idealized aperture and excitation assumptions. Numerical resolution can contribute to quantitative differences, but its contribution to a particular discrepancy must be established by a convergence check.

An exact per-panel reproduction configuration has not yet been packaged in this repository. The two-beam example illustrates the programming method and is not an experimental-data reconstruction of Fig. 6.

The mathematical formulation and numerical-simulation description are provided in the manuscript Methods subsections **“Fiber-bundle array design and geometry”** and **“Numerical simulation of the acoustic field,”** and Supplementary Note 1, particularly **Section 1.1, “Beamforming geometry and delay synthesis,”** and **Section 1.4, “Numerical simulation of the ultrasound phased array.”**

## 7. License and contact

The repository source code is distributed under the **Apache License, Version 2.0**; see [LICENSE](LICENSE). MATLAB and k-Wave are separate dependencies subject to their own licensing terms.

Contact: Yanwei Huang — <h.yanwei@wustl.edu>.

Last update: 06/October/2026
