# VPSC–DRX compression workbench

Double-click **Start_VPSC_DRX.bat** to open the parameter editor. The initial settings reproduce the `auxiliary_final` Mg–Gd–Y–Zn–Zr calibration delivered before the pole-figure ring study. The alloy preset uses 470 °C, positive true compression along specimen Z, rates of 0.05, 0.2 and 1 s⁻¹, and final axial strain 0.95. It includes all 328 original grain means. The original program, fitted parameters and reference simulations are preserved separately in [pre_hollow_calibration](pre_hollow_calibration/README.md).

## Choose the hardening law

In **Recrystallization**, **hardening when off** selects original DD_REV, original Voce, a compatible native material file, or the calibrated K–M law. The default off-choice is the original DD_REV `Mg.sx` requested by the user. Edit its coefficients in **Original material**. Native runs use the material file's own crystal and slip/twin definitions. The DRX-enabled calibration retains its paper K–M law.

[Hardening controls and original-VPSC equivalence](HARDENING_AND_DRX.md) describes the exact comparisons and remaining material limitations.

## Run a simulation

1. In **Test & grains**, choose the rates, temperature, final strain, increment, seed and initial texture. Select **Slip & hardening** to edit a family and its individual slip systems. **Recrystallization** contains the DRX switch and material parameters.
2. In **Texture settings**, specify the strain states, PF planes, specimen directions for IPFs and MTEX kernel halfwidth. For a 0.2 s⁻¹ run through strain 0.95, enter `[0.2]` for the rates and `[0, 0.12, 0.4, 0.6, 0.95]` for the strains. For other rates with different snapshot schedules, save a second settings file.
3. In **Outputs**, tick **Figure (PNG)** and **Source data** independently for each quantity. Then select **Save & run simulation**. A new timestamped results directory holds the settings used, exact solver input decks, native outputs, figures and run status. **Open results** opens the run folder. Figures are in `figures` and numerical exports in `source_data/rate_*`. **Save as** keeps a separate parameter set. **Restore calibrated values** returns the editor to the archived preset.

The completed [example figures](results/calibrated_example/figures) can be opened immediately. Each rate subdirectory contains the individual system activity, grain states and simulation histories.

## Export numerical data

For a finished calculation, select quantities in **Outputs**, choose **Export from a completed run…**, and select the folder containing `settings_used.json`. The workbench reads the strain states and material from that run and writes a fresh folder under `exports`. Tick only the **Source data** column for a data-only export. For new calculations, **Save as…** saves the output selection with the material settings.

| Source-data selection | Content per strain rate |
|---|---|
| Axial stress–strain | `stress_strain.csv` contains positive axial true compressive strain and axial stress in MPa at every saved increment, including the initial solver record. |
| DRX fraction, mean size and total density | `microstructure.csv` uses the stress plot’s axial strain axis. K–M runs supply the DRX fraction (0–1), both mean diameters (µm), mean density (m⁻²) and additional native history fields. Original VPSC reports zero DRX, the material-file hardening diameter when available, and total density for density-based laws. |
| Dislocation density by family | `density_by_family.csv` contains the recorded densities in m⁻². K–M runs also include the family CRSS in MPa. |
| Slip and twin activity | `family_activity.csv` and `system_activity.csv` contain shares of total absolute shear as fractions (multiply by 100 for the plotted percentages). `system_key.csv` maps the column numbers to families, planes and directions for the actual material. |
| PF or IPF | `texture/*.csv` contains the MTEX 1° grid, with projected x/y coordinates, unit-vector components and density in MRD. Directions are in specimen coordinates for PFs and crystal coordinates for IPFs. `texture_index.csv` lists the indices, target and actual strains, and the refined peak MRD for each map. |
| Grain states and orientations | `grains/MG_SNAPSHOTS.csv` maps requested to actual strains and filenames. Each grain file carries native Bunge Euler angles in degrees and current volume fractions. K–M files also carry grain size, density and generation. |
| DRX cycles | Native generation-volume and cycle-kinetics CSVs. The `strain` column in `DRX_CYCLES.csv` is equivalent strain. |
| All native files | `native/` contains byte-for-byte copies of the rate folder’s `.OUT` and solver/reporting CSV files, including the native equivalent-stress history where available. |

The export’s `output_manifest.json` lists files, selections and any unavailable data. The raw calculation records stay in the original `rate_*` folders. The bundled preset selects every figure and the stress, microstructure, density, activity and grain-state tables. PF/IPF grids, cycle data and the full native copies can be selected independently.

## Interpret the figures

| Output | Definition |
|---|---|
| `stress_strain.png` | Positive true axial compression stress −σzz against positive true axial compressive strain. |
| `microstructure.png` | Historical DRX volume fraction, volume-weighted mean diameter, multiplicity-weighted number mean diameter, and volume-mean total dislocation density. The historical fraction includes recrystallized grains that suffered further deformation. |
| `dislocation_modes_rate_*.png` | Volume-mean density for each slip family and their sum, in m⁻². |
| `slip_activity_rate_*.png` | Contributions of each family and each indexed system. At a given increment, a system contributes Σᵢ wᵢ abs(γ̇ᵢₛ) divided by the corresponding sum over all systems. Each panel for an individual family retains the aggregate denominator. |
| `texture_rate_*_state_*.png` | MTEX PF of each requested plane family and IPF of each requested specimen direction. The default is the basal (0001) PF and compression-axis Z IPF. Colors represent multiples of a random distribution, on a shared scale across the run. |
| `figures/texture_metrics.json` | MTEX version, peak MRD, current number of grain records and requested versus actual axial strain. An off-grid request uses the first completed increment reaching the target within solver tolerance. It does not interpolate orientations. |

At each saved strain, MTEX Python 1.0.0 reads the native Bunge orientations and **current volume fractions**. PFs include symmetry-equivalent plane poles on the upper hemisphere, with Z at the center, X right, Y up and equal-area projection. IPFs map a specimen direction into crystal coordinates, use crystal symmetry, and display the fundamental sector in stereographic projection. Both use an antipodal de la Vallée Poussin kernel. The preset has 5° halfwidth at half maximum. This was the earlier project's approximate comparison to AZtec's 10° setting on this dataset. No intensity rescaling is applied. [MTEX Python project](https://github.com/mtex-toolbox/pymtex)

## Editable physical inputs

| Section | Fields and units |
|---|---|
| Crystal | Symmetry, lattice lengths in Å and angles in degrees. Available settings are cubic, hexagonal, trigonal with hexagonal axes, tetragonal and orthorhombic. Purely geometric native axes are normalized to the length of a. Other controls include reference Burgers vector (m), forest coefficient χ, reference shear modulus (Pa), its temperature slope (K⁻¹), reference temperature (K), auxiliary Poisson ratio, slip exponent and recovery reference rate (s⁻¹). |
| Slip & hardening | One density per grain per family. Each record has friction resistance (MPa), Burgers ratio, χ ratio, storage k₁ (m⁻¹), recovery k₂ and its rate exponent, and initial/newborn densities (m⁻²). Systems within a family share its resistance. Use the actual plane and direction indices for each bidirectional variant. Hexagonal rows contain eight integers, e.g. `0 0 0 1   2 -1 -1 0`. Cubic rows contain six. The workbench checks the zone law and rejects duplicate systems. |
| DRX | Switch, critical total density (m⁻²), nucleation coefficient (m⁻² sʳ where r = q−1), rate exponent q, mobility (m⁴ J⁻¹ s⁻¹), embryo diameter (m), line-energy factor, and random or Gaussian birth orientation with angular spread (degrees). With DRX off, the selected hardening law controls the calculation. The paper_km choice suspends nucleation and radius evolution while retaining K–M hardening. |
| Test & grains | Isothermal temperature (°C), positive rates (s⁻¹), final true axial strain, maximum increment, RNG seed, grain CSV and orientation convention. The native law uses one phase, one element, fixed spherical shape and Z compression. |
| Texture settings | Up to 16 sorted strain states per run, PF planes, specimen XYZ vectors for IPFs, and kernel halfwidth (degrees). |

**Grain representation.** The archived choice, `measured_cohorts`, uses individual equivalent-circle diameters and EBSD grain areas as volume proxies. A numerical grain record can represent a fractional number of equal-sized spheres. The separately selectable equal-sphere modes assign the same radius, shape and weight 1/N to every orientation. The suffix selects the arithmetic or area-weighted mean used to choose that common diameter. These options change the initial texture weights and the simulated microstructure.

**Euler convention.** `archived` passes the exported Bunge angles directly into VPSC and uses the original interpretation in the plots. This is the setting that reproduces the pre-hollow results. For new Oxford/AZtec Mg maps, choose `aztec` to reconcile the EBSD crystal x ∥ a* frame with VPSC x ∥ a. The converter subtracts 30° from φ₂ on input and adds it back for MTEX. `vpsc` accepts native angles and applies only the plotting conversion. All three choices keep specimen Z as the compression axis. Native state CSVs always retain native angles.

**Scope.** Only the DRX-enabled Mg preset has been calibrated here. Changing crystal structure, slip families, temperature, or parameters produces a new model input. The three fitted Mg rates are the experimental basis at 470 °C. In particular, changing temperature updates the shear modulus and process temperature. The supplied mobility and nucleation coefficients have no calibrated cross-temperature law. The calibrated DRX model accepts eight families and 48 systems, all slip. Native VPSC uses the selected material file and can include its twinning modes. Its original limit is 16,000 grains. The 328-grain CSV contains grain averages, not all indexed EBSD pixels.

## Command line and source

```bat
Run_VPSC_DRX.bat --validate
Run_VPSC_DRX.bat --settings settings.json
Run_VPSC_DRX.bat --plots-only results\calibrated_example
Run_VPSC_DRX.bat --export-from results\calibrated_example --figures --data stress microstructure activity
```

The Windows release includes the calibrated executable, the separately compiled original VPSC, and the calibrated five Fortran sources under `source`, with the mode-count and diagnostic updates in `vpsc8drx.f90`. To rebuild, run `source/build.bat` with MinGW gfortran available on PATH or set `VPSC_GFORTRAN` to the compiler executable. The original, unmodified sources and executable are under `pre_hollow_calibration/original_native`. [Model and crystal audit](MODEL_AUDIT.md)

The launcher uses a local virtual environment if present, then this machine's Anaconda interpreter, then Python on PATH. `VPSC_PYTHON` overrides the interpreter. Windows x64 Python 3.12+ with NumPy 2+, SciPy, Matplotlib and Tk is required. The Windows MTEX Python and FINUFFT packages are in `dependencies`. On another machine, install the Python packages in `requirements.txt` if needed. All project data for the default run is in this directory.

A completed run records `status = complete` in `run_status.json`. A failed run retains its inputs and logs. A native abort at Eq. 10 probability ≥ 1 calls for a smaller maximum increment. Use `--plots-only` to rebuild the same figures from completed native output. The script `verify_release.py` checks the editable controls and compares the full calibrated example against its frozen archive.

