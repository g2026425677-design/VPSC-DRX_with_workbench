# Model and crystal audit

The original-VPSC off-mode and the user-selectable hardening controls were added on 2026-11-05. They are documented and verified separately in [Hardening and DRX](HARDENING_AND_DRX.md). The equations below describe the calibrated K–M branch.

## Which calibration was restored

The original delivery identifies `auxiliary_final` as the selected calibration. Its material dictionary matches the baseline saved before the hollow-PF investigation. All 328 measured grain orientations, original diameters, area fractions, native compression conditions, calibration seed and fitted material coefficients were restored. The default uses the historical orientation convention. The editor exposes the subsequently verified AZtec frame conversion separately.

The reference model omits the previously questioned original-grain resistance offset, the global CRSS rate multiplier, and the shear-age-to-GOS conversion. The three initial family densities can describe prior stored dislocation content. No stress fluctuations or target pole-figure shape are imposed. See the archived [equation-by-equation comparison](pre_hollow_calibration/FORMULA_AUDIT.md) and the [original calibration results](pre_hollow_calibration/calibration_results.md).

## Equations actually evaluated

For grain i and slip family α, the native VPSC power law uses the macroscopic equivalent strain rate as the reference rate. The preset uses exponent n = 20. The family shear increment is Δγᵢα = Δt Σₛ∈α |γ̇ᵢₛ|. Evolution uses dρᵢα/dγᵢα = k₁α √ρᵢα − k₂α ρᵢα. The family resistance is τf,α + χα μ bα √ρᵢα. The total grain density entering DRX is ρᵢ = Σα ρᵢα. Basal, prismatic and pyramidal-I slip each have their own density in each grain. Individual variants in a family share that family's resistance. [Zhou et al., 2017](https://doi.org/10.1016/j.ijplas.2017.01.001)

Nucleation requires ρᵢ ≥ ρc. Each eligible record has probability pᵢ = C ε̇ᵠ 4πrᵢ² Δt. A nucleus has the specified radius, low initial densities and the next generation number. The growth law is ṙᵢ = M λ (Σⱼ wⱼ ρⱼ − ρᵢ), with λ = f μ b² and f = 0.5 in the preset. The native code updates radii, removes consumed grains and renormalizes volume fractions. The Gaussian birth rule samples a surviving grain on which nucleation has occurred and rotates its orientation about a random axis. The alternative random rule samples a uniform orientation. [Zhou et al., 2017](https://doi.org/10.1016/j.ijplas.2017.01.001)

The workbench reads the configured family definitions instead of assuming 3+3+12 slots. It uses the same equations on those active slots. It exports each system's volume-weighted share of absolute slip rate at the native recording point. Neither the export nor MTEX plotting feeds any quantity back into the solver.

## Material transfers and remaining assumptions

- The source paper demonstrated the DRX model for FCC copper. The preset transfers the model to an effective single-phase HCP Mg matrix. The calibrated three-family fit does not independently identify each family's density from macroscopic stress, EBSD GND, or GOS.
- At the single fitted temperature, k₂α = k₂α,ref (ε̇axial / 0.2 s⁻¹)⁻ᵐᵅ. This substitutes prescribed axial rate for the shear-rate argument in the paper's recovery parameterization. The KM update still uses the actual shear increment of its own family. C and M are fitted effective coefficients at 470 °C.
- Growth uses one scalar a-type Burgers vector for the line energy. Slip resistance uses the appropriate family Burgers magnitude. The source copper model does not specify a unique HCP multimode line-energy mixture.
- EBSD supplies a 2D equivalent-circle diameter and area. The archived model interprets these as sphere diameter and volume-fraction proxy. Its fixed cohort multiplicity satisfies wᵢ ∝ mᵢrᵢ³. A newborn inherits its parent's multiplicity. This is a numerical representation choice. The initial orientation of each record is the grain mean, so it loses the intragranular spread of a pixel map. The larger final EBSD maps supplied orientation references only.
- DRX fraction is the volume of all grains of generation 1 and higher. They can accumulate dislocations again. The EBSD GOS < 2° area fraction is a reference and can be smaller. Similarly, the modeled total dislocation density and EBSD GND measure different quantities.
- The original grain exporter passed its Euler numbers directly to native VPSC. The archived preset reproduces those reference outputs. For matching the Oxford crystal basis of new data, use `aztec` in the editor. This changes the coordinate conversion. No parameter was fitted to the hollow.

## Magnesium crystal constants

| Input | What was checked |
|---|---|
| Crystal and lattice | The original AZtec `.h5oina` Mg phase specifies space group 194, a = b = 3.2089 Å, c = 5.2101 Å and angles 90°, 90°, 120°. A PNNL Mg-matrix EBSD study independently reports the same reference cell. The resulting c/a = 1.62364. The solver's rounded 1.624 differs by 0.022%. The editable preset retains the exact archived ratio by setting c = 1.624a. These are phase indexing values rather than an in-situ measurement of the user's alloy at 470 °C. [PNNL report, methods, p. 6](https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-32134.pdf) |
| Slip systems | The archived deck has 3 basal {0001}〈11−20〉, 3 prismatic {10−10}〈11−20〉, and 12 first-order pyramidal {10−11}〈11−23〉 variants. The editor checks all four-index constraints and the plane/direction zone law for every variant. Whether pyramidal-I alone accounts for nonbasal c+a slip in this particular alloy is a material assumption. [Pyramidal I slip in c-axis compressed Mg single crystals](https://www.sciencedirect.com/science/article/abs/pii/S1359646215004066) |
| Burgers magnitude | The preset uses a-type b = 0.320 nm and c+a magnitude b√(1+(c/a)²) = 0.6103 nm. The reference Mg cell gives 0.32089 and 0.61190 nm respectively. [Dobeš and Dymáček](https://doi.org/10.1016/j.jma.2020.03.001) |
| Shear modulus | The scalar pure-Mg correlation G(T) = 16.6[1−0.00053(T−300)] GPa gives 12.7012 GPa at 743.15 K. It supplies the forest term and line energy. Its transfer to this RE alloy is approximate. [Dobeš and Dymáček](https://www.sciengine.com/doi/pdfView/52B7713D9AA74334BC29B71E8B7CBD6E) |
| Elastic tensor | The preset builds the auxiliary native stiffness from G and assumed ν = 0.29, giving C₁₁ = 42.9420, C₁₂ = 17.5397 and C₄₄ = 12.7012 GPa. Real Mg has five independent hexagonal elastic constants. The auxiliary isotropic matrix supplies VPSC's required elastic input. There is no measured 470 °C single-crystal tensor or per-mode CRSS for the exact alloy in this package. [Slutsky and Garland, ultrasound measurements](https://doi.org/10.1103/PhysRev.107.972) |

## Verification

The [archive comparison](verification/archived_comparison.json) compares all numeric values in the stress, DRX, generation and mode histories and all five full grain snapshots at each rate. The [editable-input checks](verification/editable_controls.json) exercise different strengths, the DRX switch, off-grid snapshot requests, an equal-sphere input and a cubic FCC slip definition. The cubic case checks the software path only. The archive reproduces its own numerical history. The experimental stress RMSEs for the archived Mg fit remain 2.229, 5.311 and 11.033 MPa at 0.05, 0.2 and 1 s⁻¹, respectively.

