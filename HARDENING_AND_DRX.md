# Hardening-law selection and original VPSC comparison

Updated 2026-11-05.

## What the DRX switch now means

| DRX | Hardening when off | Actual calculation |
|---|---|---|
| On | Any saved off-choice | Calibrated mode-resolved K–M resistance and the paper's nucleation/growth scheme. |
| Off | `original_dd_rev` | Original VPSC, with its own DD_REV material file, IHARDLAW = 23. This is the default off-choice selected by the user. |
| Off | `original_voce` | Original VPSC, with the selected Voce material file, IHARDLAW = 0. |
| Off | `original_material` | Original VPSC, using the hardening law and coefficients read from the selected native .sx. The workbench supports single-phase native law 0, 1, 20 or 23. |
| Off | `paper_km` | Calibrated K–M slip hardening, with nucleation and grain-boundary growth disabled. This isolates the DRX population mechanism. |

In the editor's **Recrystallization** tab, set **enabled** and **hardening when off**. Choosing DD_REV or Voce selects that original Mg preset. Use **Original material** to inspect and edit its coefficients. **Browse .sx** selects another compatible native material. **Save & run simulation** saves edited material text into a separate user-material file and passes those exact bytes to VPSC. The material's hardening-law number must have the matching parameter block below it.

For native VPSC runs, the **Crystal** and **Slip & hardening** tabs are disabled. The .sx file supplies the crystal, elastic constants, active slip/twin systems, rate sensitivity and all hardening coefficients. The **Test & grains** and **Figures** controls continue to set the texture, Z-compression test, output strains and MTEX smoothing. A frozen copy of the selected material is saved as `original_material.sx` within each result folder before any rate starts.

## Why K–M also changes stress

The paper's Kocks–Mecking equation evolves each slip family's density,

\[
\frac{d\rho_\alpha}{d\gamma_\alpha}
 = k_{1\alpha}\sqrt{\rho_\alpha}-k_{2\alpha}\rho_\alpha .
\]

The calibrated slip resistance uses that state,

\[
\tau_\alpha =
\tau_{0\alpha}+\chi_\alpha\mu b_\alpha\sqrt{\rho_\alpha}.
\]

Consequently, a density change affects slip activation, stress and lattice rotation, even with DRX disabled. The grain's summed density also enters nucleation and growth. These are two uses of the density state. The earlier population-only switch left this K–M resistance active.

The selected original `Mg.sx` uses VPSC's reversible dislocation-density law 23, including forest and reversible density components and its own coupling parameters. Matching that original calculation requires running that law with its material input. Simply disabling K–M-based nucleation does not establish this equivalence.

## Exact original model used

The reference source was copied byte for byte from the user's folder

`E:\直博一年级\宝钛Mg合金仿真项目\VPSC-DRX_reproduction\VPSC_code-main (1)\VPSC_code-main (1)\VPSC_code-main`

The unmodified source is under `source/stock_vpsc`. The default material is an unchanged copy of its `ex13_dd_rev/Mg/Mg.sx`, now under `data/original_vpsc/Mg.sx`. The Voce alternative is an unchanged copy of its `ex06_torsion/Mg_voce.sx`.

A normal native-mode run first executes `bin/vpsc_original.exe`, compiled from those three unmodified sources with the original double-precision and legacy flags. Its standard files are the authoritative outputs. A companion build then runs the same inputs with a read-only output observer linked around the original activity-writing call. The observer exports full-precision orientations, individual slip/twin shares and density for plotting. It does not assign any VPSC state. Each run verifies byte-for-byte equality of the companion's `STR_STR.OUT`, `ACT_PH1.OUT`, `TEX_PH1.OUT` and, when present, `DIS_DEN.OUT` before accepting its additional diagnostics.

This two-run verification adds runtime. A discrepancy fails the result and retains the original native files.

## Verification results

The direct reference used all 328 measured orientations and their area weights, fixed spherical shape, Z compression, interaction = 3, effective exponent = 5, increment 0.0005, 470 °C and final axial strain 0.95. The original material file was supplied unchanged.

| Rate, s⁻¹ | Final original axial stress, MPa | Full history records | Direct-reference equality |
|---|---|---|---|
| 0.05 | 466.678694 | 1901 | Exact for stress, activity, texture and density files |
| 0.2 | 466.678901 | 1901 | Exact for stress, activity, texture and density files |
| 1 | 466.658813 | 1901 | Exact for stress, activity, texture and density files |

The complete file hashes and changed-parameter checks are in [verification.json](verification/original_vpsc/verification.json). The GUI end-to-end test is in [editor_test.json](verification/original_vpsc/editor_test.json). The calibrated DRX-on executable is unchanged, and regenerated physical input files at all three rates match the prior calibrated example byte for byte.

The original Mg.sx contains example material coefficients, including several labeled for 20 °C and a fixed 20 µm hardening grain size. Its outputs above demonstrate software equivalence. Those coefficients have not been recalibrated to the user's rare-earth alloy at 470 °C. The **Original material** editor makes them user-defined.

For native DD_REV, the density figure reads the original forest-plus-reversible density state and is checked against its native `DIS_DEN.OUT`. The fixed 20 µm value is labeled **hardening grain size from .sx**. Original VPSC does not evolve a DRX grain population in this mode. The native Voce law has no dislocation-density state, so its figure reports that limitation explicitly. Activity figures include the selected native twinning modes.

## Reproduce the checks

Run `python verify_stock.py` from the workbench folder. It verifies the recorded full direct-reference comparison, runs edited DD_REV, Voce and K–M-off cases, checks material changes reach the native input and confirms the calibrated DRX-on input decks remain identical. Use `source/build_original.bat` to rebuild the original and read-only observer executables.

