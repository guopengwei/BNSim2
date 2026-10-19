# Parameter records and dimensional corrections

`parameters.csv` stores every configured agent/model value with a unit, status
and provenance note. `metadata.json` stores species/material coefficients,
configuration, solver controls, units, seed, code fingerprint and compiler.
Actual overrides are exported. “Literature default” means the value was checked
in a source; it does not mean this combined model has been validated empirically.
Other biological parameters are provisional. Numerical settings are numerical
choices, not measurements. Model arrays and their ordering are public in
[`model.hpp`](../include/bnsim/model.hpp).

## Default biological values

| Parameters | Defaults in internal units | Origin/status |
| --- | --- | --- |
| QS basal, order A1/R1/RsaL/A2/R2 | all 1e-10 fmol/um^3/s | Legacy 1e-4, **assuming uM/s**, provisional |
| QS induced, same order | all 2e-9 fmol/um^3/s | Legacy 2e-3, **assuming uM/s**, provisional |
| KL1, KR1, KS, KS1, KL2, KL3, KL4 | [1e-15,1e-15,1e-12,1e-12,1e-9,1e-12,1e-12] fmol/um^3 | Legacy numeric thresholds interpreted as uM; provisional and exceptionally small |
| k0, k3, k6, k7, k10 | 0.01 /s | Legacy free-species decay; provisional |
| k1, k4, k8, k11 | 1e5 um^3/fmol/s | Legacy 0.1 interpreted as 1/(uM s), provisional |
| k2, k5, k9, k12 | 0.1 /s | Legacy unbinding/dimer dissociation; provisional |
| Membrane rates A1, A2, Q1, Q2, Q3 | all 0.5 /s | Legacy A1 rate extended to both signals and enabled inhibitors; provisional |
| QSI1 binding/unbinding | 1e5 um^3/fmol/s, 0.1 /s | Legacy mechanism, dimensional conversion; provisional |
| QSI2 kh, saturation g | 3e5 um^3/fmol/s, 2e6 um^3/fmol | Legacy 0.3, 2 in uM coordinates; provisional |
| QSI3 suppression threshold | 1e-15 fmol/um^3 | Legacy KL1; provisional |
| QS activation threshold | 1e-9 fmol/um^3 | New explicit closure, uses legacy growth QS scale; provisional |
| Growth umax, pmax, maintenance | 0.002, 0.0001, 0.0001 /s | Legacy cell constructor; provisional corrected closure |
| Growth Ks, Kc, yield | 1 fg/um^3, 1e-9 fmol/um^3, 0.85 fg/fg | Legacy values with explicit dimensions; provisional |
| Biomass/EPS density | 150, 75 fg/um^3 | Legacy biomass and commented EPS density; provisional |
| Initial biomass, EPS | 150, 0 fg | Example initialization choice |
| Division/death/shedding masses | 300, 1, 150 fg | New mass thresholds; provisional |
| Motility speed, Dr | 20 um/s, 0.062 rad^2/s | Checked against Sneddon 2012 Methods |
| Tumble gamma shape/scale | 4, 18.32 degrees | Checked against Sneddon 2012 Methods; conditional range [0,180] is a modeling choice |
| Adaptation kr, kb, m0, alpha, N | 0.2/s, 0.2/s, 1, 1.7, 6 | Legacy reduced model; provisional |
| Initial m | 2, dimensionless | Legacy initialization; concentration-dependent equilibration is not automatic |
| Ligand KA, KI | 0.003, 0.0000182 fmol/um^3 | Legacy 3 and 0.0182 **interpreted as mM**; provisional |
| Ya, Kd | 5e-6, 3.06e-6 fmol/um^3 | Legacy uM values; provisional in the reduced model |
| Motor g, w | 40, 1.3 /s | Legacy single-motor response; provisional |

The undocumented QS concentration scale cannot be recovered uniquely from this
repository. Interpreting its concentrations as uM preserves the legacy numeric
reaction structure after unit conversion. The extremely small resulting
thresholds are deliberately disclosed rather than presented as measured values.
Users should supply a dimensionally consistent, independently justified parameter
set for scientific prediction. `reference_qs` in the tests uses a separate,
well-scaled validation set; those values are not hidden simulator defaults.

## Published-model comparison

[Wei et al., Supplementary Table 2](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fsrep35136/MediaObjects/41598_2016_BFsrep35136_MOESM1_ESM.pdf)
labels zero-order concentration production and bimolecular binding constants
with /s. For concentration ODEs these require concentration/s and
1/(concentration s), respectively. This implementation corrects the dimensional
types but treats the absolute concentration scale as provisional. The table's
diffusion coefficient 1e-8 cm^2/s converts to 1 um^2/s, not the repository's 890.
Published decay, growth and EPS-density values also differ from the legacy code;
they are not silently substituted. These inconsistencies are reasons to avoid
claiming parameter calibration transfer.

[Sneddon et al., Methods](https://emonet.biology.yale.edu/assets/papers/past/sneddon_emonet_2012_paper.pdf)
was used for the speed, diffusion and angle defaults above. Its full receptor
composition, noisy adaptation and multiple flagella differ from this reduced
model. No per-step biochemical noise is retained.

## Numerical and scenario choices

| Configuration | Default / units | Status |
| --- | --- | --- |
| Extent; grid shape | [100,100,100] um; [10,10,10] | Numerical domain choice |
| Seed; workers | 0; 1 | Reproducibility controls |
| Closed walls; diffusion; mechanics | no flux; enabled; enabled | Boundary/numerical choices |
| Layered mode; boundary thickness | disabled; 10 um when enabled | Provisional geometric abstraction |
| Relative ODE tolerance | 1e-8 | Numerical choice |
| Molecular / mass ODE absolute tolerances | 1e-18 fmol / 1e-10 fg | Numerical choice |
| Minimum ODE step; attempts | 1e-14 s; 100000 | Failure limits |
| Motility maximum internal interval | 0.005 s | Numerical choice |
| Contact tolerance; maximum sweeps | 1e-8 um; 2000 | Numerical/failure choices |
| Maximum voxels / agents / transport substeps | 1000000 each | Resource limits |
| Species diffusion/decay | 0 unless configured | Explicit scenario input |

The chemotaxis example uses a prescribed, non-diffusing gradient with physical
voxel centers: `[0.01 + 0.0001*(Ly-y)] mM`. Its domain is [200,100,50] um and
grid [20,10,5]. The QS example uses a [60,60,60] um domain, [6,6,12] grid,
10-um boundary layer, bulk substrate 1 fg/um^3, zero bulk signal/inhibitor,
diffusion 30/18 um^2/s in boundary/biofilm, and signal decay 0.01/s. These are
illustrative scenario values, not the published experiment. All three inhibitors
are enabled but initially absent; dosing can be supplied through field/reservoir
setters. The initial QS cell placement is a lattice near z=1 um.
