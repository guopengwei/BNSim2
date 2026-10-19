# Scientific correctness issue register

Baseline: Git revision `203f29c`. Historical paths below refer to that revision,
which can be inspected with `git show 203f29c:path`. The supplied audit plan
reported compiler, diagnostic-program, ASan and TSan evidence. Its temporary
programs/logs were not present in this checkout; those runtime observations are
identified as **reported audit evidence**, not claimed as rerun artifacts.
The formulas and implementation defects were independently confirmed by reading
the baseline. Current regression checks are committed in
[`tests/test_main.cpp`](../tests/test_main.cpp) and
[`tests/check_examples.py`](../tests/check_examples.py). The validation report
records actual post-correction execution.

Status “corrected” refers to the implementation and named checks, not empirical
validation. Every defect identified in the supplied plan is mapped below.

## Implementation defects

| ID | Evidence at baseline | Correction | Regression check | Status |
| --- | --- | --- | --- | --- |
| I01 | `core/universe.cpp`, decay computes `(1-k)*C*dt`; audit reports 1 -> 0.1 at dt=0.1 even with k=0 | Exact exponential, independent of diffusion flag | `zero_decay_invariance`, `analytic_decay_without_diffusion` | corrected |
| I02 | `environment_thread` updates neighbors in place, visits both interface directions and clips at zero; audit reports material creation | Immutable face fluxes, one face once, multidimensional positive substeps | `closed_conservation_and_positivity`, `harmonic_interfaces_and_impermeable_material` | corrected |
| I03 | Diffusion arithmetic-average interface D is inconsistent with a finite-volume discontinuity | Harmonic interface D; impermeable side gives zero flux | `harmonic_interfaces_and_impermeable_material` | corrected |
| I04 | `Grid` defaults to `bulk`, conflating material and fixed amounts; bulk sources unaccounted | Independent material/reservoir state and signed transfer ledger | `reservoir_and_decay_accounting`, `explicit_withdrawals_after_reservoir_release` | corrected |
| I05 | `prepare_multithreading` indexes `threadNumber-2`; audit reports heap overflow with one worker | Bounded standard-thread partitioning including empty work | `empty_world_and_final_time`, sanitizer suites, CLI population 0/threads 8 | corrected |
| I06 | Shared field reads/writes, global RNG and concurrent insertion; audit reports TSan races | Disjoint results, deterministic reductions, staged births, per-ID RNG | `identical_across_worker_counts`, ThreadSanitizer | corrected |
| I07 | Thread creation failure returns without reliable joins; exceptions are swallowed in position/collision loops | Join all launched workers; deterministic exception propagation; transactional steps | `worker_exceptions_join`, `integrator_failure_is_transactional`, `transport_failure_rolls_back` | corrected |
| I08 | `spacegrid.cpp:37` uses z*y*z; integer grid spacing truncates geometry | Floating spacings, x*y*z volume, checked sizes | `geometry_and_units`, noncubic conservation/convergence | corrected |
| I09 | `Agent::registerGridPos/updateGridPos` round toward nearest voxel and add pending motion after position commit | Floor committed position only; membership rebuilt between stages | `movement_membership_uses_committed_position_once`, `geometry_and_units` | corrected |
| I10 | `consumeChemical` subtracts an ambiguous quantity and clamps; cells independently claim substrate | Amount/concentration separation and joint voxel integration | `substrate_limited_growth_and_losses`, `shared_voxel_limited_exchange` | corrected |
| I11 | Chemotaxis adaptation multiplies timestep inside and outside kr/kb | Deterministic RK4 adaptation with one time factor | `attractant_response_and_adaptation` | corrected |
| I12 | Motor rates compared directly to uniform probabilities, without dt | Persistent exponential hazard, continuous event times, multiple switches | `continuous_time_motor_statistics` | corrected |
| I13 | Uniform nonnegative rotational angles and normalized cube samples give incorrect rotational/heading statistics | Isotropic sphere initialization and independent tangent Gaussians | `isotropic_initial_headings`, `gaussian_rotation_correlation` | corrected |
| I14 | Tumble rotation repeated every step; pi/360 conversion; negative shifted angles possible | One bounded gamma angle per tumble exit, pi/180, uniform azimuth | `bounded_gamma_statistics`, `stationary_tumbles_and_single_reorientation` | corrected |
| I15 | `ChemotacticBacteria::update` always adds swimming displacement | Integrate displacement only in run intervals | `stationary_tumbles_and_single_reorientation`, `ensemble_drift_up_attractant_gradient` | corrected |
| I16 | `SimpleMetabolism` has uninitialized members and uses `_u_PG` without assignment | Fully initialized validated parameters and explicit QS allocation | `growth_dilution_and_radius`, `eps_production_and_allocation_budget`, UBSan | corrected |
| I17 | QSLux reads inhibitor IDs but leaves sensed values at zero; disabled exchanges | Enabled species mappings and exchange for all three inhibitors | `membrane_exchange_analytical`, `circuits_feedback_activation_and_inhibitors` | corrected |
| I18 | QSI1 unbinding fails to return free receptor in the R1 equation | Stoichiometric reaction flux updates all partners | `reversible_binding_invariants` | corrected |
| I19 | Rhl, RsaL feedback and activation logic are commented out | Restore both circuits, negative feedback and concentration-based active state | `circuits_feedback_activation_and_inhibitors`, reference trajectory | corrected |
| I20 | QS exchange uses radius cubed, conflating capsule and cell volume; only A1 exchanged | Amount-based cell interior, both signals and inhibitors, equal/opposite exchange | `membrane_exchange_analytical`, `growth_dilution_and_radius` | corrected |
| I21 | Per-step biochemical perturbations depend on dt and schedule | Deterministic adaptive kinetics; only motility stochastic | `independent_reference_trajectories`, `identical_across_worker_counts` | corrected |
| I22 | Growth leaves radius unchanged; `pow(V,0.33)` approximates cube root; yield/EPS ignored | Density-based volumes, `cbrt`, implemented growth/EPS and explicit losses | `growth_dilution_and_radius`, `eps_production_and_allocation_budget`, starvation test | corrected |
| I23 | Division gives each cell 45% of mass and new default network concentrations/parameters | Exact halves of all masses/amounts; inherited parameters and concentrations | `conservative_division_and_inheritance` | corrected |
| I24 | EPS release treats capsule thickness as particle radius; EPS mass creation disabled | Mass transfer and density-derived radius | `shedding_uses_mass_and_density`, `eps_production_and_allocation_budget` | corrected |
| I25 | Lifecycle inactivity is disconnected from actual growth geometry and lacks integrated accounting | Mass death threshold; inactive particles retain mechanics and amounts | `dead_particles_inactive_and_retained`, mixed-world checks | corrected |
| I26 | Collision search uses only six face neighbors; normalized zero separation gives NaNs | Full radius-aware voxel range, deterministic coincident direction, convergence/failure | `coincident_diagonal_and_large_contacts`, `confinement_failure_rolls_back` | corrected |
| I27 | Position clamp ignores radius and loses reflecting trajectories; axis helpers confuse vertical axes | Radius-aware reflection/projected contacts; z-based layers | `wall_contacts_and_reflection`, `vertical_layers_and_reservoir_reclassification` | corrected |
| I28 | Global mutable CONFIG, overlapping raw ownership, leaked thread arrays, unsafe custom container growth | Instance configuration and standard value containers, single ownership | `independent_universes_and_ids`, ASan/leak suite | corrected |
| I29 | Unchecked identifiers/coordinates/size products/parameters; constructor registration uses uninitialized movement | Checked construction, mapping, bounds and resource limits; explicit initial state | `reject_invalid_inputs`, `geometry_and_units` | corrected |
| I30 | Conflicting `PAQS_H_` guards, undefined methods, nonstandard stack arrays; QS main commented out | New public API, complete methods, portable C++17 library and two executables | Warning-as-error builds, all CTest groups and both CLI runs | corrected |
| I31 | `sprintf` treats output path as a format string into 100 bytes; int-time filenames overwrite | Filesystem paths, unique run directories, checked streams, step indices and 17-digit times | `examples_and_outputs`, long percent-containing directory | corrected |
| I32 | “Concentration2D” sums depth without labeling; exports cast mixed agents to QS unconditionally | Full 3D fields, labeled depth averages and typed per-agent output | `examples_and_outputs`, mixed-agent numerical tests | corrected |
| I33 | Chemotaxis gradient uses voxel index as a physical coordinate; loop overshoots final duration | Physical voxel centers and shortened final step | CLI gradient check, `empty_world_and_final_time` | corrected |
| I34 | No run seed/configuration/provenance/version or reproducible validation | Complete metadata, numerical references, tests, Linux/macOS and sanitizer CI | `examples_and_outputs`, `seed_controls_motility_only`, reference provenance hashes | corrected |
| I35 | Found during implementation: non-SSP DP internal stages can remain negative for zero-state cascades at any h | Error-controlled positive SSPRK3 startup, then DP5(4); no clipping | `zero_state_reaction_cascade_startup`, initially empty QS examples | corrected |

## Modeling assumptions (explicit, not implementation bugs)

| ID | Assumption/evidence | Implementation and checks | Remaining limitation |
| --- | --- | --- | --- |
| M01 | Legacy chemical concentrations lack a documented unit scale; published production/binding dimensions are incomplete | Explicit assumed uM QS scale and dimensionally converted rates; parameter records; unit tests | Absolute biological scale requires independent justification; see `parameters.md` |
| M02 | Legacy growth closure was unfinished; requested corrected closure fixes allocation/yield interpretation | Defined biomass/EPS/substrate equations and maintenance/yield ledgers; conservation tests | No explicit carbon/amino-acid cost of QS synthesis; not an elemental physiological model |
| M03 | Reduced receptor model differs from Sneddon's full stochastic multi-motor model | Deterministic receptor adaptation, single-motor hazard, checked motility statistics and ensemble response | No claim to reproduce multi-flagella coordination or experimental drift magnitudes |
| M04 | Requested legacy inhibitor identities differ from parts of Wei 2016 | QSI1 competition, QSI2 intracellular degradation, QSI3 receptor-production suppression | Published LasI inhibition/extracellular enzyme experiments are not directly reproduced |
| M05 | Legacy Rhl/RsaL and complex stability are model choices | Restored commented mechanisms; binding/reference tests | QS circuitry and parameter values need biological calibration |
| M06 | Uniform layer height and fixed voxel carrier volume approximate the environment | Greatest z surface, voxel-center layers, reservoir accounting | No resolved pore volume, local interface geometry or detailed flow |
| M07 | Spherical particles and geometric exclusion approximate mechanics | Full-neighbor contact resolution and clear convergence failure | No elastic, hydrodynamic or adhesion-force calibration |
| M08 | Transport Euler steps, held ligand/rates, stage-boundary lifecycle events and startup integration introduce discretization errors | Refinement checks, reference solutions, exact fixed-rate dwell statistics | New scenarios require independent timestep/grid/tolerance studies |
| M09 | Random library implementations vary between platforms | Identical state on 1/2/8 workers in the same build; seeds/fingerprints recorded | Cross-platform reproducibility is numerical/statistical, not bitwise |

## Missing empirical validation (open research work)

| ID | Missing evidence | What is established here | What would close it |
| --- | --- | --- | --- |
| V01 | No matched concentration/time-course dataset for both QS circuits and all inhibitors | Correct stoichiometry, qualitative responses and numerical-reference agreement | Independent quantitative calibration and held-out dose/time-course data |
| V02 | No matched growth/EPS/substrate/maintenance measurements | Conservative corrected closure and deterministic integration | Strain/medium-specific biomass, EPS and substrate balances with uncertainty |
| V03 | No matched single-cell motor/trajectory dataset for the reduced model | Specified distributions, correlation, adaptation and directional response | Trajectory and motor-dwell validation across gradients and conditions |
| V04 | No biofilm geometry/transport/adhesion calibration dataset | Numerically defined layer/reservoir and collision mechanics | Spatially resolved transport, morphology and mechanical measurements |

These open empirical items do not block running or numerically validating the
simulator. They block claims of calibrated biological prediction. No historical
fit or experimental validation has been inferred from passing software tests.
