# BNSim

BNSim is a C++17 simulator for bacterial chemotaxis, quorum sensing, growth,
and extracellular polymeric substance (EPS) production in three-dimensional
environments. It provides a reusable library and two runnable examples.

Chemistry and growth use deterministic integration. Motility uses seeded
random processes for swimming direction, rotational diffusion, and motor
switching. The biological models retain corrected legacy mechanisms;
**their provisional parameters have not been experimentally calibrated**.

## Build and run

Requirements: a C++17 compiler, CMake 3.16 or newer, platform thread support,
and Python 3 for the output tests. The library and examples have no third-party
numerical dependencies. Run these commands from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
(cd build && ctest --output-on-failure)
```

Run either example:

```sh
./build/chemotaxis_example --seed 0 --threads 8 --duration 2.35 --dt 0.1 --population 16 --output runs/chemotaxis
./build/qs_biofilm_example --seed 0 --threads 8 --duration 2.35 --dt 0.1 --population 16 --output runs/qs_biofilm
```

| Example | Scenario |
| --- | --- |
| `chemotaxis_example` | Swimming cells in a prescribed, static aspartate gradient initialized at physical voxel centers. |
| `qs_biofilm_example` | Substrate-dependent growth, Las/Rhl signaling, and EPS production with a boundary layer and bulk reservoir above the biomass. All three inhibitor mechanisms are enabled; initial inhibitor concentrations are zero. |

Both programs accept the same options and provide `--help`:

| Option | Default | Meaning |
| --- | --- | --- |
| `--seed N` | `0` | Unsigned 64-bit random seed. |
| `--threads N` | `1` | Worker count, from 1 to 4096. |
| `--duration SECONDS` | `1` | Finite, nonnegative end time. |
| `--dt SECONDS` | `0.1` | Finite, positive simulation step; the final step is shortened to reach the end time. |
| `--population N` | `20` | Initial cell count; zero is valid. Dense populations must fit mechanically in the example domain. |
| `--output DIRECTORY` | `runs` | Parent directory for a new run. |

Unknown, duplicate, missing, or invalid arguments produce an error and a nonzero
exit code. Each invocation creates a separate numbered run directory, so repeat
runs do not overwrite earlier results. To build without tests or their Python
dependency, add `-DBUILD_TESTING=OFF` when configuring CMake.

## Models and numerical behavior

- **Transport:** conservative finite-volume diffusion with harmonic interface
  diffusivities, automatic positivity-preserving substeps, and exponential
  first-order decay. Closed boxes have no-flux chemical boundaries; fixed
  concentration reservoirs are explicit and their exchanges are recorded.
- **Biochemistry:** Las and Rhl circuits, RsaL feedback, reversible binding,
  and the three legacy inhibitor mechanisms: receptor competition, signal
  degradation, and receptor-production suppression. Each voxel's extracellular
  pool and resident cells are integrated together with adaptive error control.
- **Growth and lifecycle:** substrate-limited biomass and QS-controlled EPS
  production, explicit maintenance and yield losses, growth dilution,
  conservative division, and EPS shedding by mass. Dead particles remain
  mechanically present and metabolically inactive.
- **Motility and mechanics:** a reduced receptor/adaptation model, continuous-time
  motor switching, stationary tumbles, rotational diffusion, radius-aware
  reflecting walls, and overlap relaxation with explicit convergence checks.

The [model specification](docs/model.md) defines the equations, stage order,
integration methods, boundaries, and numerical limitations. The layered biofilm
environment uses vertical `z` and a configured boundary layer above the greatest
particle surface height. This abstraction and the growth/EPS closure remain
modeling assumptions.

## Output files

Each run is written beneath the selected output directory as `run-00000000`,
`run-00000001`, and so on. Snapshots include the initial state and every completed
step. The eight-digit filename suffix is a **step index**, while file contents
store timestamps at full double precision.

| File | Contents |
| --- | --- |
| `metadata.json` | Configuration, seed, units, species definitions, code revision, source fingerprint, compiler, and build information. |
| `parameters.csv` | Actual agent parameter values, dimensions, provenance, and verification status. |
| `agents_00000000.csv` | Stable IDs, positions, biomass/EPS, radii, intracellular amounts and concentrations, and motility state. |
| `field_00000000.csv` | Full three-dimensional chemical fields, voxel centers and volumes, material properties, and reservoir targets. |
| `depth_average_00000000.csv` | Explicitly labeled depth-averaged concentrations in the horizontal plane. |
| `diagnostics_00000000.json` | Chemical inventory accounting, reservoir exchange, losses, lifecycle events, and solver diagnostics. |

Agent kinds are `0` for QS bacteria, `1` for chemotactic bacteria, and `2` for
EPS particles. The `alive` field distinguishes active cells from retained dead
material; EPS particles are also inactive. These files are analysis exports,
not restart checkpoints: they do not serialize random engine states.

## Use the library

Link the CMake target `bnsim` and include headers from `include/bnsim/`.
Configuration belongs to each `Universe`; species IDs follow their order in
the constructor. Enable cell networks explicitly through `AgentSpec`.

```cpp
#include <bnsim/simulation.hpp>
#include <iostream>

int main() {
    using namespace BNSim;

    SimulationConfig config;
    config.extent = {40, 30, 20}; // micrometers
    config.shape = {8, 6, 4};
    config.seed = 0;
    config.workers = 2;

    Universe world(config, {{"Aspartate", Quantity::Molar, 0, 0, 0, 0}});
    for (std::size_t i = 0; i < world.grids().size(); ++i)
        world.set_concentration(i, 0, from_millimolar(0.01));

    AgentSpec cell;
    cell.kind = AgentKind::ChemotacticBacterium;
    cell.position = {20, 15, 10};
    cell.motility = MotilityParameters{};
    cell.division = false;
    cell.shedding = false;
    const AgentId id = world.add_agent(cell);

    world.run_until(1.25, 0.1);
    const Agent& result = world.agent(id);
    std::cout << "Cell " << result.id << " at " << world.time() << " s\n";
}
```

State access is immutable. Retain agent IDs and reacquire views after mutations,
which can invalidate references. A failed step restores the previous state,
including time, diagnostics, and random streams. Internal workers are managed
by the simulator; application code must not call the same `Universe`
concurrently. Separate instances share no mutable simulation state.

Internal quantities are explicit:

| Quantity | Unit |
| --- | --- |
| Length and time | micrometers (`um`) and seconds (`s`) |
| Biomass, EPS, and substrate mass | femtograms (`fg`) |
| Signaling chemical amount | femtomoles (`fmol`) |
| Concentration | `fg/um^3` or `fmol/um^3`, according to species |

Use the conversion helpers at input boundaries. For example, `1 mM` is
`1e-3 fmol/um^3`, and `1 g/L` is `1 fg/um^3`. Intracellular state stores amounts;
concentrations are derived from the current cell volume.

## Validation and reproducibility

The test suite covers analytical transport, conservation, independent
biochemical reference trajectories, lifecycle accounting, motility statistics,
collisions, invalid inputs, and example output files. For a fixed seed, complete
simulation states are checked for exact equality across 1, 2, and 8 workers on
the same build and platform. Cross-platform comparisons require numerical
tolerances.

The [recorded local validation](docs/validation_results.md) passed all seven
CTest entries in Release, AddressSanitizer/UndefinedBehaviorSanitizer with leak
detection, and ThreadSanitizer builds. The
[CI workflow](.github/workflows/ci.yml) defines Linux and macOS builds and Linux
sanitizer checks. See [validation instructions](docs/validation.md) to reproduce
these checks or regenerate the optional SciPy reference solutions.

Numerical verification establishes properties of the implementation. Biological
predictive accuracy still requires suitable experimental data and calibration.

## Documentation and migration

| Document | Purpose |
| --- | --- |
| [Model specification](docs/model.md) | Equations, numerical methods, boundaries, and limitations. |
| [Parameter reference](docs/parameters.md) | Units, defaults, provenance, provisional values, and dimensional corrections. |
| [Issue register](docs/issue_register.md) | Defect evidence, corrections, regression checks, modeling assumptions, and open empirical validation. |
| [Validation guide](docs/validation.md) | Test, sanitizer, and independent reference generation commands. |
| [Validation results](docs/validation_results.md) | Observed outcomes, toolchain, and scope of verification. |
| [Migration guide](docs/migration.md) | Mapping from the original global configuration and raw-pointer API. |

This version intentionally changes the legacy API and numerical behavior.
Historical trajectories are not preserved. The original implementation remains
available in Git history; use the migration guide when adapting existing code.

## Attribution and license

Originally developed by the CMU SLD Group. Distributed under the
[GNU General Public License, version 2](LICENSE).

Please cite the original BNSim publication when using this work:

Wei, Guopeng, Paul Bogdan, and Radu Marculescu. “Efficient Modeling and Simulation
of Bacteria-Based Nanonetworks with BNSim.” *IEEE Journal on Selected Areas in
Communications* 31.12 (2013): 868–878.

The model documentation also compares mechanisms with
[Sneddon et al. (2012)](https://doi.org/10.1073/pnas.1113706109) and
[Wei et al. (2016)](https://doi.org/10.1038/srep35136), and records where this
implementation differs from those publications.
