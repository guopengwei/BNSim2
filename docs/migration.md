# Migrating from revision 203f29c

The unsafe legacy API and generated Eclipse `Debug/` build files were removed.
The original implementation remains available in Git history. The change is
intentionally API-breaking; there are no wrappers that preserve implicit global
state, inconsistent quantities or historical trajectories.

| Old operation | New operation |
| --- | --- |
| Assign `CONFIG::*` statics | Populate an instance `SimulationConfig` before construction |
| Default `Universe`, then register `MoleculeInfo*` | `Universe(config, vector<Species>)`; species order supplies stable IDs |
| Global `CONFIG::universe` inside a network | Pass parameters, amounts, volume and local environment explicitly |
| `new QSBacteria` / raw `addAgent` | Fill `AgentSpec`, enable `qs` and/or `growth`, call `add_agent` |
| `new ChemotacticBacteria` | Set kind and `motility` in `AgentSpec` |
| `new EPS` with capsule thickness | EPS spec with biomass 0, transferred EPS mass and density |
| Mutable `Agent*` or `Grid*` access | `agent(id)`, `agents()`, `grids()` return immutable views |
| `getGrid(x,y,z)` | `index(x,y,z)` and `grids()[index]` |
| `setConc` / `consumeChemical` | Checked `set_concentration` at setup; consumption is part of the coupled ODE |
| Layer enum doubles as fixed-concentration flag | Independent `set_material` and `set_reservoir` |
| `evolute()` with implicit timestep | `step(dt)` or `run_until(end_time, dt)` |
| Radius determines mutable cell mass inconsistently | Supply mass/density; spherical radius is always derived |
| Per-cell network concentrations | `initial_concentrations` at creation; amount state thereafter |
| Global `rand()` | Universe seed and ID-derived per-agent streams |
| Legacy exporters | `RunWriter(output_parent, world, duration, dt, scenario)` and `snapshot` |

Headers are under `include/bnsim/` and all symbols remain in namespace `BNSim`.
Link the reusable CMake target `bnsim`. Both examples are separately linked from
`src/apps/example.cpp`, so there is no duplicate or disabled `main`.

An initially unspecified network is absent. A QS cell should configure both
signal species; growth needs a mass-based substrate mapping. Setting only
`kind=QSBacterium` labels a particle but does not implicitly create a network.
Inhibitor IDs are optional; missing ones disable their mechanisms. Initial cell
surfaces must fit within the domain. Overlapping initial cells are allowed and
resolved at the next mechanical stage, which can fail for impossible packings.

Convert legacy concentrations explicitly. Use `from_micromolar`,
`from_millimolar`, `from_grams_per_liter`, and the bimolecular-rate conversion
helper as appropriate. Do not mix molecular counts, molarity and femtograms.
Thresholds formerly expressed as radii must be converted to mass with the
intended density; default mass thresholds are new provisional choices.

References to vectors or agents can be invalidated by a successful step or a
creation/removal operation. Retain IDs and reacquire immutable views. Do not
invoke mutations concurrently on the same `Universe`. `remove_agent` is an
explicit external removal and updates inventory ledgers. Death itself retains
the particle and all chemical amounts.

Outputs use step indices, so `0.1` and `0.2` seconds no longer overwrite the same
integer-time file. Field CSVs include voxel centers, volumes, species units and
reservoir targets; projections explicitly represent depth averages. Per-cell
files expose amounts and concentrations with distinct headers. Recreate a
scenario from code, its seed and configuration to repeat a run; exported CSVs
do not contain serialized RNG states and are not restart files.
