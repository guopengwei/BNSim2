// BNSim, GPL-2.0. See LICENSE.
#pragma once
#include "model.hpp"
#include <functional>

namespace BNSim {
struct IntegrationOptions {
    double relative_tolerance = 1e-8;
    double molecular_absolute_tolerance = 1e-18; // fmol
    double mass_absolute_tolerance = 1e-10;      // fg
    double minimum_step = 1e-14; // s
    std::size_t max_attempts = 100000;
};
struct SimulationConfig {
    Vec3 extent{100, 100, 100}; // um
    std::array<std::size_t, 3> shape{10, 10, 10};
    std::size_t workers = 1;
    std::uint64_t seed = 0;
    bool diffusion = true, mechanics = true, layered = false;
    double boundary_layer = 10; // um above greatest particle surface z
    std::vector<double> bulk_concentrations; // species units per um^3
    double mechanical_tolerance = 1e-8; // um
    std::size_t mechanical_iterations = 2000;
    std::size_t max_voxels = 1000000, max_agents = 1000000, max_transport_substeps = 1000000;
    IntegrationOptions integration;
    void validate(std::size_t species_count) const;
};

struct Grid {
    std::vector<double> amounts; // fg or fmol according to species
    Material material = Material::Boundary;
    bool bulk = false;
    std::vector<std::optional<double>> fixed_concentrations; // explicit reservoir, independent of material
    std::vector<AgentId> members; // non-owning IDs, sorted
};
struct Diagnostics {
    std::vector<double> external, reservoir, decay, reaction_source, reaction_sink;
    double external_particle_mass = 0, maintenance_loss = 0, yield_loss = 0;
    double biomass_height = 0, max_overlap = 0;
    std::size_t transport_substeps = 0, ode_accepted = 0, ode_rejected = 0, ode_startup = 0, mechanical_iterations = 0;
    std::uint64_t births = 0, deaths = 0, shed_particles = 0;
};

class Universe {
public:
    Universe(SimulationConfig config, std::vector<Species> species);
    const SimulationConfig& config() const { return config_; }
    const std::vector<Species>& species() const { return species_; }
    const std::vector<Grid>& grids() const { return grids_; }
    const std::vector<Agent>& agents() const { return agents_; }
    const Diagnostics& diagnostics() const { return diagnostics_; }
    double time() const { return time_; }
    std::uint64_t step_index() const { return step_; }
    Vec3 spacing() const { return spacing_; }
    double voxel_volume() const { return volume_; }
    SpeciesId species_id(const std::string& name) const;
    std::size_t index(std::size_t x, std::size_t y, std::size_t z) const;
    std::size_t index(Vec3 position) const;
    Vec3 center(std::size_t voxel) const;
    double concentration(std::size_t voxel, SpeciesId species) const;
    const Agent& agent(AgentId id) const;
    void set_concentration(std::size_t voxel, SpeciesId species, double concentration);
    void set_material(std::size_t voxel, Material material);
    void set_reservoir(std::size_t voxel, SpeciesId species, std::optional<double> concentration);
    AgentId add_agent(const AgentSpec& spec);
    void remove_agent(AgentId id); // explicit external removal; death never removes a particle
    void step(double dt); // strong exception guarantee, including RNGs and diagnostics
    void run_until(double end_time, double dt);
    double chemical_total(SpeciesId species) const; // extracellular plus intracellular signal moieties
    double particle_mass() const;
    double expected_chemical_total(SpeciesId species) const;
    void validate_state() const;
private:
    struct Face { std::size_t a, b; double inverse_spacing_squared; };
    SimulationConfig config_;
    std::vector<Species> species_;
    std::vector<Grid> grids_;
    std::vector<Agent> agents_;
    std::vector<Face> faces_;
    Diagnostics diagnostics_;
    Vec3 spacing_;
    double volume_ = 0, time_ = 0;
    AgentId next_id_ = 0;
    std::uint64_t step_ = 0;
    std::size_t agent_index(AgentId id) const;
    AgentId add_impl(const AgentSpec&);
    std::optional<double> reservoir(std::size_t voxel, SpeciesId species) const;
    void rebuild_membership();
    void update_layers();
    void transport(double dt);
    void decay(double dt);
    void biology(double dt);
    void movement(double dt);
    void lifecycle();
    void relax();
    void advance(double dt);
};
} // namespace BNSim
