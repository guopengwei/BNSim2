// BNSim, GPL-2.0. See LICENSE and docs/model.md.
#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace BNSim {
constexpr double pi = 3.1415926535897932384626433832795;
using AgentId = std::uint64_t;
using SpeciesId = std::size_t;
struct Vec3 {
    double x = 0, y = 0, z = 0;
    double& operator[](std::size_t i) { return i == 0 ? x : (i == 1 ? y : z); }
    double operator[](std::size_t i) const { return i == 0 ? x : (i == 1 ? y : z); }
};
Vec3 operator+(Vec3 a, Vec3 b);
Vec3 operator-(Vec3 a, Vec3 b);
Vec3 operator*(Vec3 a, double b);
double dot(Vec3 a, Vec3 b);
Vec3 cross(Vec3 a, Vec3 b);
double norm(Vec3 a);
Vec3 normalized(Vec3 a);
void require_nonnegative(double value, const std::string& name);
void require_positive(double value, const std::string& name);
std::uint64_t stream_seed(std::uint64_t seed, AgentId id);

enum class Quantity { Mass, Molar };
enum class Material { Boundary, Biofilm };
enum class AgentKind { QSBacterium, ChemotacticBacterium, EPS };
// 1 M = 1 fmol/um^3; 1 g/L = 1 fg/um^3.
double from_micromolar(double value);
double from_millimolar(double value);
double from_grams_per_liter(double value);
double bimolecular_from_per_micromolar_second(double value);

struct Species {
    std::string name;
    Quantity quantity = Quantity::Molar;
    double diffusion = 0, biofilm_diffusion = 0; // um^2/s
    double decay = 0, biofilm_decay = 0;         // 1/s
    std::string provenance = "user supplied; unverified";
    void validate() const;
};

enum Chemical : std::size_t { A1, R1, RA1, C1, RsaL, A2, R2, RA2, C2, Q1, RQ1, Q2, Q3, chemical_count };
inline constexpr std::array<const char*, chemical_count> chemical_names{
    "A1", "R1", "RA1", "C1", "RsaL", "A2", "R2", "RA2", "C2", "QSI1", "RQSI1", "QSI2", "QSI3"};
using Chemicals = std::array<double, chemical_count>;

struct QSParameters {
    SpeciesId signal1 = 0, signal2 = 1;
    std::array<std::optional<SpeciesId>, 3> inhibitors{};
    // Order: A1, R1, RsaL, A2, R2. Concentration rates in fmol/um^3/s.
    std::array<double, 5> basal{1e-10, 1e-10, 1e-10, 1e-10, 1e-10};
    std::array<double, 5> induced{2e-9, 2e-9, 2e-9, 2e-9, 2e-9};
    // KL1, KR1, KS, KS1, KL2, KL3, KL4, all fmol/um^3.
    std::array<double, 7> threshold{1e-15, 1e-15, 1e-12, 1e-12, 1e-9, 1e-12, 1e-12};
    // Legacy k[0..12]; bimolecular entries 1,4,8,11 in um^3/fmol/s, others 1/s.
    std::array<double, 13> rate{.01, 1e5, .1, .01, 1e5, .1, .01, .01, 1e5, .1, .01, 1e5, .1};
    std::array<double, 5> exchange{.5, .5, .5, .5, .5}; // 1/s: A1,A2,Q1,Q2,Q3
    double inhibitor_binding = 1e5, inhibitor_unbinding = .1;
    double signal_degradation = 3e5; // um^3/fmol/s
    double inhibitor_saturation = 2e6; // um^3/fmol
    double receptor_suppression = 1e-15; // fmol/um^3
    double activation_threshold = 1e-9;  // fmol/um^3 of C1
    bool rsal_feedback = true;
    void validate(const std::vector<Species>& species) const;
};

struct GrowthParameters {
    SpeciesId substrate = 0;
    double umax = .002, pmax = .0001, maintenance = .0001; // 1/s
    double substrate_half = 1; // fg/um^3
    double qs_half = 1e-9;     // fmol/um^3
    double yield = .85;       // fg biomass+EPS / fg substrate
    void validate(const std::vector<Species>& species) const;
};

struct MotilityParameters {
    SpeciesId attractant = 0;
    double speed = 20;                // um/s
    double rotational_diffusion = .062; // rad^2/s
    double kr = .2, kb = .2;         // 1/s
    double m0 = 1, alpha = 1.7, receptors = 6;
    double active_binding = .003, inactive_binding = .0000182; // fmol/um^3, legacy mM
    double ya = 5e-6, kd = 3.06e-6; // fmol/um^3, legacy uM
    double g = 40, motor_frequency = 1.3; // dimensionless, 1/s
    double gamma_shape = 4, gamma_scale_degrees = 18.32;
    double max_step = .005; // s; rate freezing and rotational integration resolution
    void validate(const std::vector<Species>& species) const;
};

struct MotilityState {
    Vec3 heading{1, 0, 0};
    double methylation = 2;
    bool running = true;
    double hazard_remaining = 1; // unit exponential cumulative hazard
    double run_seconds = 0, tumble_seconds = 0;
    std::uint64_t switches = 0, reorientations = 0;
};

struct AgentSpec {
    AgentKind kind = AgentKind::QSBacterium;
    Vec3 position{};
    double biomass = 150, eps = 0;       // fg
    double biomass_density = 150, eps_density = 75; // fg/um^3
    double division_mass = 300, death_mass = 1, shedding_mass = 150; // fg
    bool division = true, shedding = true, alive = true;
    std::optional<QSParameters> qs;
    std::optional<GrowthParameters> growth;
    std::optional<MotilityParameters> motility;
    Chemicals initial_concentrations{}; // fmol/um^3, cell interior only
    std::optional<Vec3> initial_heading;
    double initial_methylation = 2;
    bool initial_running = true;
};

struct Agent {
    AgentId id = 0;
    AgentKind kind = AgentKind::QSBacterium;
    Vec3 position{};
    double biomass = 0, eps = 0, biomass_density = 150, eps_density = 75;
    double division_mass = 300, death_mass = 1, shedding_mass = 150;
    bool division = true, shedding = true, alive = true, activated = false;
    std::optional<QSParameters> qs;
    std::optional<GrowthParameters> growth;
    std::optional<MotilityParameters> motility;
    Chemicals amounts{}; // fmol, never concentrations
    MotilityState motion;
    std::mt19937_64 random;
    double volume() const { return biomass / biomass_density; }
    double total_volume() const { return volume() + eps / eps_density; }
    double cell_radius() const { return std::cbrt(3 * volume() / (4 * pi)); }
    double radius() const { return std::cbrt(3 * total_volume() / (4 * pi)); }
    double concentration(Chemical c) const { return volume() > 0 ? amounts[c] / volume() : 0; }
};

struct ReactionRates {
    Chemicals derivative{}; // fmol/s
    std::array<double, 2> signal_source{}, signal_sink{}; // fmol/s, includes free signal decay
};
// Pure amount-based network evaluation with explicit host volume and parameters.
ReactionRates qs_reactions(const QSParameters&, const Chemicals& amounts, double volume);
struct GrowthRates { double biomass, eps, substrate, maintenance_loss, yield_loss; };
GrowthRates growth_rates(const GrowthParameters&, double biomass, double substrate_concentration, double c1);
double receptor_activity(const MotilityParameters&, double methylation, double ligand);
std::array<double, 2> motor_rates(const MotilityParameters&, double activity); // run->tumble, tumble->run
Vec3 isotropic_direction(std::mt19937_64&);
void rotational_increment(Vec3&, double diffusion, double dt, std::mt19937_64&);
double tumble_angle(const MotilityParameters&, std::mt19937_64&);
void reflect_position(Vec3& position, Vec3& heading, Vec3 extent, double radius);
void advance_motility(Agent&, double ligand, double dt, Vec3 extent);
} // namespace BNSim
