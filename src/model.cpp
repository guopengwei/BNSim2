// BNSim, GPL-2.0. Corrected legacy reaction mechanisms; see docs/model.md.
#include "bnsim/model.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace BNSim {
Vec3 operator+(Vec3 a, Vec3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
Vec3 operator-(Vec3 a, Vec3 b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
Vec3 operator*(Vec3 a, double b) { return {a.x*b, a.y*b, a.z*b}; }
double dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec3 cross(Vec3 a, Vec3 b) { return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }
double norm(Vec3 a) { return std::hypot(a.x, a.y, a.z); }
Vec3 normalized(Vec3 a) { const auto n = norm(a); require_positive(n, "vector norm"); return a * (1/n); }
void require_nonnegative(double x, const std::string& name) {
    if (!std::isfinite(x) || x < 0) throw std::invalid_argument(name + " must be finite and nonnegative");
}
void require_positive(double x, const std::string& name) {
    if (!std::isfinite(x) || x <= 0) throw std::invalid_argument(name + " must be finite and positive");
}
std::uint64_t stream_seed(std::uint64_t seed, AgentId id) {
    auto x = seed + 0x9e3779b97f4a7c15ULL * (id + 1);
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
double from_micromolar(double x) { require_nonnegative(x, "micromolar concentration"); return x * 1e-6; }
double from_millimolar(double x) { require_nonnegative(x, "millimolar concentration"); return x * 1e-3; }
double from_grams_per_liter(double x) { require_nonnegative(x, "g/L concentration"); return x; }
double bimolecular_from_per_micromolar_second(double x) {
    require_nonnegative(x, "bimolecular rate");
    require_nonnegative(x * 1e6, "converted bimolecular rate");
    return x * 1e6;
}
void Species::validate() const {
    if (name.empty() || name.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-") != std::string::npos)
        throw std::invalid_argument("species name must contain only letters, digits, underscores or hyphens");
    if (quantity != Quantity::Mass && quantity != Quantity::Molar) throw std::invalid_argument("invalid species quantity");
    require_nonnegative(diffusion, name + " diffusion");
    require_nonnegative(biofilm_diffusion, name + " biofilm diffusion");
    require_nonnegative(decay, name + " decay");
    require_nonnegative(biofilm_decay, name + " biofilm decay");
    if (provenance.empty()) throw std::invalid_argument(name + " requires parameter provenance");
}
static void check_species(const std::vector<Species>& s, SpeciesId id, Quantity q, const char* name) {
    if (id >= s.size() || s[id].quantity != q) throw std::invalid_argument(std::string(name) + " requires a valid species with matching quantity units");
}
void QSParameters::validate(const std::vector<Species>& s) const {
    check_species(s, signal1, Quantity::Molar, "Las signal");
    check_species(s, signal2, Quantity::Molar, "Rhl signal");
    std::vector<SpeciesId> ids{signal1, signal2};
    for (const auto& id : inhibitors) if (id) { check_species(s, *id, Quantity::Molar, "inhibitor"); ids.push_back(*id); }
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) throw std::invalid_argument("QS species must be distinct");
    for (double x : basal) require_nonnegative(x, "QS basal synthesis");
    for (double x : induced) require_nonnegative(x, "QS induced synthesis");
    for (double x : threshold) require_positive(x, "QS threshold");
    for (double x : rate) require_nonnegative(x, "QS rate");
    for (double x : exchange) require_nonnegative(x, "QS exchange");
    require_nonnegative(inhibitor_binding, "QSI binding");
    require_nonnegative(inhibitor_unbinding, "QSI unbinding");
    require_nonnegative(signal_degradation, "QSI degradation");
    require_nonnegative(inhibitor_saturation, "QSI saturation");
    require_positive(receptor_suppression, "QSI suppression threshold");
    require_positive(activation_threshold, "QS activation threshold");
}
void GrowthParameters::validate(const std::vector<Species>& s) const {
    check_species(s, substrate, Quantity::Mass, "growth substrate");
    require_nonnegative(umax, "umax"); require_nonnegative(pmax, "pmax");
    require_nonnegative(maintenance, "maintenance");
    require_positive(substrate_half, "substrate half saturation");
    require_positive(qs_half, "QS half saturation");
    require_positive(yield, "yield");
    if (pmax > umax || yield > 1) throw std::invalid_argument("growth requires pmax <= umax and yield <= 1");
}
void MotilityParameters::validate(const std::vector<Species>& s) const {
    check_species(s, attractant, Quantity::Molar, "attractant");
    require_nonnegative(speed, "swim speed");
    require_nonnegative(rotational_diffusion, "rotational diffusion");
    require_nonnegative(kr, "kr"); require_nonnegative(kb, "kb");
    require_nonnegative(m0, "m0"); require_positive(alpha, "alpha"); require_positive(receptors, "receptors");
    require_positive(active_binding, "active ligand binding"); require_positive(inactive_binding, "inactive ligand binding");
    require_nonnegative(ya, "CheY scale"); require_positive(kd, "motor binding");
    require_nonnegative(g, "motor free energy"); require_nonnegative(motor_frequency, "motor frequency");
    if (g > 100) throw std::invalid_argument("motor free energy exceeds supported range (100)");
    require_positive(gamma_shape, "tumble gamma shape"); require_positive(gamma_scale_degrees, "tumble gamma scale");
    require_positive(max_step, "motility maximum step");
}
static double saturation(double x, double k) { return x == 0 ? 0 : 1 / (1 + k / x); }
ReactionRates qs_reactions(const QSParameters& p, const Chemicals& n, double v) {
    require_positive(v, "QS host volume");
    Chemicals c{};
    for (std::size_t i=0; i<chemical_count; ++i) { require_nonnegative(n[i], "intracellular amount"); c[i] = n[i]/v; }
    ReactionRates out;
    auto& d = out.derivative;
    const auto& k = p.rate;
    const auto& t = p.threshold;
    const double feedback = p.rsal_feedback ? 1/(1+c[RsaL]/t[3]) : 1;
    const double suppression = p.inhibitors[2] ? 1/(1+c[Q3]/p.receptor_suppression) : 1;
    const std::array<double,5> synthesis{
        p.basal[0]+p.induced[0]*saturation(c[C1],t[0])*feedback,
        (p.basal[1]+p.induced[1]*saturation(c[C1],t[1]))*suppression,
        p.basal[2]+p.induced[2]*saturation(c[C1],t[2]),
        p.basal[3]+p.induced[3]*saturation(c[C2],t[4]),
        p.basal[4]+p.induced[4]*saturation(c[C2],t[5])*saturation(c[C1],t[6])};
    const double binding1 = k[1]*c[A1]*c[R1]-k[2]*c[RA1];
    const double dimer1 = k[4]*c[RA1]*c[RA1]-k[5]*c[C1];
    const double binding2 = k[8]*c[A2]*c[R2]-k[9]*c[RA2];
    const double dimer2 = k[11]*c[RA2]*c[RA2]-k[12]*c[C2];
    const double competitive = p.inhibitors[0] ? p.inhibitor_binding*c[Q1]*c[R1]-p.inhibitor_unbinding*c[RQ1] : 0;
    const double degradation = p.inhibitors[1] ? p.signal_degradation*c[A1]*c[Q2]/(1+p.inhibitor_saturation*c[Q2]) : 0;
    d[A1] = synthesis[0]-k[0]*c[A1]-binding1-degradation;
    d[R1] = synthesis[1]-k[3]*c[R1]-binding1-competitive;
    d[RA1] = binding1-2*dimer1; d[C1] = dimer1;
    d[RsaL] = synthesis[2]-k[6]*c[RsaL];
    d[A2] = synthesis[3]-k[7]*c[A2]-binding2;
    d[R2] = synthesis[4]-k[10]*c[R2]-binding2;
    d[RA2] = binding2-2*dimer2; d[C2] = dimer2;
    d[Q1] = -competitive; d[RQ1] = competitive;
    for (auto& x : d) x *= v;
    out.signal_source = {synthesis[0]*v, synthesis[3]*v};
    out.signal_sink = {(k[0]*c[A1]+degradation)*v, k[7]*c[A2]*v};
    return out;
}
GrowthRates growth_rates(const GrowthParameters& p, double x, double s, double c1) {
    require_nonnegative(x, "growth biomass"); require_nonnegative(s, "growth substrate"); require_nonnegative(c1, "growth C1");
    const double f = saturation(s, p.substrate_half), allocation = p.pmax*saturation(c1,p.qs_half);
    const double production = p.umax*f*x, maintenance = p.maintenance*x;
    return {(p.umax-allocation)*f*x-maintenance, allocation*f*x, production/p.yield,
            maintenance, production*(1/p.yield-1)};
}
} // namespace BNSim
