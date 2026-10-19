// BNSim, GPL-2.0. Instance-owned state and transactional stage orchestration.
#include "bnsim/simulation.hpp"
#include "numerics.hpp"
#include <algorithm>
#include <limits>
#include <set>
#include <stdexcept>

namespace BNSim {
void SimulationConfig::validate(std::size_t count) const {
    if (workers == 0 || workers > 4096) throw std::invalid_argument("workers must be in [1,4096]");
    if (!max_voxels || !max_agents || !max_transport_substeps || !mechanical_iterations)
        throw std::invalid_argument("resource and iteration limits must be positive");
    std::size_t voxels=1;
    for (std::size_t a=0; a<3; ++a) {
        require_positive(extent[a], "world extent");
        if (shape[a] == 0 || shape[a] > max_voxels/voxels) throw std::invalid_argument("grid shape exceeds voxel limit or is zero");
        voxels*=shape[a];
        require_positive(extent[a]/static_cast<double>(shape[a]), "grid spacing");
    }
    if (count && voxels > std::numeric_limits<std::size_t>::max()/count/sizeof(double))
        throw std::invalid_argument("chemical field allocation overflows size_t");
    require_nonnegative(boundary_layer, "boundary layer thickness");
    require_positive(mechanical_tolerance, "mechanical tolerance");
    if (layered && bulk_concentrations.size() != count) throw std::invalid_argument("layered worlds require one bulk concentration per species");
    if (!layered && !bulk_concentrations.empty()) throw std::invalid_argument("bulk concentrations require layered=true");
    for (double x : bulk_concentrations) require_nonnegative(x, "bulk concentration");
    require_positive(integration.relative_tolerance, "ODE relative tolerance");
    require_positive(integration.molecular_absolute_tolerance, "ODE molecular tolerance");
    require_positive(integration.mass_absolute_tolerance, "ODE mass tolerance");
    require_positive(integration.minimum_step, "ODE minimum step");
    if (!integration.max_attempts) throw std::invalid_argument("ODE attempt limit must be positive");
}
Universe::Universe(SimulationConfig config, std::vector<Species> species)
    : config_(std::move(config)), species_(std::move(species)) {
    config_.validate(species_.size());
    std::set<std::string> names;
    for (const auto& s : species_) { s.validate(); if (!names.insert(s.name).second) throw std::invalid_argument("duplicate species: "+s.name); }
    for (std::size_t a=0; a<3; ++a) spacing_[a]=config_.extent[a]/static_cast<double>(config_.shape[a]);
    volume_=spacing_.x*spacing_.y*spacing_.z;
    require_positive(volume_, "voxel volume");
    const auto n=config_.shape[0]*config_.shape[1]*config_.shape[2], m=species_.size();
    grids_.resize(n);
    for (auto& g : grids_) { g.amounts.resize(m); g.fixed_concentrations.resize(m); }
    diagnostics_.external.resize(m); diagnostics_.reservoir.resize(m); diagnostics_.decay.resize(m);
    diagnostics_.reaction_source.resize(m); diagnostics_.reaction_sink.resize(m);
    for (std::size_t x=0; x<config_.shape[0]; ++x)
        for (std::size_t y=0; y<config_.shape[1]; ++y)
            for (std::size_t z=0; z<config_.shape[2]; ++z) {
                std::array<std::size_t,3> coordinate{x,y,z};
                for (std::size_t a=0; a<3; ++a) if (coordinate[a]+1 < config_.shape[a]) {
                    auto neighbor=coordinate; ++neighbor[a];
                    const double metric=1/spacing_[a]/spacing_[a];
                    require_positive(metric, "finite-volume face metric");
                    faces_.push_back({index(x,y,z),index(neighbor[0],neighbor[1],neighbor[2]),metric});
                }
            }
    update_layers();
}
SpeciesId Universe::species_id(const std::string& name) const {
    for (std::size_t i=0; i<species_.size(); ++i) if (species_[i].name == name) return i;
    throw std::out_of_range("unknown species: "+name);
}
std::size_t Universe::index(std::size_t x, std::size_t y, std::size_t z) const {
    if (x>=config_.shape[0] || y>=config_.shape[1] || z>=config_.shape[2]) throw std::out_of_range("grid coordinates outside world");
    return (x*config_.shape[1]+y)*config_.shape[2]+z;
}
std::size_t Universe::index(Vec3 p) const {
    std::array<std::size_t,3> c{};
    for (std::size_t a=0; a<3; ++a) {
        if (!std::isfinite(p[a]) || p[a]<0 || p[a]>config_.extent[a]) throw std::out_of_range("position outside world");
        c[a]=p[a]==config_.extent[a] ? config_.shape[a]-1 : std::min(config_.shape[a]-1,static_cast<std::size_t>(std::floor(p[a]/spacing_[a])));
    }
    return index(c[0],c[1],c[2]);
}
Vec3 Universe::center(std::size_t i) const {
    if (i>=grids_.size()) throw std::out_of_range("voxel index");
    const auto z=i%config_.shape[2], y=(i/config_.shape[2])%config_.shape[1], x=i/(config_.shape[2]*config_.shape[1]);
    return {(static_cast<double>(x)+.5)*spacing_.x,(static_cast<double>(y)+.5)*spacing_.y,(static_cast<double>(z)+.5)*spacing_.z};
}
double Universe::concentration(std::size_t i, SpeciesId s) const { return grids_.at(i).amounts.at(s)/volume_; }
std::size_t Universe::agent_index(AgentId id) const {
    const auto it=std::lower_bound(agents_.begin(),agents_.end(),id,[](const Agent& a, AgentId key) { return a.id<key; });
    if (it==agents_.end() || it->id!=id) throw std::out_of_range("unknown agent ID "+std::to_string(id));
    return static_cast<std::size_t>(it-agents_.begin());
}
const Agent& Universe::agent(AgentId id) const { return agents_[agent_index(id)]; }
std::optional<double> Universe::reservoir(std::size_t i, SpeciesId s) const {
    if (grids_[i].fixed_concentrations[s]) return grids_[i].fixed_concentrations[s];
    if (grids_[i].bulk) return config_.bulk_concentrations[s];
    return {};
}
void Universe::set_concentration(std::size_t i, SpeciesId s, double c) {
    require_nonnegative(c, "field concentration");
    auto& amount=grids_.at(i).amounts.at(s);
    if (reservoir(i,s)) throw std::invalid_argument("use set_reservoir to change a fixed concentration");
    const double n=c*volume_; require_nonnegative(n,"field amount");
    const double external=diagnostics_.external[s]+(n-amount);
    if (!std::isfinite(external)) throw std::invalid_argument("external adjustment overflows ledger");
    diagnostics_.external[s]=external; amount=n;
}
void Universe::set_material(std::size_t i, Material m) {
    if (m!=Material::Boundary && m!=Material::Biofilm) throw std::invalid_argument("unknown material");
    if (config_.layered) throw std::invalid_argument("layered material is derived from particle heights");
    grids_.at(i).material=m;
}
void Universe::set_reservoir(std::size_t i, SpeciesId s, std::optional<double> c) {
    auto& g=grids_.at(i); (void)g.amounts.at(s);
    if (c) { require_nonnegative(*c,"reservoir concentration"); require_nonnegative(*c*volume_,"reservoir amount"); }
    const auto effective=c ? c : (g.bulk ? std::optional<double>{config_.bulk_concentrations[s]} : std::nullopt);
    if (effective) {
        const double n=*effective*volume_, adjustment=diagnostics_.reservoir[s]+(n-g.amounts[s]);
        if (!std::isfinite(adjustment)) throw std::invalid_argument("reservoir adjustment overflows ledger");
        diagnostics_.reservoir[s]=adjustment; g.amounts[s]=n;
    }
    g.fixed_concentrations[s]=c;
}
static double internal_amount(const Agent& a, SpeciesId s) {
    if (!a.qs) return 0;
    const auto& q=*a.qs; const auto& n=a.amounts;
    if (s==q.signal1) return n[A1]+n[RA1]+2*n[C1];
    if (s==q.signal2) return n[A2]+n[RA2]+2*n[C2];
    if (q.inhibitors[0] && s==*q.inhibitors[0]) return n[Q1]+n[RQ1];
    if (q.inhibitors[1] && s==*q.inhibitors[1]) return n[Q2];
    if (q.inhibitors[2] && s==*q.inhibitors[2]) return n[Q3];
    return 0;
}
AgentId Universe::add_impl(const AgentSpec& p) {
    if (agents_.size()>=config_.max_agents || next_id_==std::numeric_limits<AgentId>::max()) throw std::runtime_error("agent or identifier limit reached");
    if (p.kind!=AgentKind::QSBacterium && p.kind!=AgentKind::ChemotacticBacterium && p.kind!=AgentKind::EPS) throw std::invalid_argument("unknown agent kind");
    require_nonnegative(p.biomass,"initial biomass"); require_nonnegative(p.eps,"initial EPS");
    require_positive(p.biomass_density,"biomass density"); require_positive(p.eps_density,"EPS density");
    require_positive(p.division_mass,"division mass"); require_nonnegative(p.death_mass,"death mass"); require_positive(p.shedding_mass,"shedding mass");
    if (p.death_mass>=p.division_mass/2) throw std::invalid_argument("death mass must be below daughter division mass");
    if (p.qs) p.qs->validate(species_);
    if (p.growth) p.growth->validate(species_);
    if (p.motility) p.motility->validate(species_);
    if (p.kind==AgentKind::EPS) {
        if (p.biomass!=0 || p.eps<=0 || p.qs || p.growth || p.motility) throw std::invalid_argument("EPS particles require positive EPS only and no networks");
    } else require_positive(p.biomass,"cell biomass");
    Agent a;
    a.id=next_id_; a.kind=p.kind; a.position=p.position;
    a.biomass=p.biomass; a.eps=p.eps; a.biomass_density=p.biomass_density; a.eps_density=p.eps_density;
    a.division_mass=p.division_mass; a.death_mass=p.death_mass; a.shedding_mass=p.shedding_mass;
    a.division=p.division; a.shedding=p.shedding; a.alive=p.alive && p.kind!=AgentKind::EPS;
    a.qs=p.qs; a.growth=p.growth; a.motility=p.motility;
    require_positive(a.total_volume(),"particle volume");
    for (std::size_t k=0; k<3; ++k)
        if (!std::isfinite(a.position[k]) || a.position[k]<a.radius() || a.position[k]>config_.extent[k]-a.radius() || 2*a.radius()>=config_.extent[k])
            throw std::invalid_argument("initial particle surface must fit inside radius-aware walls");
    for (std::size_t k=0; k<chemical_count; ++k) {
        require_nonnegative(p.initial_concentrations[k],"initial intracellular concentration");
        if (!p.qs && p.initial_concentrations[k]!=0) throw std::invalid_argument("intracellular chemicals require QS parameters");
        a.amounts[k]=p.initial_concentrations[k]*a.volume(); require_nonnegative(a.amounts[k],"intracellular amount");
    }
    if (p.qs) {
        if ((!p.qs->inhibitors[0] && (a.amounts[Q1]!=0 || a.amounts[RQ1]!=0)) ||
            (!p.qs->inhibitors[1] && a.amounts[Q2]!=0) || (!p.qs->inhibitors[2] && a.amounts[Q3]!=0))
            throw std::invalid_argument("initial inhibitor requires an enabled species mapping");
        a.activated=a.alive && a.concentration(C1)>=p.qs->activation_threshold;
    }
    if (!std::isfinite(p.initial_methylation)) throw std::invalid_argument("initial methylation must be finite");
    a.random.seed(stream_seed(config_.seed,a.id));
    a.motion.heading=p.initial_heading ? normalized(*p.initial_heading) : isotropic_direction(a.random);
    a.motion.methylation=p.initial_methylation; a.motion.running=p.initial_running;
    std::exponential_distribution<double> exponential(1);
    do { a.motion.hazard_remaining=exponential(a.random); } while (a.motion.hazard_remaining==0);
    agents_.push_back(a); ++next_id_;
    for (std::size_t s=0; s<species_.size(); ++s) diagnostics_.external[s]+=internal_amount(a,s);
    diagnostics_.external_particle_mass+=a.biomass+a.eps;
    return a.id;
}
AgentId Universe::add_agent(const AgentSpec& p) {
    Universe next=*this;
    const auto id=next.add_impl(p); next.rebuild_membership(); next.update_layers(); next.validate_state();
    *this=std::move(next); return id;
}
void Universe::remove_agent(AgentId id) {
    Universe next=*this;
    const auto i=next.agent_index(id); const auto& a=next.agents_[i];
    for (std::size_t s=0; s<species_.size(); ++s) next.diagnostics_.external[s]-=internal_amount(a,s);
    next.diagnostics_.external_particle_mass-=a.biomass+a.eps;
    next.agents_.erase(next.agents_.begin()+static_cast<std::ptrdiff_t>(i));
    next.rebuild_membership(); next.update_layers(); next.validate_state(); *this=std::move(next);
}
void Universe::rebuild_membership() {
    for (auto& g : grids_) g.members.clear();
    for (const auto& a : agents_) grids_[index(a.position)].members.push_back(a.id);
}
void Universe::update_layers() {
    diagnostics_.biomass_height=0;
    for (const auto& a : agents_) diagnostics_.biomass_height=std::max(diagnostics_.biomass_height,a.position.z+a.radius());
    if (!config_.layered) return;
    for (std::size_t i=0; i<grids_.size(); ++i) {
        auto& g=grids_[i]; const auto z=center(i).z;
        g.material=z<=diagnostics_.biomass_height ? Material::Biofilm : Material::Boundary;
        g.bulk=z>diagnostics_.biomass_height+config_.boundary_layer;
        for (std::size_t s=0; s<species_.size(); ++s) if (const auto c=reservoir(i,s)) {
            const double n=*c*volume_; require_nonnegative(n,"bulk reservoir amount");
            diagnostics_.reservoir[s]+=n-g.amounts[s]; g.amounts[s]=n;
        }
    }
}
double Universe::chemical_total(SpeciesId s) const {
    (void)species_.at(s); double total=0;
    for (const auto& g : grids_) total+=g.amounts[s];
    for (const auto& a : agents_) total+=internal_amount(a,s);
    return total;
}
double Universe::particle_mass() const { double total=0; for (const auto& a : agents_) total+=a.biomass+a.eps; return total; }
double Universe::expected_chemical_total(SpeciesId s) const {
    return diagnostics_.external.at(s)+diagnostics_.reservoir.at(s)+diagnostics_.reaction_source.at(s)-diagnostics_.reaction_sink.at(s)-diagnostics_.decay.at(s);
}
void Universe::validate_state() const {
    std::vector<unsigned> membership(agents_.size());
    for (std::size_t i=0; i<grids_.size(); ++i) {
        for (std::size_t s=0; s<species_.size(); ++s) {
            require_nonnegative(grids_[i].amounts[s],"voxel "+std::to_string(i)+" "+species_[s].name);
            if (const auto c=reservoir(i,s); c && grids_[i].amounts[s]!=*c*volume_) throw std::runtime_error("reservoir concentration changed");
        }
        for (auto id : grids_[i].members) {
            const auto j=agent_index(id); ++membership[j];
            if (index(agents_[j].position)!=i) throw std::runtime_error("stale grid membership");
        }
    }
    for (std::size_t i=0; i<agents_.size(); ++i) {
        const auto& a=agents_[i];
        if (membership[i]!=1) throw std::runtime_error("agent must belong to exactly one voxel");
        require_nonnegative(a.biomass,"agent biomass"); require_nonnegative(a.eps,"agent EPS");
        require_positive(a.total_volume(),"agent volume");
        for (auto n : a.amounts) require_nonnegative(n,"intracellular amount");
        for (std::size_t k=0; k<3; ++k)
            if (!std::isfinite(a.position[k]) || a.position[k]<a.radius()-config_.mechanical_tolerance || a.position[k]>config_.extent[k]-a.radius()+config_.mechanical_tolerance)
                throw std::runtime_error("agent surface outside walls");
        if (a.motility && (!std::isfinite(a.motion.methylation) || std::abs(norm(a.motion.heading)-1)>1e-12 || !(a.motion.hazard_remaining>0)))
            throw std::runtime_error("invalid motility state");
    }
    for (const auto* v : {&diagnostics_.external,&diagnostics_.reservoir,&diagnostics_.decay,&diagnostics_.reaction_source,&diagnostics_.reaction_sink})
        for (auto x : *v) if (!std::isfinite(x)) throw std::runtime_error("nonfinite quantity ledger");
    for (auto x : {diagnostics_.external_particle_mass,diagnostics_.maintenance_loss,diagnostics_.yield_loss})
        if (!std::isfinite(x)) throw std::runtime_error("nonfinite growth ledger");
}
void Universe::movement(double dt) {
    auto next=agents_;
    detail::parallel_for(agents_.size(),config_.workers,[&](std::size_t i) {
        if (agents_[i].alive && agents_[i].motility)
            advance_motility(next[i],concentration(index(agents_[i].position),agents_[i].motility->attractant),dt,config_.extent);
    });
    agents_.swap(next);
}
void Universe::advance(double dt) {
    transport(dt/2); biology(dt); movement(dt); lifecycle(); relax(); rebuild_membership(); update_layers(); transport(dt/2);
    validate_state(); time_+=dt; ++step_;
}
void Universe::step(double dt) {
    require_positive(dt,"timestep");
    if (dt/2==0 || !std::isfinite(time_+dt) || time_+dt==time_ || step_==std::numeric_limits<std::uint64_t>::max()) throw std::invalid_argument("timestep cannot advance representable time/index");
    Universe next=*this;
    try { next.advance(dt); }
    catch (const std::exception& e) { throw std::runtime_error("step "+std::to_string(step_)+" at t="+std::to_string(time_)+": "+e.what()); }
    *this=std::move(next);
}
void Universe::run_until(double end_time, double dt) {
    require_nonnegative(end_time,"end time"); require_positive(dt,"timestep");
    if (end_time<time_) throw std::invalid_argument("end time precedes simulation time");
    while (time_<end_time) {
        const bool final=dt>=end_time-time_;
        step(std::min(dt,end_time-time_));
        if (final) time_=end_time;
    }
}
} // namespace BNSim
