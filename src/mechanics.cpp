// Deterministic lifecycle and spatial contact relaxation, GPL-2.0.
#include "bnsim/simulation.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace BNSim {
static Vec3 pair_direction(AgentId a, AgentId b) {
    std::mt19937_64 stream(stream_seed(a,b));
    return isotropic_direction(stream); // deterministic geometry, independent of motility RNG
}
static void project_walls(Agent& a, Vec3 extent) {
    const double r=a.radius(); require_positive(r,"mechanical radius");
    for (std::size_t k=0; k<3; ++k) {
        if (2*r>=extent[k]) throw std::runtime_error("particle diameter cannot fit inside world");
        if (!std::isfinite(a.position[k])) throw std::runtime_error("nonfinite mechanical position");
        a.position[k]=std::clamp(a.position[k],r,extent[k]-r);
    }
}
void Universe::lifecycle() {
    auto new_id=[&]() {
        if (agents_.size()>=config_.max_agents || next_id_==std::numeric_limits<AgentId>::max()) throw std::runtime_error("lifecycle agent/ID limit reached");
        return next_id_++;
    };
    // Children are appended in parent ID order. Oversized daughters divide again
    // within this stage, bounded by max_agents, without biochemical randomness.
    for (std::size_t i=0; i<agents_.size(); ++i) {
        if (!agents_[i].alive || agents_[i].kind==AgentKind::EPS) continue;
        if (agents_[i].biomass<=agents_[i].death_mass) {
            agents_[i].alive=false; agents_[i].activated=false; ++diagnostics_.deaths; continue;
        }
        while (agents_[i].division && agents_[i].biomass>=agents_[i].division_mass) {
            Agent child=agents_[i]; child.id=new_id();
            agents_[i].biomass*=.5; agents_[i].eps*=.5;
            for (auto& n : agents_[i].amounts) n*=.5;
            child.biomass=agents_[i].biomass; child.eps=agents_[i].eps; child.amounts=agents_[i].amounts;
            child.random.seed(stream_seed(config_.seed,child.id));
            std::exponential_distribution<double> exponential(1);
            do { child.motion.hazard_remaining=exponential(child.random); } while (child.motion.hazard_remaining==0);
            child.motion.run_seconds=child.motion.tumble_seconds=0; child.motion.switches=child.motion.reorientations=0;
            const auto displacement=pair_direction(agents_[i].id,child.id)*child.radius();
            child.position=child.position+displacement;
            agents_[i].position=agents_[i].position-displacement;
            agents_.push_back(std::move(child)); ++diagnostics_.births;
        }
        if (agents_[i].shedding && agents_[i].eps>=agents_[i].shedding_mass) {
            Agent particle; particle.id=new_id(); particle.kind=AgentKind::EPS; particle.alive=false;
            particle.division=particle.shedding=false; particle.eps=agents_[i].eps; particle.eps_density=agents_[i].eps_density;
            agents_[i].eps=0;
            particle.position=agents_[i].position+pair_direction(agents_[i].id,particle.id)*(agents_[i].radius()+particle.radius());
            particle.random.seed(stream_seed(config_.seed,particle.id));
            agents_.push_back(std::move(particle)); ++diagnostics_.shed_particles;
        }
    }
    for (auto& a : agents_) project_walls(a,config_.extent);
}
void Universe::relax() {
    diagnostics_.mechanical_iterations=0; diagnostics_.max_overlap=0;
    if (!config_.mechanics || agents_.empty()) return;
    double largest=0; for (const auto& a : agents_) largest=std::max(largest,a.radius());
    for (std::size_t iteration=0; iteration<config_.mechanical_iterations; ++iteration) {
        rebuild_membership(); double overlap_max=0;
        // Bins refer to positions at the beginning of this sweep. Corrections
        // can create new contacts; the following sweep discovers those.
        std::vector<Vec3> positions; positions.reserve(agents_.size());
        for (const auto& a : agents_) positions.push_back(a.position);
        for (std::size_t i=0; i<agents_.size(); ++i) {
            const double reach=agents_[i].radius()+largest;
            std::array<std::size_t,3> low{},high{};
            for (std::size_t k=0; k<3; ++k) {
                low[k]=static_cast<std::size_t>(std::floor(std::max(0.,positions[i][k]-reach)/spacing_[k]));
                high[k]=std::min(config_.shape[k]-1,static_cast<std::size_t>(std::floor(std::min(config_.extent[k],positions[i][k]+reach)/spacing_[k])));
            }
            for (std::size_t x=low[0]; x<=high[0]; ++x)
                for (std::size_t y=low[1]; y<=high[1]; ++y)
                    for (std::size_t z=low[2]; z<=high[2]; ++z)
                        for (const auto id : grids_[index(x,y,z)].members) if (id>agents_[i].id) {
                            auto& a=agents_[i]; auto& b=agents_[agent_index(id)];
                            const auto delta=b.position-a.position; const double distance=norm(delta);
                            const double overlap=a.radius()+b.radius()-distance;
                            overlap_max=std::max(overlap_max,overlap);
                            if (overlap>config_.mechanical_tolerance) {
                                const auto direction=distance>1e-14 ? delta*(1/distance) : pair_direction(a.id,b.id);
                                a.position=a.position-direction*(overlap/2); b.position=b.position+direction*(overlap/2);
                                project_walls(a,config_.extent); project_walls(b,config_.extent);
                            }
                        }
        }
        diagnostics_.mechanical_iterations=iteration+1; diagnostics_.max_overlap=overlap_max;
        if (overlap_max<=config_.mechanical_tolerance) { rebuild_membership(); return; }
    }
    throw std::runtime_error("mechanical relaxation failed after "+std::to_string(config_.mechanical_iterations)+" sweeps; maximum overlap="+std::to_string(diagnostics_.max_overlap)+" um");
}
} // namespace BNSim
