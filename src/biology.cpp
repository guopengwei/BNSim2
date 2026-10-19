// Coupled voxel/cell amount integration, GPL-2.0.
#include "bnsim/simulation.hpp"
#include "numerics.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace BNSim {
void Universe::biology(double dt) {
    constexpr std::size_t stride=2+chemical_count;
    const auto m=species_.size();
    struct Result { std::vector<std::size_t> cells; std::vector<double> y; detail::OdeStats stats; };
    std::vector<Result> results(grids_.size());
    detail::parallel_for(grids_.size(),config_.workers,[&](std::size_t voxel) {
        auto& result=results[voxel];
        for (const auto id : grids_[voxel].members) {
            const auto i=agent_index(id); const auto& a=agents_[i];
            if (a.alive && (a.qs || a.growth)) result.cells.push_back(i);
        }
        if (result.cells.empty()) return;
        const auto base=m+stride*result.cells.size();
        const auto source=base, sink=base+m, inflow=base+2*m, outflow=base+3*m, losses=base+4*m;
        result.y.resize(losses+2);
        std::copy(grids_[voxel].amounts.begin(),grids_[voxel].amounts.end(),result.y.begin());
        std::vector<double> atol(result.y.size(),config_.integration.molecular_absolute_tolerance);
        for (std::size_t s=0; s<m; ++s) if (species_[s].quantity==Quantity::Mass)
            for (const auto offset : {std::size_t{0},source,sink,inflow,outflow}) atol[offset+s]=config_.integration.mass_absolute_tolerance;
        atol[losses]=atol[losses+1]=config_.integration.mass_absolute_tolerance;
        for (std::size_t j=0; j<result.cells.size(); ++j) {
            const auto& a=agents_[result.cells[j]]; const auto start=m+j*stride;
            result.y[start]=a.biomass; result.y[start+1]=a.eps;
            atol[start]=atol[start+1]=config_.integration.mass_absolute_tolerance;
            std::copy(a.amounts.begin(),a.amounts.end(),result.y.begin()+static_cast<std::ptrdiff_t>(start+2));
        }
        const auto rhs=[&](const std::vector<double>& y, std::vector<double>& d) {
            for (std::size_t j=0; j<result.cells.size(); ++j) {
                const auto& a=agents_[result.cells[j]]; const auto start=m+j*stride;
                const double volume=y[start]/a.biomass_density;
                if (!(volume>0)) { d[start]=std::numeric_limits<double>::quiet_NaN(); return; }
                Chemicals n{};
                std::copy_n(y.begin()+static_cast<std::ptrdiff_t>(start+2),chemical_count,n.begin());
                if (a.qs) {
                    const auto& p=*a.qs;
                    const auto rates=qs_reactions(p,n,volume);
                    for (std::size_t k=0; k<chemical_count; ++k) d[start+2+k]+=rates.derivative[k];
                    d[source+p.signal1]+=rates.signal_source[0]; d[source+p.signal2]+=rates.signal_source[1];
                    d[sink+p.signal1]+=rates.signal_sink[0]; d[sink+p.signal2]+=rates.signal_sink[1];
                    const std::array<std::optional<SpeciesId>,5> mapped{p.signal1,p.signal2,p.inhibitors[0],p.inhibitors[1],p.inhibitors[2]};
                    constexpr std::array<Chemical,5> internal{A1,A2,Q1,Q2,Q3};
                    for (std::size_t k=0; k<5; ++k) if (mapped[k]) {
                        const auto s=*mapped[k];
                        const double flux=p.exchange[k]*(n[internal[k]]-y[s]/volume_*volume);
                        d[start+2+internal[k]]-=flux; d[s]+=flux;
                    }
                }
                if (a.growth) {
                    const auto& p=*a.growth;
                    const auto rates=growth_rates(p,y[start],y[p.substrate]/volume_,n[C1]/volume);
                    d[start]+=rates.biomass; d[start+1]+=rates.eps; d[p.substrate]-=rates.substrate;
                    d[sink+p.substrate]+=rates.substrate;
                    d[losses]+=rates.maintenance_loss; d[losses+1]+=rates.yield_loss;
                }
            }
            for (std::size_t s=0; s<m; ++s) if (reservoir(voxel,s)) {
                if (d[s]>0) d[outflow+s]=d[s]; else d[inflow+s]=-d[s];
                d[s]=0;
            }
        };
        try { result.stats=detail::integrate(result.y,dt,atol,config_.integration,rhs); }
        catch (const std::exception& e) { throw std::runtime_error("voxel "+std::to_string(voxel)+" chemistry: "+e.what()); }
    });
    // Deterministic commit and reduction; no cell or field was modified by workers.
    for (std::size_t voxel=0; voxel<results.size(); ++voxel) {
        const auto& result=results[voxel]; if (result.cells.empty()) continue;
        const auto base=m+stride*result.cells.size(); const auto& y=result.y;
        for (std::size_t s=0; s<m; ++s) {
            grids_[voxel].amounts[s]=y[s];
            diagnostics_.reaction_source[s]+=y[base+s]; diagnostics_.reaction_sink[s]+=y[base+m+s];
            diagnostics_.reservoir[s]+=y[base+2*m+s]-y[base+3*m+s];
        }
        diagnostics_.maintenance_loss+=y[base+4*m]; diagnostics_.yield_loss+=y[base+4*m+1];
        diagnostics_.ode_accepted+=result.stats.accepted; diagnostics_.ode_rejected+=result.stats.rejected;
        diagnostics_.ode_startup+=result.stats.startup;
        for (std::size_t j=0; j<result.cells.size(); ++j) {
            auto& a=agents_[result.cells[j]]; const auto start=m+j*stride;
            a.biomass=y[start]; a.eps=y[start+1];
            std::copy_n(y.begin()+static_cast<std::ptrdiff_t>(start+2),chemical_count,a.amounts.begin());
            a.activated=a.qs && a.concentration(C1)>=a.qs->activation_threshold;
        }
    }
}
} // namespace BNSim
