// Conservative finite-volume transport, GPL-2.0.
#include "bnsim/simulation.hpp"
#include "numerics.hpp"
#include <algorithm>
#include <stdexcept>

namespace BNSim {
void Universe::decay(double dt) {
    for (std::size_t s=0; s<species_.size(); ++s)
        for (std::size_t i=0; i<grids_.size(); ++i) {
            auto& g=grids_[i];
            const double k=g.material==Material::Biofilm ? species_[s].biofilm_decay : species_[s].decay;
            const double loss=g.amounts[s]*(-std::expm1(-k*dt));
            diagnostics_.decay[s]+=loss;
            if (reservoir(i,s)) diagnostics_.reservoir[s]+=loss;
            else g.amounts[s]*=std::exp(-k*dt);
        }
}
void Universe::transport(double dt) {
    decay(dt/2);
    if (config_.diffusion && !faces_.empty()) for (std::size_t s=0; s<species_.size(); ++s) {
        std::vector<double> coefficient(faces_.size()), outgoing(grids_.size()), flux(faces_.size()), delta(grids_.size());
        for (std::size_t f=0; f<faces_.size(); ++f) {
            const auto& face=faces_[f];
            const auto diffusion=[&](std::size_t i) { return grids_[i].material==Material::Biofilm ? species_[s].biofilm_diffusion : species_[s].diffusion; };
            const double a=diffusion(face.a), b=diffusion(face.b);
            // Stable harmonic mean, including impermeable interfaces.
            const double small=std::min(a,b), large=std::max(a,b);
            const double harmonic=small==0 ? 0 : small*(2/(1+small/large));
            coefficient[f]=harmonic*face.inverse_spacing_squared;
            require_nonnegative(coefficient[f],"transport face coefficient");
            outgoing[face.a]+=coefficient[f]; outgoing[face.b]+=coefficient[f];
        }
        const double maximum=*std::max_element(outgoing.begin(),outgoing.end());
        require_nonnegative(maximum,"multidimensional diffusion outflow rate");
        const double required=std::max(1.,std::ceil(dt*maximum/.9));
        if (!std::isfinite(required) || required>static_cast<double>(config_.max_transport_substeps))
            throw std::runtime_error("diffusion positivity substeps exceed configured limit");
        const auto steps=static_cast<std::size_t>(required);
        const double h=dt/static_cast<double>(steps);
        for (std::size_t sub=0; sub<steps; ++sub) {
            // Field is immutable while independent face fluxes are computed.
            detail::parallel_for(faces_.size(),config_.workers,[&](std::size_t f) {
                const auto& face=faces_[f];
                flux[f]=h*coefficient[f]*(grids_[face.a].amounts[s]-grids_[face.b].amounts[s]);
            });
            std::fill(delta.begin(),delta.end(),0);
            for (std::size_t f=0; f<faces_.size(); ++f) {
                delta[faces_[f].a]-=flux[f]; delta[faces_[f].b]+=flux[f];
            }
            for (std::size_t i=0; i<grids_.size(); ++i) {
                if (reservoir(i,s)) diagnostics_.reservoir[s]-=delta[i];
                else {
                    const double next=grids_[i].amounts[s]+delta[i];
                    require_nonnegative(next,"transport amount at voxel "+std::to_string(i));
                    grids_[i].amounts[s]=next;
                }
            }
            ++diagnostics_.transport_substeps;
        }
    }
    decay(dt/2);
}
} // namespace BNSim
