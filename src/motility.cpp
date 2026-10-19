// Reduced single-motor chemotaxis, GPL-2.0. See docs/model.md for limitations.
#include "bnsim/model.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace BNSim {
static double uniform_open(std::mt19937_64& rng) {
    return (static_cast<double>(rng() >> 12) + .5) / 4503599627370496.;
}
static double exponential_hazard(std::mt19937_64& rng) { return -std::log(uniform_open(rng)); }
Vec3 isotropic_direction(std::mt19937_64& rng) {
    const double z = 2*uniform_open(rng)-1, phi = 2*pi*uniform_open(rng);
    const double r = std::sqrt(std::max(0., 1-z*z));
    return {r*std::cos(phi), r*std::sin(phi), z};
}
static Vec3 tangent(Vec3 h) {
    return normalized(cross(h, std::abs(h.z) < .9 ? Vec3{0,0,1} : Vec3{0,1,0}));
}
void rotational_increment(Vec3& heading, double diffusion, double dt, std::mt19937_64& rng) {
    require_nonnegative(diffusion, "rotational diffusion"); require_nonnegative(dt, "rotation dt");
    if (diffusion == 0 || dt == 0) return;
    const Vec3 u = tangent(heading), v = cross(heading,u);
    std::normal_distribution<double> gaussian(0, std::sqrt(2*diffusion*dt));
    const double a = gaussian(rng), b = gaussian(rng), angle = std::hypot(a,b);
    if (angle > 0) heading = normalized(heading*std::cos(angle)+(u*a+v*b)*(std::sin(angle)/angle));
}
double tumble_angle(const MotilityParameters& p, std::mt19937_64& rng) {
    require_positive(p.gamma_shape,"tumble gamma shape"); require_positive(p.gamma_scale_degrees,"tumble gamma scale");
    std::gamma_distribution<double> gamma(p.gamma_shape, p.gamma_scale_degrees);
    for (int i=0; i<10000; ++i) {
        const double degrees = gamma(rng);
        if (std::isfinite(degrees) && degrees >= 0 && degrees <= 180) return degrees*pi/180;
    }
    throw std::runtime_error("tumble gamma rejection limit reached; check shape and scale");
}
double receptor_activity(const MotilityParameters& p, double m, double ligand) {
    require_nonnegative(ligand, "attractant concentration");
    if (!std::isfinite(m)) throw std::invalid_argument("methylation must be finite");
    // Log differences avoid overflow for very high concentrations.
    const auto log_ratio=[&](double binding) {
        return ligand>binding ? std::log(ligand)-std::log(binding)+std::log1p(binding/ligand) : std::log1p(ligand/binding);
    };
    const double binding = log_ratio(p.active_binding)-log_ratio(p.inactive_binding);
    const double energy = p.receptors*(p.alpha*(p.m0-m)-binding);
    if (energy >= 0) { const auto e=std::exp(-energy); return e/(1+e); }
    return 1/(1+std::exp(energy));
}
std::array<double,2> motor_rates(const MotilityParameters& p, double activity) {
    if (!std::isfinite(activity) || activity < 0 || activity > 1) throw std::invalid_argument("activity must be in [0,1]");
    const double y = p.ya*activity, fraction=y==0 ? 0 : 1/(1+p.kd/y), energy = p.g*(fraction-.5)/2;
    return {p.motor_frequency*std::exp(energy), p.motor_frequency*std::exp(-energy)};
}
void reflect_position(Vec3& pos, Vec3& heading, Vec3 extent, double radius) {
    require_nonnegative(radius, "wall radius");
    for (std::size_t axis=0; axis<3; ++axis) {
        const double width = extent[axis]-2*radius;
        require_positive(width, "available wall interval");
        if (!std::isfinite(pos[axis])) throw std::runtime_error("nonfinite position during reflection");
        double q = std::fmod(pos[axis]-radius, 2*width);
        if (q < 0) q += 2*width;
        if (q > width) { pos[axis] = radius+2*width-q; heading[axis] = -heading[axis]; }
        else pos[axis] = radius+q;
        if ((q == 0 && heading[axis] < 0) || (q == width && heading[axis] > 0)) heading[axis] = -heading[axis];
    }
}
static void adapt(MotilityState& state, const MotilityParameters& p, double ligand, double dt) {
    const auto f = [&](double m) { const auto a=receptor_activity(p,m,ligand); return p.kr*(1-a)-p.kb*a; };
    const auto m=state.methylation, a=f(m), b=f(m+dt*a/2), c=f(m+dt*b/2), d=f(m+dt*c);
    state.methylation += dt*(a+2*b+2*c+d)/6;
}
void advance_motility(Agent& cell, double ligand, double dt, Vec3 extent) {
    if (!cell.alive || !cell.motility) return;
    const auto& p=*cell.motility;
    auto& s=cell.motion;
    double elapsed=0;
    std::size_t attempts=0;
    while (elapsed < dt) {
        if (++attempts > 10000000) throw std::runtime_error("motility interval/event limit reached");
        const double slice=std::min(p.max_step,dt-elapsed);
        const auto rates=motor_rates(p,receptor_activity(p,s.methylation,ligand));
        const double rate=rates[s.running ? 0 : 1];
        const double wait=rate > 0 ? s.hazard_remaining/rate : std::numeric_limits<double>::infinity();
        const bool event=wait <= slice;
        const double h=event ? wait : slice;
        if (!(h > 0) || elapsed+h == elapsed) throw std::runtime_error("motility event time underflow");
        if (s.running) {
            const Vec3 old=s.heading;
            rotational_increment(s.heading,p.rotational_diffusion,h,cell.random);
            const Vec3 midpoint=norm(old+s.heading) > 1e-10 ? normalized(old+s.heading) : old;
            cell.position=cell.position+midpoint*(p.speed*h);
            reflect_position(cell.position,s.heading,extent,cell.radius());
            s.run_seconds+=h;
        } else s.tumble_seconds+=h;
        adapt(s,p,ligand,h);
        elapsed = h == dt-elapsed ? dt : elapsed+h;
        s.hazard_remaining-=rate*h;
        if (event) {
            s.running=!s.running; ++s.switches;
            s.hazard_remaining=exponential_hazard(cell.random);
            if (s.running) {
                const double angle=tumble_angle(p,cell.random), azimuth=2*pi*uniform_open(cell.random);
                const auto u=tangent(s.heading), v=cross(s.heading,u);
                s.heading=normalized(s.heading*std::cos(angle)+(u*std::cos(azimuth)+v*std::sin(azimuth))*std::sin(angle));
                ++s.reorientations;
            }
        }
    }
}
} // namespace BNSim
