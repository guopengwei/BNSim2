// Checked, reproducible chemotaxis and QS/biofilm examples, GPL-2.0.
#include "bnsim/io.hpp"
#include <algorithm>
#include <charconv>
#include <iostream>
#include <set>
#include <stdexcept>

using namespace BNSim;
namespace {
std::uint64_t integer(const std::string& s) {
    std::uint64_t value=0;
    const auto result=std::from_chars(s.data(),s.data()+s.size(),value);
    if (result.ec!=std::errc{} || result.ptr!=s.data()+s.size()) throw std::invalid_argument("invalid unsigned integer: "+s);
    return value;
}
double real(const std::string& s) {
    std::size_t end=0; const double value=std::stod(s,&end);
    if (end!=s.size() || !std::isfinite(value)) throw std::invalid_argument("invalid finite number: "+s);
    return value;
}
}
int main(int argc, char** argv) {
    try {
        SimulationConfig config;
        double duration=1, dt=.1; std::size_t population=20; std::filesystem::path output_path="runs";
        std::set<std::string> seen;
        for (int i=1; i<argc; ++i) {
            const std::string arg=argv[i];
            if (arg=="--help") {
                std::cout << "Options: --seed N --threads N --duration SECONDS --dt SECONDS --population N --output DIRECTORY\n"; return 0;
            }
            if (!seen.insert(arg).second) throw std::invalid_argument("duplicate option: "+arg);
            if (i+1==argc) throw std::invalid_argument("missing value for "+arg);
            const std::string value=argv[++i];
            if (arg=="--seed") config.seed=integer(value);
            else if (arg=="--threads") { const auto n=integer(value); if (n>4096) throw std::invalid_argument("too many threads"); config.workers=static_cast<std::size_t>(n); }
            else if (arg=="--duration") duration=real(value);
            else if (arg=="--dt") dt=real(value);
            else if (arg=="--population") { const auto n=integer(value); if (n>config.max_agents) throw std::invalid_argument("population exceeds agent limit"); population=static_cast<std::size_t>(n); }
            else if (arg=="--output") output_path=value;
            else throw std::invalid_argument("unknown option: "+arg);
        }
        require_nonnegative(duration,"duration"); require_positive(dt,"dt");
#ifdef BNSIM_QS_EXAMPLE
        const std::string scenario="qs_biofilm";
        config.extent={60,60,60}; config.shape={6,6,12}; config.layered=true; config.boundary_layer=10;
        config.bulk_concentrations={1,0,0,0,0,0};
        std::vector<Species> species{{"substrate",Quantity::Mass,30,18,0,0,"illustrative transport; provisional"},
            {"AI1",Quantity::Molar,30,18,.01,.01,"illustrative transport; legacy decay; provisional"},
            {"AI2",Quantity::Molar,30,18,.01,.01,"illustrative transport; legacy decay; provisional"},
            {"QSI1",Quantity::Molar,30,18,0,0,"illustrative; provisional"},
            {"QSI2",Quantity::Molar,30,18,0,0,"illustrative; provisional"},
            {"QSI3",Quantity::Molar,30,18,0,0,"illustrative; provisional"}};
#else
        const std::string scenario="chemotaxis";
        config.extent={200,100,50}; config.shape={20,10,5}; config.diffusion=false;
        std::vector<Species> species{{"Aspartate",Quantity::Molar,0,0,0,0,"prescribed static gradient; illustrative"}};
#endif
        Universe universe(config,species);
#ifndef BNSIM_QS_EXAMPLE
        std::mt19937_64 setup(stream_seed(config.seed,123456789));
        std::uniform_real_distribution<double> uniform(.1,.9);
#endif
        for (std::size_t i=0; i<population; ++i) {
            AgentSpec cell;
#ifdef BNSIM_QS_EXAMPLE
            const auto side=static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(population))));
            cell.position={config.extent.x*static_cast<double>(i%side+1)/static_cast<double>(side+1),
                           config.extent.y*static_cast<double>(i/side+1)/static_cast<double>(side+1),1};
            cell.qs=QSParameters{}; cell.qs->signal1=1; cell.qs->signal2=2; cell.qs->inhibitors={3,4,5};
            cell.growth=GrowthParameters{};
#else
            cell.kind=AgentKind::ChemotacticBacterium;
            cell.position={config.extent.x*uniform(setup),config.extent.y*uniform(setup),config.extent.z*uniform(setup)};
            cell.motility=MotilityParameters{}; cell.division=false; cell.shedding=false;
#endif
            universe.add_agent(cell);
        }
        for (std::size_t i=0; i<universe.grids().size(); ++i) {
#ifdef BNSIM_QS_EXAMPLE
            if (!universe.grids()[i].bulk) universe.set_concentration(i,0,1);
#else
            universe.set_concentration(i,0,from_millimolar(.01+.0001*(config.extent.y-universe.center(i).y)));
#endif
        }
        RunWriter writer(output_path,universe,duration,dt,scenario); writer.snapshot(universe);
        while (universe.time()<duration) {
            universe.run_until(std::min(duration,universe.time()+dt),dt); writer.snapshot(universe);
        }
        std::cout << "Completed " << scenario << " at " << universe.time() << " s; " << universe.agents().size()
                  << " agents; output: " << writer.directory() << '\n';
        return 0;
    } catch (const std::exception& e) { std::cerr << "BNSim: " << e.what() << '\n'; return 1; }
}
