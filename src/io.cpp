// BNSim, GPL-2.0. Outputs use full precision and never overwrite a snapshot.
#include "bnsim/io.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifndef BNSIM_CODE_VERSION
#define BNSIM_CODE_VERSION "unknown"
#endif
#ifndef BNSIM_SOURCE_SHA256
#define BNSIM_SOURCE_SHA256 "unknown"
#endif
#ifndef BNSIM_COMPILER
#define BNSIM_COMPILER "unknown"
#endif
#ifndef BNSIM_BUILD_TYPE
#define BNSIM_BUILD_TYPE "unknown"
#endif

namespace BNSim {
static std::string quoted(const std::string& s) {
    std::ostringstream out; out << '"';
    for (const unsigned char c : s) {
        if (c=='"' || c=='\\') out << '\\' << c;
        else if (c<32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<unsigned>(c);
        else out << c;
    }
    out << '"'; return out.str();
}
static std::ofstream output(const std::filesystem::path& path) {
    if (std::filesystem::exists(path)) throw std::runtime_error("refusing to overwrite "+path.string());
    std::ofstream out; out.exceptions(std::ios::badbit|std::ios::failbit);
    out.open(path); out << std::setprecision(17); return out;
}
static std::string suffix(std::uint64_t index) {
    std::ostringstream s; s << std::setw(8) << std::setfill('0') << index; return s.str();
}
static const char* amount_unit(const Species& s) { return s.quantity==Quantity::Mass ? "fg" : "fmol"; }
RunWriter::RunWriter(const std::filesystem::path& parent, const Universe& u, double duration, double dt, const std::string& scenario) {
    require_nonnegative(duration,"output duration"); require_positive(dt,"output timestep");
    if (parent.empty()) throw std::invalid_argument("output directory is empty");
    std::filesystem::create_directories(parent);
    bool made=false;
    for (std::uint64_t i=0; i<1000000; ++i) {
        directory_=parent/("run-"+suffix(i));
        if (std::filesystem::create_directory(directory_)) { made=true; break; }
    }
    if (!made) throw std::runtime_error("could not allocate a unique run directory");
    metadata(u,duration,dt,scenario);
}
void RunWriter::metadata(const Universe& u, double duration, double dt, const std::string& scenario) const {
    const auto& c=u.config(); const auto& o=c.integration;
    auto out=output(directory_/"metadata.json");
    out << "{\n  \"schema\": 1,\n  \"scenario\": " << quoted(scenario)
        << ",\n  \"code_version\": " << quoted(BNSIM_CODE_VERSION) << ",\n  \"source_sha256\": " << quoted(BNSIM_SOURCE_SHA256)
        << ",\n  \"compiler\": " << quoted(BNSIM_COMPILER) << ",\n  \"build_type\": " << quoted(BNSIM_BUILD_TYPE)
        << ",\n  \"seed\": " << c.seed << ",\n  \"workers\": " << c.workers
        << ",\n  \"duration_s\": " << duration << ",\n  \"timestep_s\": " << dt
        << ",\n  \"initial_time_s\": " << u.time() << ",\n  \"initial_population\": " << u.agents().size()
        << ",\n  \"units\": {\"length\":\"um\",\"time\":\"s\",\"mass\":\"fg\",\"chemical_amount\":\"fmol\",\"concentration\":\"species amount / um^3\"},"
        << "\n  \"extent_um\": [" << c.extent.x << ',' << c.extent.y << ',' << c.extent.z << "],"
        << "\n  \"shape\": [" << c.shape[0] << ',' << c.shape[1] << ',' << c.shape[2] << "],"
        << "\n  \"diffusion\": " << (c.diffusion ? "true" : "false") << ",\n  \"mechanics\": " << (c.mechanics ? "true" : "false")
        << ",\n  \"layered\": " << (c.layered ? "true" : "false") << ",\n  \"boundary_layer_um\": " << c.boundary_layer
        << ",\n  \"mechanical_tolerance_um\": " << c.mechanical_tolerance << ",\n  \"mechanical_iterations\": " << c.mechanical_iterations
        << ",\n  \"max_voxels\": " << c.max_voxels << ",\n  \"max_agents\": " << c.max_agents << ",\n  \"max_transport_substeps\": " << c.max_transport_substeps
        << ",\n  \"integration\": {\"rtol\":" << o.relative_tolerance << ",\"molecular_atol_fmol\":" << o.molecular_absolute_tolerance
        << ",\"mass_atol_fg\":" << o.mass_absolute_tolerance << ",\"minimum_step_s\":" << o.minimum_step << ",\"max_attempts\":" << o.max_attempts << "},"
        << "\n  \"bulk_concentrations\": [";
    for (std::size_t i=0; i<c.bulk_concentrations.size(); ++i) out << (i ? "," : "") << c.bulk_concentrations[i];
    out << "],\n  \"species\": [";
    for (std::size_t i=0; i<u.species().size(); ++i) {
        const auto& s=u.species()[i];
        out << (i ? "," : "") << "{\"id\":" << i << ",\"name\":" << quoted(s.name) << ",\"amount_unit\":" << quoted(amount_unit(s))
            << ",\"diffusion_um2_per_s\":" << s.diffusion << ",\"biofilm_diffusion_um2_per_s\":" << s.biofilm_diffusion
            << ",\"decay_per_s\":" << s.decay << ",\"biofilm_decay_per_s\":" << s.biofilm_decay << ",\"provenance\":" << quoted(s.provenance) << '}';
    }
    out << "],\n  \"parameters\":\"parameters.csv\",\n  \"model_status\":\"corrected legacy; provisional; not experimentally calibrated\","
        << "\n  \"boundaries\":\"radius-aware reflecting mechanical walls; no-flux chemicals except recorded reservoirs\","
        << "\n  \"stages\":\"half transport/decay; coupled local biology; motility/lifecycle; mechanics/layers; half transport/decay\","
        << "\n  \"sources\":[\"https://doi.org/10.1038/srep35136\",\"https://doi.org/10.1073/pnas.1113706109\"],"
        << "\n  \"reproducibility\":\"same build/platform; initial snapshots are data exports, not RNG checkpoints\"\n}\n";
    out.close();
    auto p=output(directory_/"parameters.csv");
    p << "agent_id,parameter,value,unit,status,provenance\n";
    const auto row=[&](AgentId id,const std::string& name,double value,const char* unit,const char* source="configured value; default origins in docs/parameters.md",const char* status="provisional") {
        p << id << ',' << name << ',' << value << ',' << unit << ',' << status << ',' << source << '\n';
    };
    for (const auto& a : u.agents()) {
        row(a.id,"kind",static_cast<int>(a.kind),"enum");
        row(a.id,"biomass_density",a.biomass_density,"fg/um^3"); row(a.id,"eps_density",a.eps_density,"fg/um^3");
        row(a.id,"division_mass",a.division_mass,"fg"); row(a.id,"death_mass",a.death_mass,"fg"); row(a.id,"shedding_mass",a.shedding_mass,"fg");
        row(a.id,"division_enabled",a.division,"bool"); row(a.id,"shedding_enabled",a.shedding,"bool");
        if (a.qs) {
            const auto& q=*a.qs;
            row(a.id,"qs.signal1",static_cast<double>(q.signal1),"species_id"); row(a.id,"qs.signal2",static_cast<double>(q.signal2),"species_id");
            for (std::size_t k=0; k<3; ++k) row(a.id,"qs.inhibitor"+std::to_string(k+1),q.inhibitors[k] ? static_cast<double>(*q.inhibitors[k]) : -1,"species_id_or_minus1");
            for (std::size_t k=0; k<5; ++k) {
                row(a.id,"qs.basal"+std::to_string(k),q.basal[k],"fmol/um^3/s");
                row(a.id,"qs.induced"+std::to_string(k),q.induced[k],"fmol/um^3/s");
                row(a.id,"qs.exchange"+std::to_string(k),q.exchange[k],"1/s");
            }
            for (std::size_t k=0; k<7; ++k) row(a.id,"qs.threshold"+std::to_string(k),q.threshold[k],"fmol/um^3");
            for (std::size_t k=0; k<13; ++k) row(a.id,"qs.rate"+std::to_string(k),q.rate[k],(k==1||k==4||k==8||k==11) ? "um^3/fmol/s" : "1/s");
            row(a.id,"qs.inhibitor_binding",q.inhibitor_binding,"um^3/fmol/s"); row(a.id,"qs.inhibitor_unbinding",q.inhibitor_unbinding,"1/s");
            row(a.id,"qs.signal_degradation",q.signal_degradation,"um^3/fmol/s"); row(a.id,"qs.inhibitor_saturation",q.inhibitor_saturation,"um^3/fmol");
            row(a.id,"qs.receptor_suppression",q.receptor_suppression,"fmol/um^3"); row(a.id,"qs.activation_threshold",q.activation_threshold,"fmol/um^3");
            row(a.id,"qs.rsal_feedback",q.rsal_feedback,"bool");
        }
        if (a.growth) {
            const auto& g=*a.growth;
            row(a.id,"growth.substrate",static_cast<double>(g.substrate),"species_id");
            row(a.id,"growth.umax",g.umax,"1/s"); row(a.id,"growth.pmax",g.pmax,"1/s"); row(a.id,"growth.maintenance",g.maintenance,"1/s");
            row(a.id,"growth.substrate_half",g.substrate_half,"fg/um^3"); row(a.id,"growth.qs_half",g.qs_half,"fmol/um^3"); row(a.id,"growth.yield",g.yield,"fg/fg");
        }
        if (a.motility) {
            const auto& m=*a.motility;
            row(a.id,"motility.attractant",static_cast<double>(m.attractant),"species_id");
            row(a.id,"motility.speed",m.speed,"um/s","default from Sneddon 2012 Methods",m.speed==20 ? "literature default" : "user override; provisional");
            row(a.id,"motility.rotational_diffusion",m.rotational_diffusion,"rad^2/s","default from Sneddon 2012 Methods",m.rotational_diffusion==.062 ? "literature default" : "user override; provisional");
            row(a.id,"motility.kr",m.kr,"1/s"); row(a.id,"motility.kb",m.kb,"1/s"); row(a.id,"motility.m0",m.m0,"1");
            row(a.id,"motility.alpha",m.alpha,"1"); row(a.id,"motility.receptors",m.receptors,"1");
            row(a.id,"motility.active_binding",m.active_binding,"fmol/um^3"); row(a.id,"motility.inactive_binding",m.inactive_binding,"fmol/um^3");
            row(a.id,"motility.ya",m.ya,"fmol/um^3"); row(a.id,"motility.kd",m.kd,"fmol/um^3"); row(a.id,"motility.g",m.g,"1"); row(a.id,"motility.motor_frequency",m.motor_frequency,"1/s");
            row(a.id,"motility.gamma_shape",m.gamma_shape,"1","default from Sneddon 2012 Methods",m.gamma_shape==4 ? "literature default" : "user override; provisional");
            row(a.id,"motility.gamma_scale_degrees",m.gamma_scale_degrees,"degrees","default from Sneddon 2012 Methods",m.gamma_scale_degrees==18.32 ? "literature default" : "user override; provisional");
            row(a.id,"motility.max_step",m.max_step,"s","numerical resolution; docs/model.md","numerical choice");
        }
    }
    p.close();
}
void RunWriter::snapshot(const Universe& u) const {
    const auto key=suffix(u.step_index());
    auto agents=output(directory_/("agents_"+key+".csv"));
    agents << "step,time_s,agent_id,kind,alive,activated,x_um,y_um,z_um,biomass_fg,eps_fg,cell_volume_um3,total_radius_um,heading_x,heading_y,heading_z,methylation,running,run_seconds,tumble_seconds,switches,reorientations";
    for (auto name : chemical_names) agents << ',' << name << "_fmol," << name << "_fmol_per_um3";
    agents << '\n';
    for (const auto& a : u.agents()) {
        agents << u.step_index() << ',' << u.time() << ',' << a.id << ',' << static_cast<int>(a.kind) << ',' << a.alive << ',' << a.activated
            << ',' << a.position.x << ',' << a.position.y << ',' << a.position.z << ',' << a.biomass << ',' << a.eps << ',' << a.volume() << ',' << a.radius()
            << ',' << a.motion.heading.x << ',' << a.motion.heading.y << ',' << a.motion.heading.z << ',' << a.motion.methylation << ',' << a.motion.running
            << ',' << a.motion.run_seconds << ',' << a.motion.tumble_seconds << ',' << a.motion.switches << ',' << a.motion.reorientations;
        for (std::size_t k=0; k<chemical_count; ++k) agents << ',' << a.amounts[k] << ',' << a.concentration(static_cast<Chemical>(k));
        agents << '\n';
    }
    agents.close();
    auto fields=output(directory_/("field_"+key+".csv"));
    fields << "step,time_s,voxel_id,x_um,y_um,z_um,volume_um3,material,bulk,species_id,species,amount,amount_unit,concentration,concentration_unit,fixed_reservoir_concentration\n";
    for (std::size_t i=0; i<u.grids().size(); ++i) {
        const auto center=u.center(i); const auto& g=u.grids()[i];
        for (std::size_t s=0; s<u.species().size(); ++s) {
            const auto& species=u.species()[s]; const auto unit=amount_unit(species);
            fields << u.step_index() << ',' << u.time() << ',' << i << ',' << center.x << ',' << center.y << ',' << center.z << ',' << u.voxel_volume()
                << ',' << (g.material==Material::Biofilm ? "biofilm" : "boundary") << ',' << g.bulk << ',' << s << ',' << species.name
                << ',' << g.amounts[s] << ',' << unit << ',' << u.concentration(i,s) << ',' << unit << "/um^3,";
            if (g.fixed_concentrations[s]) fields << *g.fixed_concentrations[s]; else if (g.bulk) fields << u.config().bulk_concentrations[s];
            fields << '\n';
        }
    }
    fields.close();
    auto projection=output(directory_/("depth_average_"+key+".csv"));
    projection << "step,time_s,x_um,y_um,species,depth_average_concentration,concentration_unit\n";
    for (std::size_t x=0; x<u.config().shape[0]; ++x) for (std::size_t y=0; y<u.config().shape[1]; ++y)
        for (std::size_t s=0; s<u.species().size(); ++s) {
            double sum=0; for (std::size_t z=0; z<u.config().shape[2]; ++z) sum+=u.concentration(u.index(x,y,z),s);
            const auto center=u.center(u.index(x,y,0));
            projection << u.step_index() << ',' << u.time() << ',' << center.x << ',' << center.y << ',' << u.species()[s].name << ','
                << sum/static_cast<double>(u.config().shape[2]) << ',' << amount_unit(u.species()[s]) << "/um^3\n";
        }
    projection.close();
    auto d=output(directory_/("diagnostics_"+key+".json")); const auto& stats=u.diagnostics();
    d << "{\"step\":" << u.step_index() << ",\"time_s\":" << u.time() << ",\"particle_mass_fg\":" << u.particle_mass()
        << ",\"external_particle_mass_fg\":" << stats.external_particle_mass << ",\"maintenance_loss_fg\":" << stats.maintenance_loss
        << ",\"yield_loss_fg\":" << stats.yield_loss << ",\"biomass_height_um\":" << stats.biomass_height
        << ",\"max_overlap_um\":" << stats.max_overlap << ",\"mechanical_iterations\":" << stats.mechanical_iterations
        << ",\"transport_substeps\":" << stats.transport_substeps << ",\"ode_accepted\":" << stats.ode_accepted << ",\"ode_rejected\":" << stats.ode_rejected << ",\"ode_startup\":" << stats.ode_startup
        << ",\"births\":" << stats.births << ",\"deaths\":" << stats.deaths << ",\"shed_particles\":" << stats.shed_particles << ",\"species\":[";
    for (std::size_t s=0; s<u.species().size(); ++s)
        d << (s ? "," : "") << "{\"name\":" << quoted(u.species()[s].name) << ",\"unit\":" << quoted(amount_unit(u.species()[s]))
          << ",\"total\":" << u.chemical_total(s) << ",\"external\":" << stats.external[s] << ",\"reservoir\":" << stats.reservoir[s]
          << ",\"decay\":" << stats.decay[s] << ",\"reaction_source\":" << stats.reaction_source[s] << ",\"reaction_sink\":" << stats.reaction_sink[s]
          << ",\"budget_residual\":" << u.chemical_total(s)-u.expected_chemical_total(s) << '}';
    d << "]}\n"; d.close();
}
} // namespace BNSim
