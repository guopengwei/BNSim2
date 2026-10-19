// Scientific regression checks. Assertions stay enabled in release builds.
#include "bnsim/io.hpp"
#include "numerics.hpp"
#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

using namespace BNSim;
namespace {
struct Test { std::string group, name; std::function<void()> run; };
std::vector<Test>& tests() { static std::vector<Test> all; return all; }
struct Register { Register(const char* group,const char* name,std::function<void()> f) { tests().push_back({group,name,std::move(f)}); } };
#define TEST(group,name) void name(); Register register_##name(#group,#name,name); void name()
#define CHECK(expression) do { if (!(expression)) throw std::runtime_error(std::string("check failed: ")+ #expression + " at line " + std::to_string(__LINE__)); } while (false)
void near(double a,double b,double absolute=1e-12,double relative=1e-8) {
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b)>absolute+relative*std::max(std::abs(a),std::abs(b))) {
        std::ostringstream out; out << std::setprecision(17) << a << " != " << b << " (atol=" << absolute << ", rtol=" << relative << ')';
        throw std::runtime_error(out.str());
    }
}
template<class F> void throws(F f) { bool caught=false; try { f(); } catch (const std::exception&) { caught=true; } CHECK(caught); }
SimulationConfig box() { SimulationConfig c; c.extent={10,10,10}; c.shape={1,1,1}; c.mechanics=false; return c; }
std::vector<Species> chemicals() { return {{"substrate",Quantity::Mass},{"A1"},{"A2"},{"QSI1"},{"QSI2"},{"QSI3"},{"attractant"}}; }
QSParameters inert_qs() {
    QSParameters q; q.signal1=1; q.signal2=2; q.inhibitors={3,4,5}; q.basal.fill(0); q.induced.fill(0); q.rate.fill(0); q.exchange.fill(0);
    q.inhibitor_binding=q.inhibitor_unbinding=q.signal_degradation=0; return q;
}
AgentSpec cell(Vec3 position={5,5,5}) {
    AgentSpec a; a.position=position; a.division=false; a.shedding=false; return a;
}
void budgets(const Universe& u, double relative=1e-8) {
    const auto& d=u.diagnostics();
    for (std::size_t s=0; s<u.species().size(); ++s) {
        const double throughput=std::abs(d.external[s])+std::abs(d.reservoir[s])+d.reaction_source[s]+d.reaction_sink[s]+d.decay[s];
        near(u.chemical_total(s),u.expected_chemical_total(s),1e-17+32*std::numeric_limits<double>::epsilon()*throughput,relative);
    }
}
std::string state(const Universe& u) {
    std::ostringstream out; out << std::hexfloat << u.time() << ' ' << u.step_index();
    for (const auto& g : u.grids()) {
        out << ' ' << static_cast<int>(g.material) << ' ' << g.bulk;
        for (auto n : g.amounts) out << ' ' << n;
        for (const auto& c : g.fixed_concentrations) out << ' ' << c.has_value() << ' ' << c.value_or(0);
        for (auto id : g.members) out << ' ' << id;
    }
    for (const auto& a : u.agents()) {
        out << ' ' << a.id << ' ' << static_cast<int>(a.kind) << ' ' << a.position.x << ' ' << a.position.y << ' ' << a.position.z
            << ' ' << a.biomass << ' ' << a.eps << ' ' << a.biomass_density << ' ' << a.eps_density << ' ' << a.alive << ' ' << a.activated
            << ' ' << a.division << ' ' << a.shedding << ' ' << a.division_mass << ' ' << a.death_mass << ' ' << a.shedding_mass;
        for (auto n : a.amounts) out << ' ' << n;
        out << ' ' << a.motion.heading.x << ' ' << a.motion.heading.y << ' ' << a.motion.heading.z << ' ' << a.motion.methylation
            << ' ' << a.motion.running << ' ' << a.motion.hazard_remaining << ' ' << a.motion.run_seconds << ' ' << a.motion.tumble_seconds
            << ' ' << a.motion.switches << ' ' << a.motion.reorientations << ' ' << a.random;
    }
    const auto& d=u.diagnostics();
    for (const auto* v : {&d.external,&d.reservoir,&d.decay,&d.reaction_source,&d.reaction_sink}) for (auto n : *v) out << ' ' << n;
    out << ' ' << d.external_particle_mass << ' ' << d.maintenance_loss << ' ' << d.yield_loss << ' ' << d.biomass_height
        << ' ' << d.max_overlap << ' ' << d.transport_substeps << ' ' << d.ode_accepted << ' ' << d.ode_rejected << ' ' << d.ode_startup
        << ' ' << d.mechanical_iterations << ' ' << d.births << ' ' << d.deaths << ' ' << d.shed_particles;
    return out.str();
}

TEST(validation, geometry_and_units) {
    auto c=box(); c.extent={10,11,13}; c.shape={3,5,7}; Universe u(c,{{"signal"}});
    near(u.voxel_volume(),10.*11*13/105);
    CHECK(u.index(Vec3{u.spacing().x*.9,0,0})==0);
    CHECK(u.index(Vec3{u.spacing().x*1.01,0,0})==u.index(1,0,0));
    CHECK(u.index(c.extent)==u.grids().size()-1);
    near(u.center(0).x,5./3); near(from_micromolar(1),1e-6); near(from_millimolar(1),1e-3);
    near(from_grams_per_liter(150),150); near(bimolecular_from_per_micromolar_second(.1),1e5);
    throws([&]{u.index(Vec3{-1,0,0});}); throws([&]{u.index(3,0,0);});
    throws([&]{u.index(Vec3{std::numeric_limits<double>::quiet_NaN(),0,0});});
}
TEST(validation, reject_invalid_inputs) {
    auto c=box(); c.shape[0]=0; throws([&]{Universe u(c,{});});
    c=box(); c.extent.z=std::numeric_limits<double>::infinity(); throws([&]{Universe u(c,{});});
    c=box(); c.shape.fill(std::numeric_limits<std::size_t>::max()); throws([&]{Universe u(c,{});});
    c=box(); c.workers=0; throws([&]{Universe u(c,{});});
    c=box(); c.layered=true; throws([&]{Universe u(c,{{"s"}});});
    throws([&]{Universe u(box(),{{"s"},{"s"}});});
    throws([&]{Universe u(box(),{{"s",Quantity::Molar,-1}});});
    Universe u(box(),chemicals());
    throws([&]{u.set_concentration(0,0,-1);}); throws([&]{u.set_concentration(0,99,1);});
    throws([&]{u.set_concentration(0,0,std::numeric_limits<double>::infinity());});
    throws([&]{u.species_id("missing");}); throws([&]{u.step(0);}); throws([&]{u.step(-1);});
    throws([&]{u.step(std::numeric_limits<double>::quiet_NaN());});
    auto a=cell(); a.qs=QSParameters{}; throws([&]{u.add_agent(a);}); // mass species as signal
    a=cell(); a.position={0,0,0}; throws([&]{u.add_agent(a);});
    a=cell(); a.growth=GrowthParameters{}; a.growth->pmax=1; throws([&]{u.add_agent(a);});
    a=cell(); a.growth=GrowthParameters{}; a.growth->yield=1.01; throws([&]{u.add_agent(a);});
    CHECK(u.agents().empty());
}
TEST(validation, empty_world_and_final_time) {
    for (auto workers : {std::size_t{1},std::size_t{8}}) {
        auto c=box(); c.workers=workers; Universe u(c,{}); u.run_until(.25,.1);
        CHECK(u.time()==.25); CHECK(u.step_index()==3); u.validate_state();
        throws([&]{u.run_until(.1,.1);});
    }
}
TEST(validation, independent_universes_and_ids) {
    Universe a(box(),{}),b(box(),{}); const auto id=a.add_agent(cell());
    CHECK(id==0); CHECK(b.add_agent(cell())==0); a.remove_agent(id); CHECK(a.add_agent(cell())==1);
    CHECK(b.agents().size()==1); CHECK(a.diagnostics().external_particle_mass==150);
    throws([&]{a.agent(id);}); a.validate_state();
}
TEST(validation, worker_exceptions_join) {
    std::array<int,17> done{};
    throws([&]{detail::parallel_for(done.size(),8,[&](std::size_t i){ done[i]=1; if (i==3) throw std::runtime_error("worker failure"); });});
    CHECK(std::all_of(done.begin(),done.end(),[](int x){return x==1;}));
}

TEST(transport, zero_decay_invariance) {
    auto c=box(); c.shape={3,2,4}; c.workers=8;
    Universe u(c,{{"s",Quantity::Molar,100,2,0,0}});
    for (std::size_t i=0;i<u.grids().size();++i) {u.set_concentration(i,0,1); if(i%2)u.set_material(i,Material::Biofilm);}
    u.step(.1); for (std::size_t i=0;i<u.grids().size();++i) CHECK(u.concentration(i,0)==1);
    budgets(u);
}
TEST(transport, analytic_decay_without_diffusion) {
    auto c=box(); c.diffusion=false; Universe u(c,{{"s",Quantity::Molar,0,0,.3,.3}});
    u.set_concentration(0,0,2); u.run_until(2.3,.17);
    near(u.concentration(0,0),2*std::exp(-.3*2.3)); budgets(u);
}
TEST(transport, closed_conservation_and_positivity) {
    auto c=box(); c.extent={7,5,3}; c.shape={5,3,2};
    Universe u(c,{{"s",Quantity::Mass,50,2,0,0}}); u.set_concentration(0,0,1);
    for(std::size_t i=0;i<u.grids().size();i+=3)u.set_material(i,Material::Biofilm);
    const double before=u.chemical_total(0); u.step(.7);
    near(u.chemical_total(0),before,1e-12,1e-12); budgets(u);
    CHECK(u.diagnostics().transport_substeps>10);
    for (const auto& g:u.grids())CHECK(g.amounts[0]>=0);
}
TEST(transport, harmonic_interfaces_and_impermeable_material) {
    auto c=box(); c.extent={2,1,1}; c.shape={2,1,1};
    Universe u(c,{{"s",Quantity::Molar,1,9,0,0}});u.set_material(1,Material::Biofilm);u.set_concentration(0,0,1);
    u.step(.01); const double interface=1.8;
    near(u.concentration(0,0),.5+.5*std::pow(1-2*interface*.005,2),1e-15,1e-14);
    Universe wall(c,{{"s",Quantity::Molar,1,0,0,0}});wall.set_material(1,Material::Biofilm);wall.set_concentration(0,0,1);wall.step(1);
    CHECK(wall.concentration(0,0)==1); CHECK(wall.concentration(1,0)==0);
}
TEST(transport, reservoir_and_decay_accounting) {
    auto c=box(); c.shape={2,1,1};
    Universe u(c,{{"s",Quantity::Mass,2,2,.2,.2}}); u.set_reservoir(0,0,1);u.run_until(2,.1);
    CHECK(u.concentration(0,0)==1);CHECK(u.concentration(1,0)>0);CHECK(u.diagnostics().reservoir[0]>u.voxel_volume());budgets(u);
    throws([&]{u.set_concentration(0,0,2);});u.set_reservoir(0,0,std::nullopt);u.step(.1);CHECK(u.concentration(0,0)<1);budgets(u);
}
TEST(transport, explicit_withdrawals_after_reservoir_release) {
    Universe u(box(),{{"s",Quantity::Mass}});u.set_reservoir(0,0,1);u.set_reservoir(0,0,std::nullopt);
    u.set_concentration(0,0,0);u.set_concentration(0,0,.25);budgets(u);
    CHECK(u.diagnostics().external[0]<0);
}
double cosine_error(std::size_t n,double dt,bool discrete) {
    auto c=box(); c.extent={10,3,7}; c.shape={n,1,1}; Universe u(c,{{"s",Quantity::Molar,1,1,0,0}});
    for(std::size_t i=0;i<n;++i)u.set_concentration(i,0,1+.5*std::cos(pi*u.center(i).x/10));
    u.run_until(.4,dt);
    const double rate=discrete ? 4*std::pow(std::sin(pi/(2*static_cast<double>(n)))/u.spacing().x,2) : pi*pi/100;
    double error=0;
    for(std::size_t i=0;i<n;++i)error=std::max(error,std::abs(u.concentration(i,0)-(1+.5*std::cos(pi*u.center(i).x/10)*std::exp(-rate*.4))));
    return error;
}
TEST(transport, temporal_and_spatial_convergence) {
    const double t1=cosine_error(12,.1,true),t2=cosine_error(12,.05,true),t3=cosine_error(12,.025,true);
    CHECK(t1/t2>1.9);CHECK(t2/t3>1.9);
    const double x1=cosine_error(10,.001,false),x2=cosine_error(20,.00025,false),x3=cosine_error(40,.0000625,false);
    CHECK(x1/x2>3.8);CHECK(x2/x3>3.8);
}
TEST(transport, transport_failure_rolls_back) {
    auto c=box();c.shape={2,2,2};c.max_transport_substeps=1;Universe u(c,{{"s",Quantity::Molar,1000,1000,1,1}});
    u.set_concentration(0,0,1);const auto before=state(u);throws([&]{u.step(10);});CHECK(state(u)==before);
}

TEST(biology, reversible_binding_invariants) {
    Universe u(box(),chemicals());auto a=cell();a.qs=inert_qs();auto& q=*a.qs;
    q.rate[1]=q.rate[4]=q.rate[8]=q.rate[11]=1e5;q.rate[2]=q.rate[5]=q.rate[9]=q.rate[12]=.1;
    q.inhibitor_binding=1e5;q.inhibitor_unbinding=.1;
    a.initial_concentrations={2e-6,1e-6,.2e-6,.1e-6,0,1.3e-6,.9e-6,.3e-6,.1e-6,.5e-6,.2e-6,0,0};
    const auto id=u.add_agent(a);const auto n=u.agent(id).amounts;u.run_until(12,1);const auto z=u.agent(id).amounts;
    near(z[A1]+z[RA1]+2*z[C1],n[A1]+n[RA1]+2*n[C1],1e-19);
    near(z[R1]+z[RA1]+2*z[C1]+z[RQ1],n[R1]+n[RA1]+2*n[C1]+n[RQ1],1e-19);
    near(z[Q1]+z[RQ1],n[Q1]+n[RQ1],1e-19);
    near(z[A2]+z[RA2]+2*z[C2],n[A2]+n[RA2]+2*n[C2],1e-19);
    near(z[R2]+z[RA2]+2*z[C2],n[R2]+n[RA2]+2*n[C2],1e-19);budgets(u);
}
TEST(biology, membrane_exchange_analytical) {
    Universe u(box(),chemicals());auto a=cell();a.qs=inert_qs();a.qs->exchange.fill(.5);
    for(auto k:{A1,A2,Q1,Q2,Q3})a.initial_concentrations[k]=1e-6;
    const auto id=u.add_agent(a);u.run_until(2,.3);
    const double equilibrium=1e-6/(1001),expected=equilibrium+(1e-6-equilibrium)*std::exp(-.5*1.001*2);
    for(auto k:{A1,A2,Q1,Q2,Q3})near(u.agent(id).amounts[k],expected,1e-17,1e-8);
    budgets(u);
}
TEST(biology, shared_voxel_limited_exchange) {
    Universe u(box(),chemicals());auto a=cell({3,5,5});a.qs=inert_qs();a.qs->exchange.fill(10);u.add_agent(a);
    a.position={7,5,5};u.add_agent(a);u.set_concentration(0,1,1e-9);u.step(3);
    near(u.agents()[0].amounts[A1],u.agents()[1].amounts[A1],1e-20);budgets(u);
    CHECK(u.chemical_total(1)<=1.000001e-6);CHECK(u.diagnostics().ode_rejected>0);
}
TEST(biology, starvation_and_maintenance_budget) {
    Universe u(box(),chemicals());auto a=cell();a.growth=GrowthParameters{};a.growth->maintenance=.1;u.add_agent(a);u.step(2);
    near(u.agents()[0].biomass,150*std::exp(-.2),1e-7);CHECK(u.agents()[0].eps==0);CHECK(u.chemical_total(0)==0);
    near(u.particle_mass()+u.diagnostics().maintenance_loss,150,1e-7);
}
TEST(biology, substrate_limited_growth_and_losses) {
    Universe u(box(),chemicals());u.set_concentration(0,0,.001);auto a=cell({3,5,5});a.growth=GrowthParameters{};
    a.growth->umax=1;a.growth->pmax=0;a.growth->maintenance=.01;a.growth->substrate_half=.1;
    u.add_agent(a);a.position={7,5,5};u.add_agent(a);u.step(10);budgets(u);
    CHECK(u.chemical_total(0)>=0);CHECK(u.diagnostics().reaction_sink[0]<=1+1e-12);
    near(u.chemical_total(0)+u.particle_mass()+u.diagnostics().maintenance_loss+u.diagnostics().yield_loss,301,1e-6);
}
TEST(biology, growth_dilution_and_radius) {
    Universe u(box(),chemicals());u.set_reservoir(0,0,1);auto a=cell();a.qs=inert_qs();a.initial_concentrations[A1]=1e-6;
    a.growth=GrowthParameters{};a.growth->umax=.2;a.growth->pmax=0;a.growth->maintenance=0;const auto id=u.add_agent(a);
    const double radius=u.agent(id).radius();u.step(1);
    near(u.agent(id).biomass,150*std::exp(.1),1e-7);near(u.agent(id).amounts[A1],1e-6,1e-20);
    near(u.agent(id).concentration(A1),1e-6*std::exp(-.1),1e-17);
    near(u.agent(id).radius()/radius,std::exp(.1/3),1e-9);budgets(u);
}
TEST(biology, fixed_reservoir_during_zero_state_startup) {
    Universe u(box(),chemicals());const double target=.123456789123456;
    u.set_reservoir(0,0,target);auto a=cell();a.qs=QSParameters{};a.qs->signal1=1;a.qs->signal2=2;a.growth=GrowthParameters{};
    u.add_agent(a);const double amount=u.grids()[0].amounts[0];u.step(.7);
    CHECK(u.grids()[0].amounts[0]==amount);CHECK(u.diagnostics().ode_startup>0);budgets(u);
}
TEST(biology, circuits_feedback_activation_and_inhibitors) {
    auto make=[](int inhibitor,bool feedback=true){
        Universe u(box(),chemicals());auto a=cell();a.qs=QSParameters{};auto& q=*a.qs;q.signal1=1;q.signal2=2;q.inhibitors={3,4,5};q.exchange.fill(0);
        q.basal[1]=1e-7;q.induced[1]=2e-7;q.activation_threshold=1e-12;q.rsal_feedback=feedback;
        a.initial_concentrations[A1]=a.initial_concentrations[R1]=a.initial_concentrations[A2]=a.initial_concentrations[R2]=1e-6;
        if(inhibitor>=0)a.initial_concentrations[std::array<Chemical,3>{Q1,Q2,Q3}[static_cast<std::size_t>(inhibitor)]]=1e-5;
        u.add_agent(a);u.run_until(20,1);budgets(u);return u;
    };
    const auto normal=make(-1),competition=make(0),degradation=make(1),suppression=make(2),no_feedback=make(-1,false);
    const auto& a=normal.agents()[0];CHECK(a.amounts[C1]>0);CHECK(a.amounts[C2]>0);CHECK(a.amounts[RsaL]>0);CHECK(a.activated);
    CHECK(competition.agents()[0].amounts[C1]<a.amounts[C1]);
    CHECK(degradation.agents()[0].amounts[A1]<a.amounts[A1]);
    CHECK(suppression.agents()[0].amounts[R1]<a.amounts[R1]);
    CHECK(no_feedback.diagnostics().reaction_source[1]>normal.diagnostics().reaction_source[1]);
}
TEST(biology, eps_production_and_allocation_budget) {
    Universe u(box(),chemicals());u.set_concentration(0,0,1);auto a=cell();a.qs=inert_qs();a.initial_concentrations[C1]=1e-6;
    a.growth=GrowthParameters{};a.growth->umax=.2;a.growth->pmax=.1;u.add_agent(a);u.step(2);
    CHECK(u.agents()[0].eps>0);
    near(u.chemical_total(0)+u.particle_mass()+u.diagnostics().maintenance_loss+u.diagnostics().yield_loss,1150,1e-6);budgets(u);
}
TEST(biology, conservative_division_and_inheritance) {
    auto c=box();c.extent={30,30,30};Universe u(c,chemicals());auto a=cell({15,15,15});a.biomass=600;a.eps=40;a.division=true;
    a.qs=inert_qs();a.qs->activation_threshold=3e-7;a.initial_concentrations[A1]=1e-6;a.initial_concentrations[C1]=.5e-6;
    u.add_agent(a);const double total=u.chemical_total(1);u.step(.01);CHECK(u.agents().size()==4);CHECK(u.diagnostics().births==3);
    near(u.particle_mass(),640);near(u.chemical_total(1),total,1e-20);budgets(u);
    for(const auto& daughter:u.agents()) {
        near(daughter.biomass,150);near(daughter.eps,10);near(daughter.concentration(A1),1e-6,1e-20);
        CHECK(daughter.qs->activation_threshold==3e-7);CHECK(daughter.qs->rate==a.qs->rate);
    }
}
TEST(biology, shedding_uses_mass_and_density) {
    Universe u(box(),chemicals());auto a=cell();a.eps=200;a.shedding=true;a.shedding_mass=100;u.add_agent(a);u.step(.1);
    CHECK(u.agents().size()==2);CHECK(u.agents()[0].eps==0);CHECK(u.agents()[1].kind==AgentKind::EPS);
    near(u.agents()[1].radius(),std::cbrt(3*(200./75)/(4*pi)));near(u.particle_mass(),350);CHECK(u.diagnostics().shed_particles==1);
}
TEST(biology, dead_particles_inactive_and_retained) {
    Universe u(box(),chemicals());auto a=cell();a.biomass=.5;a.qs=inert_qs();a.qs->basal[0]=1e-9;u.add_agent(a);u.step(.1);
    CHECK(!u.agents()[0].alive);const auto mass=u.agents()[0].biomass;const auto amounts=u.agents()[0].amounts;u.step(2);
    CHECK(u.agents().size()==1);CHECK(u.agents()[0].biomass==mass);CHECK(u.agents()[0].amounts==amounts);
}
TEST(biology, integrator_failure_is_transactional) {
    auto c=box();c.shape={2,2,2};c.workers=8;c.integration.max_attempts=1;Universe u(c,chemicals());
    auto a=cell({2,2,2});a.qs=inert_qs();a.qs->rate[1]=1e14;a.initial_concentrations[A1]=a.initial_concentrations[R1]=1e-6;u.add_agent(a);
    const auto before=state(u);throws([&]{u.step(1);});CHECK(state(u)==before);
}
TEST(biology, zero_state_reaction_cascade_startup) {
    std::vector<double> y(3),atol(3,1e-14);IntegrationOptions options;options.relative_tolerance=1e-9;
    const auto result=detail::integrate(y,1.,atol,options,[](const auto& state,auto& d){d[0]=1;d[1]=state[0]*state[0];d[2]=state[1]*state[1];});
    near(y[0],1,1e-12);near(y[1],1./3,1e-10);near(y[2],1./63,1e-10);
    CHECK(result.startup>0);
}

QSParameters reference_qs() {
    auto q=inert_qs();q.basal.fill(2e-9);q.induced.fill(2e-8);q.threshold.fill(.2e-6);
    q.rate={.03,.7e6,.2,.04,.3e6,.1,.02,.02,.6e6,.25,.03,.2e6,.08};q.exchange={.2,.25,.3,.35,.4};
    q.inhibitor_binding=.5e6;q.inhibitor_unbinding=.12;q.signal_degradation=.2e6;q.inhibitor_saturation=1e6;q.receptor_suppression=.3e-6;return q;
}
TEST(biology, independent_reference_trajectories) {
    Universe u(box(),chemicals());const std::array<double,6> ext{.4,.5e-6,.2e-6,.8e-6,.4e-6,.6e-6};
    for(std::size_t i=0;i<ext.size();++i)u.set_concentration(0,i,ext[i]);
    auto a=cell();a.eps=5;a.qs=reference_qs();a.growth=GrowthParameters{};
    a.growth->umax=.03;a.growth->pmax=.01;a.growth->maintenance=.002;a.growth->substrate_half=.2;a.growth->qs_half=.1e-6;a.growth->yield=.8;
    a.initial_concentrations={1e-6,.7e-6,.2e-6,.1e-6,.05e-6,.8e-6,.6e-6,.1e-6,.03e-6,.4e-6,.05e-6,.2e-6,.3e-6};u.add_agent(a);
    std::ifstream in(std::string(BNSIM_TEST_DATA)+"/coupled.csv");CHECK(in.good());std::string line;std::getline(in,line);
    std::size_t rows=0;
    while(std::getline(in,line)) {
        std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double time;row>>time;u.run_until(time,.37);double expected;
        for(std::size_t s=0;s<6;++s){CHECK(bool(row>>expected));near(u.grids()[0].amounts[s],expected,1e-15,3e-8);}
        CHECK(bool(row>>expected));near(u.agents()[0].biomass,expected,1e-6,3e-8);
        CHECK(bool(row>>expected));near(u.agents()[0].eps,expected,1e-6,3e-8);
        for(auto n:u.agents()[0].amounts){CHECK(bool(row>>expected));near(n,expected,1e-15,3e-8);}
        budgets(u,1e-7);++rows;
    }
    CHECK(rows==6);
}

void check_contacts(const Universe& u) {
    u.validate_state();
    for(std::size_t i=0;i<u.agents().size();++i)for(std::size_t j=i+1;j<u.agents().size();++j)
        CHECK(norm(u.agents()[i].position-u.agents()[j].position)+u.config().mechanical_tolerance*1.01>=u.agents()[i].radius()+u.agents()[j].radius());
}
TEST(mechanics, coincident_diagonal_and_large_contacts) {
    auto c=box();c.extent={30,30,30};c.shape={15,15,15};c.mechanics=true;Universe u(c,{});
    auto a=cell({10,10,10});u.add_agent(a);u.add_agent(a);a.position={11,11,11};u.add_agent(a);
    a.position={15,15,15};a.biomass=150*(4*pi/3)*27;u.add_agent(a);a.position={19,19,15};u.add_agent(a);u.step(.1);check_contacts(u);
}
TEST(mechanics, wall_contacts_and_reflection) {
    Vec3 p{35,-17,21},heading{1,0,0};reflect_position(p,heading,{10,10,10},1);
    CHECK(p.x>=1&&p.x<=9);CHECK(p.y>=1&&p.y<=9);CHECK(p.z>=1&&p.z<=9);CHECK(heading.x==1);
    auto c=box();c.mechanics=true;Universe u(c,{});auto a=cell({1,1,1});u.add_agent(a);u.add_agent(a);u.step(.1);check_contacts(u);
}
TEST(mechanics, movement_membership_uses_committed_position_once) {
    auto c=box();c.extent={20,10,10};c.shape={10,1,1};Universe u(c,{{"attractant"}});
    auto a=cell({2,5,5});a.motility=MotilityParameters{};a.motility->motor_frequency=0;a.motility->rotational_diffusion=0;a.motility->speed=5;a.initial_heading=Vec3{1,0,0};
    const auto id=u.add_agent(a);u.step(1);near(u.agent(id).position.x,7);CHECK(u.grids()[u.index(Vec3{7,5,5})].members.front()==id);u.validate_state();
}
TEST(mechanics, confinement_failure_rolls_back) {
    auto c=box();c.extent={2,2,2};c.mechanics=true;c.mechanical_iterations=30;Universe u(c,{});auto a=cell({1,1,1});
    a.biomass=150*4*pi*.9*.9*.9/3;u.add_agent(a);u.add_agent(a);const auto before=state(u);throws([&]{u.step(.1);});CHECK(before==state(u));
}
TEST(mechanics, vertical_layers_and_reservoir_reclassification) {
    auto c=box();c.extent={20,20,20};c.shape={1,1,10};c.layered=true;c.bulk_concentrations={2};c.boundary_layer=2;Universe u(c,{{"s",Quantity::Mass}});
    auto a=cell({10,10,3});const auto id=u.add_agent(a);CHECK(u.grids()[0].material==Material::Biofilm);CHECK(!u.grids()[2].bulk);CHECK(u.grids()[3].bulk);
    const double initial=u.diagnostics().reservoir[0];u.set_concentration(2,0,.5);u.remove_agent(id);
    CHECK(u.grids()[2].bulk);CHECK(u.concentration(2,0)==2);CHECK(u.diagnostics().reservoir[0]>initial);budgets(u);
}

Agent motor_agent(std::uint64_t seed=0) {
    Agent a;a.biomass=150;a.position={50,50,50};a.motility=MotilityParameters{};a.motility->g=0;a.motility->kr=a.motility->kb=0;
    a.motility->speed=0;a.motility->rotational_diffusion=0;a.motility->max_step=1;
    a.random.seed(seed);std::exponential_distribution<double> e(1);a.motion.hazard_remaining=e(a.random);return a;
}
TEST(motility, isotropic_initial_headings) {
    std::mt19937_64 random(0);Vec3 mean{},square{};const int n=30000;
    for(int i=0;i<n;++i){const auto v=isotropic_direction(random);mean=mean+v;for(std::size_t k=0;k<3;++k)square[k]+=v[k]*v[k];near(norm(v),1);}
    for(std::size_t k=0;k<3;++k){near(mean[k]/n,0,.015,0);near(square[k]/n,1./3,.012,0);}
}
TEST(motility, gaussian_rotation_correlation) {
    std::mt19937_64 random(4);double mean=0;const int n=10000;
    for(int i=0;i<n;++i){Vec3 h{0,0,1};for(int j=0;j<50;++j)rotational_increment(h,.2,.01,random);mean+=h.z;}
    near(mean/n,std::exp(-2*.2*.5),.012,0);
}
TEST(motility, bounded_gamma_statistics) {
    MotilityParameters p;std::mt19937_64 random(8);double sum=0,squares=0;const int n=50000;
    for(int i=0;i<n;++i){double angle=tumble_angle(p,random)*180/pi;CHECK(angle>=0&&angle<=180);sum+=angle;squares+=angle*angle;}
    // Independent numerical quadrature of theta^3 exp(-theta/scale) on [0,180].
    double z=0,first=0,second=0;const double dx=.01;
    for(int i=0;i<18000;++i){const double x=(i+.5)*dx,w=x*x*x*std::exp(-x/p.gamma_scale_degrees)*dx;z+=w;first+=w*x;second+=w*x*x;}
    near(sum/n,first/z,.7,0);near(squares/n-std::pow(sum/n,2),second/z-std::pow(first/z,2),25,0);
}
TEST(motility, continuous_time_motor_statistics) {
    auto a=motor_agent(5);advance_motility(a,0,20000,{100,100,100});
    const double runs=static_cast<double>((a.motion.switches+1)/2),tumbles=static_cast<double>(a.motion.switches/2);
    near(a.motion.run_seconds/runs,1/1.3,.025,0);near(a.motion.tumble_seconds/tumbles,1/1.3,.025,0);
    CHECK(a.motion.reorientations==a.motion.switches/2);near(a.motion.run_seconds+a.motion.tumble_seconds,20000,1e-7);
    int unchanged=0,running=0;const int n=10000;
    for(int i=0;i<n;++i){auto b=motor_agent(static_cast<std::uint64_t>(i));advance_motility(b,0,.7,{100,100,100});unchanged+=b.motion.switches==0;running+=b.motion.running;}
    near(static_cast<double>(unchanged)/n,std::exp(-1.3*.7),.02,0);
    near(static_cast<double>(running)/n,.5+.5*std::exp(-2*1.3*.7),.02,0);
}
TEST(motility, stationary_tumbles_and_single_reorientation) {
    auto a=motor_agent();a.motion.running=false;a.motility->speed=20;a.motility->motor_frequency=0;advance_motility(a,0,1,{100,100,100});
    CHECK(a.position.x==50&&a.position.y==50&&a.position.z==50);CHECK(a.motion.reorientations==0);CHECK(a.motion.run_seconds==0);
    a.motility->motor_frequency=1;a.motion.hazard_remaining=1;advance_motility(a,0,1,{100,100,100});
    CHECK(a.position.x==50&&a.position.y==50&&a.position.z==50);CHECK(a.motion.reorientations==1);CHECK(a.motion.running);
}
TEST(motility, attractant_response_and_adaptation) {
    MotilityParameters p;const double low=receptor_activity(p,2,0),high=receptor_activity(p,2,1e-3);
    CHECK(high<low);CHECK(motor_rates(p,high)[0]<motor_rates(p,low)[0]);
    auto a=motor_agent();a.motility=p;a.motility->speed=0;a.motility->motor_frequency=0;advance_motility(a,1e-3,500,{100,100,100});
    near(receptor_activity(p,a.motion.methylation,1e-3),.5,1e-8);
    const auto tiny=receptor_activity(p,1e300,1e300);CHECK(std::isfinite(tiny));
}
TEST(motility, ensemble_drift_up_attractant_gradient) {
    const int count=256;
    std::array<double,2> sum{},squares{};
    for (int setting=0;setting<2;++setting) for (int i=0;i<count;++i) {
        auto a=motor_agent(stream_seed(42,static_cast<AgentId>(i)));
        a.position={1000,1000,1000};a.motility=MotilityParameters{};
        a.motion.heading=isotropic_direction(a.random);
        const auto& p=*a.motility;
        const double initial_ligand=5e-5;
        a.motion.methylation=p.m0-(std::log1p(initial_ligand/p.active_binding)-std::log1p(initial_ligand/p.inactive_binding))/p.alpha;
        for(int step=0;step<400;++step) {
            const double ligand=initial_ligand*std::exp(setting ? (a.position.x-1000)/100 : 0);
            advance_motility(a,ligand,.1,{2000,2000,2000});
        }
        const double dx=a.position.x-1000;sum[setting]+=dx;squares[setting]+=dx*dx;
    }
    std::array<double,2> mean{},error{};
    for(std::size_t k=0;k<2;++k){mean[k]=sum[k]/count;error[k]=std::sqrt((squares[k]/count-mean[k]*mean[k])/count);}
    std::cout << "  drift: uniform=" << mean[0] << " +/- " << error[0] << "; gradient=" << mean[1] << " +/- " << error[1] << " um (standard error)\n";
    CHECK(std::abs(mean[0])<4*error[0]);CHECK(mean[1]>3*error[1]);
}

Universe mixed_world(std::size_t workers) {
    auto c=box();c.extent={40,20,30};c.shape={4,2,3};c.workers=workers;c.mechanics=true;c.seed=123;
    auto species=chemicals();for(auto& s:species){s.diffusion=2;s.biofilm_diffusion=.8;}
    Universe u(c,species);for(std::size_t i=0;i<u.grids().size();++i){u.set_concentration(i,0,.5);u.set_concentration(i,6,.0001+u.center(i).x*1e-5);}
    for(std::size_t i=0;i<8;++i){auto a=cell({5+8.*static_cast<double>(i%4),5+8.*static_cast<double>(i/4),10});
        a.qs=reference_qs();a.growth=GrowthParameters{};a.growth->umax=.6;a.growth->pmax=.1;a.division=true;a.division_mass=170;a.shedding=true;a.shedding_mass=1;
        a.initial_concentrations[C1]=1e-6;a.eps=.5;a.motility=MotilityParameters{};a.motility->attractant=6;a.motility->speed=1;u.add_agent(a);}
    auto dead=cell({20,15,20});dead.alive=false;u.add_agent(dead);
    AgentSpec eps;eps.kind=AgentKind::EPS;eps.biomass=0;eps.eps=150;eps.position={30,15,20};u.add_agent(eps);return u;
}
TEST(reproducibility, identical_across_worker_counts) {
    auto a=mixed_world(1),b=mixed_world(2),c=mixed_world(8);
    for(int i=0;i<5;++i){a.step(.2);b.step(.2);c.step(.2);CHECK(state(a)==state(b));CHECK(state(a)==state(c));}
    CHECK(a.diagnostics().births>0);CHECK(a.diagnostics().shed_particles>0);check_contacts(a);budgets(a,1e-7);
}
TEST(reproducibility, seed_controls_motility_only) {
    auto c=box();c.seed=1;Universe a(c,{{"s"}});c.seed=2;Universe b(c,{{"s"}});auto p=cell();p.motility=MotilityParameters{};
    a.add_agent(p);b.add_agent(p);CHECK(dot(a.agents()[0].motion.heading,b.agents()[0].motion.heading)<.999999);
    auto x=mixed_world(1),y=mixed_world(1);x.step(.1);y.step(.1);CHECK(state(x)==state(y));
}
} // namespace
int main(int argc,char** argv) {
    std::size_t count=0,failed=0;
    for(const auto& test:tests())if(argc==1||test.group==argv[1]){
        ++count;try{test.run();std::cout<<"PASS "<<test.name<<'\n';}catch(const std::exception& e){++failed;std::cerr<<"FAIL "<<test.name<<": "<<e.what()<<'\n';}
    }
    if(!count){std::cerr<<"No matching test group\n";return 2;}
    std::cout<<count-failed<<'/'<<count<<" passed\n";return failed ? 1 : 0;
}
