#include "../src/pf_bridge.h"
#include "../src/pf_restart.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <stdexcept>
using namespace pf_bridge;
void require(bool v,const char *msg){if(!v) throw std::runtime_error(msg);}
double norm(Vector a){return std::sqrt(dot(a,a));}
Pairs graph(const std::vector<Sample> &s) {
    Pairs pairs;
    for(size_t i=0;i<s.size();++i) for(size_t j=i+1;j<s.size();++j)
        if(!(s[i].solid&&s[j].solid)) pairs.emplace_back(i,j);
    return pairs;
}
int main() {
    Geometry g; g.sphere_radius=.3;g.gap=.2;g.liquid_radius=.32;
    g.validate();
    require(std::abs(g.half_separation()-.4)<1e-14,"Sphere spacing");
    require(g.in_solid(g.center(0))&&g.in_solid(g.center(1)),"Both solid cores");
    require(g.in_liquid({0,.08,0})&&g.in_liquid({0,0,0}),"Connected initialized liquid neck");
    require(!g.in_liquid(g.center(0))&&!g.in_liquid({0,.5,0}),"Exclude solid and outside fluid");
    require(g.body_of({-.11,0,0})==0&&g.body_of({.11,0,0})==1,"Map solid sites to correct body");
    require(std::abs(g.clearance({0,0,0},0)-.1)<1e-14,"Geometric gap clearance");

    Parameters p;p.dx=.1;p.rho0=1;p.sound=10;p.mu=0;p.sff=.001;p.sfs=.003;
    pf_rigid::Sphere left,right;
    left.center=g.center(0);right.center=g.center(1);
    left.radius=right.radius=g.sphere_radius;
    left.mass=right.mass=2.0;
    std::vector<pf_rigid::Sphere> bodies{left,right};
    std::vector<Sample> state={
        {{0,.02,0},{},1.,.001,false},
        {{-.025,.115,0},{},1.,.001,false},
        {{.025,-.115,0},{},1.,.001,false},
        {{-.1,0,0},{},1.,.001,true},
        {{-.15,.16,0},{},1.,.001,true},
        {{.1,0,0},{},1.,.001,true},
        {{.15,-.16,0},{},1.,.001,true}
    };
    const std::vector<pf_rigid::Site> sites={
        {3,0,{.3,0,0}}, {4,0,{.25,.16,0}},
        {5,1,{-.3,0,0}}, {6,1,{-.25,-.16,0}}
    };
    pf_rigid::validate_sites(state,bodies,sites);
    pf_rigid::sync_sites(state,bodies,sites);
    auto pairs=graph(state);
    auto rates=evaluate(state,pairs,p);
    auto decomposition=diagnostics(state,pairs,p,g,bodies,sites);
    require(decomposition.neck_particles>0&&decomposition.liquid_min_clearance>0&&
            decomposition.wetting_left>0&&decomposition.wetting_right>0&&decomposition.liquid_connected,
            "No initially penetrated liquid and connected wetting neck");
    // Two separate fluid patches can both touch a sphere while not forming
    // a spanning liquid neck. The connectivity flag must distinguish them.
    auto separated=state;
    separated[0].x={-.18,.25,0};
    separated[1].x={.18,.25,0};
    separated[2].x={.18,.26,0};
    auto disconnected=diagnostics(separated,graph(separated),p,g,bodies,sites);
    require(!disconnected.liquid_connected,"Disconnected fluid patches are not a bridge");
    auto loads=pf_rigid::loads_from_fluid(state,rates,bodies,sites);
    // At rho=rho0 and zero viscosity, the total solid loads are exactly
    // the direct fluid-solid PF contribution. Liquid-fluid PF still acts internally.
    for(size_t j=0;j<2;++j) {
        require(norm(subtract(loads[j].force,decomposition.pf_solid[j].force))<1e-13,
                "Separately reported PF force matches direct pair sum");
        require(norm(subtract(loads[j].torque,decomposition.pf_solid[j].torque))<1e-13,
                "Separately reported PF torque matches direct pair sum");
    }
    // Relative rigid-body translations must generate the energy gradient.
    auto energy=[&](std::vector<pf_rigid::Sphere> b) {
        auto s=state;pf_rigid::sync_sites(s,b,sites);
        return evaluate(s,graph(s),p).pair_energy;
    };
    for(size_t j=0;j<2;++j) {
        auto plus=bodies,minus=bodies;
        constexpr double eps=1e-6;
        plus[j].center[0]+=eps;minus[j].center[0]-=eps;
        double gradient=-(energy(plus)-energy(minus))/(2*eps);
        require(std::abs(gradient-loads[j].force[0])<1e-9,
                "Two-body force agrees with capillary pair energy gradient");
    }
    // Two freely moving bodies + liquid exchange equal-and-opposite impulses.
    auto state0=state;
    bodies[0].velocity={-.01,0,0}; bodies[1].velocity={.01,0,0};
    pf_rigid::sync_sites(state,bodies,sites);
    const auto P0=pf_rigid::total_momentum(state,bodies);
    const auto L0=pf_rigid::total_momentum(state,bodies,true);
    for(int k=0;k<12;++k) (void)pf_rigid::advance(state,bodies,sites,p,1e-5,graph,nullptr,
                                                   SerialEvaluator{},10.);
    require(norm(subtract(pf_rigid::total_momentum(state,bodies),P0))<2e-13,
            "Two free spheres + liquid conserve linear momentum");
    require(norm(subtract(pf_rigid::total_momentum(state,bodies,true),L0))<1e-9,
            "Two free spheres + liquid conserve angular momentum");
    pf_restart::Checkpoint saved;
    saved.kind="bridge";saved.parameters=p;saved.radius=g.liquid_radius;
    saved.bridge_gap=g.gap;saved.contact_stiffness=10.;
    saved.nf=3;saved.ns=4;saved.spheres=bodies;saved.state=state;
    saved.referenceP=P0;saved.referenceL=L0;
    std::stringstream stream;pf_restart::write(stream,saved);
    auto restored=pf_restart::read(stream);
    require(restored.spheres.size()==2&&restored.bridge_gap==g.gap&&
            restored.contact_stiffness==10.,"Version 3 two-sphere checkpoint roundtrip");
    std::vector<pf_rigid::Sphere> expected{left,right};
    auto initial=state0;
    pf_restart::check_compatibility(restored,"bridge",p,g.liquid_radius,initial,3,
                                    &expected,&sites,g.gap,10.);
    bool rejected=false;
    try {pf_restart::check_compatibility(restored,"bridge",p,g.liquid_radius,initial,3,
                                         &expected,&sites,g.gap+.01,10.);}
    catch(const std::runtime_error&){rejected=true;}
    require(rejected,"Reject changed bridge gap on checkpoint resume");
    // Exact split/resume on the same deterministic graph.
    auto whole=state,split=state;
    auto whole_b=bodies,split_b=bodies;
    for(int k=0;k<8;++k) pf_rigid::advance(whole,whole_b,sites,p,1e-5,graph,nullptr,SerialEvaluator{},10.);
    for(int k=0;k<3;++k) pf_rigid::advance(split,split_b,sites,p,1e-5,graph,nullptr,SerialEvaluator{},10.);
    auto cp=saved;cp.state=split;cp.spheres=split_b;
    std::stringstream mid;pf_restart::write(mid,cp);
    auto cp_read=pf_restart::read(mid);
    split=cp_read.state;split_b=cp_read.spheres;
    for(int k=3;k<8;++k) pf_rigid::advance(split,split_b,sites,p,1e-5,graph,nullptr,SerialEvaluator{},10.);
    for(size_t i=0;i<whole.size();++i)
        require(whole[i].x==split[i].x&&whole[i].v==split[i].v&&whole[i].rho==split[i].rho,
                "Two sphere split/restart liquid state");
    for(int i=0;i<2;++i)
        require(whole_b[i].center==split_b[i].center&&whole_b[i].velocity==split_b[i].velocity&&
                whole_b[i].orientation==split_b[i].orientation,"Two sphere split/restart DEM state");

    // Fixed spheres' external supports must receive the negative impulse of
    // the fluid. Do not count the fixed bodies as freely moving momentum.
    auto fixed_state=state0;
    auto fixed_b=std::vector<pf_rigid::Sphere>{left,right};
    fixed_b[0].fixed=true;fixed_b[1].fixed=true;
    pf_rigid::sync_sites(fixed_state,fixed_b,sites);
    auto fluid_start=pf_rigid::total_momentum(fixed_state,fixed_b);
    auto angular_start=pf_rigid::total_momentum(fixed_state,fixed_b,true);
    Vector supportP{},supportL{};
    for(int i=0;i<8;++i) {
        auto r=pf_rigid::advance(fixed_state,fixed_b,sites,p,1e-5,graph,nullptr,SerialEvaluator{},10.);
        for(auto &site:sites) {
            auto f=scale(r.force[site.particle],1e-5);
            supportP=add(supportP,f);
            supportL=add(supportL,cross(fixed_state[site.particle].x,f));
        }
    }
    require(norm(subtract(add(pf_rigid::total_momentum(fixed_state,fixed_b),supportP),fluid_start))<1e-12,
            "Fixed supports close total linear momentum ledger");
    require(norm(subtract(add(pf_rigid::total_momentum(fixed_state,fixed_b,true),supportL),angular_start))<1e-10,
            "Fixed supports close total angular momentum ledger");
    std::cout<<"Two-sphere bridge geometry, PF loads, conservative free/fixed DEM, v3 checkpoint passed\n";
}
