#include "../src/pf_restart.h"
#include <iostream>
#include <sstream>
#include <cmath>
using namespace pf_rigid;
void require(bool condition,const char *what) {if(!condition) throw std::runtime_error(what);}
double norm(Vector a){return std::sqrt(dot(a,a));}
Pairs graph(const std::vector<Sample> &s){Pairs out;for(size_t i=0;i<s.size();++i)
  for(size_t j=i+1;j<s.size();++j) if(!(s[i].solid&&s[j].solid)) out.emplace_back(i,j);return out;}
int main() {
    Parameters p;p.dx=.1;p.sound=10;p.mu=.01;p.sff=.001;p.sfs=.003;
    Sphere b;b.center={0,0,0};b.radius=.2;b.mass=2.;b.velocity={0,.012,0};
    b.angular_momentum={0,0,.03};
    std::vector<Sphere> bodies={b};
    std::vector<Sample> state={
        {{.27,.03,.05},{-.01,0,.02},1.003,.001,false},
        {{.32,-.13,0},{.03,.001,0},.996,.001,false},
        {{.19,.0,.0},{},1,.001,true},
        {{0,.19,0},{},1,.001,true},
        {{-.19,0,0},{},1,.001,true}};
    std::vector<Site> sites={{2,0,{.19,0,0}},{3,0,{0,.19,0}},{4,0,{-.19,0,0}}};
    validate_sites(state,bodies,sites);
    sync_sites(state,bodies,sites);
    require(std::abs(state[2].v[1]-(b.velocity[1]+b.omega()[2]*.19))<1e-15,"Rigid rotation surface velocity");
    auto rates=evaluate(state,graph(state),p);
    auto loads=loads_from_fluid(state,rates,bodies,sites);
    Vector Fa{},Ta{};
    for(size_t i=2;i<state.size();++i) {
        Fa=add(Fa,rates.force[i]);
        Ta=add(Ta,cross(subtract(state[i].x,b.center),rates.force[i]));
    }
    require(norm(subtract(Fa,loads[0].force))<1e-14,"Force reduction");
    require(norm(subtract(Ta,loads[0].torque))<1e-14,"Torque reduction");
    const auto P0=total_momentum(state,bodies);
    const auto L0=total_momentum(state,bodies,true);
    for(int i=0;i<10;++i) (void)advance(state,bodies,sites,p,1e-5,graph);
    require(norm(subtract(total_momentum(state,bodies),P0))<5e-15,"Coupled linear momentum");
    require(norm(subtract(total_momentum(state,bodies,true),L0))<1e-9,"Coupled angular momentum");
    for(auto &site:sites) {
        const auto arm=rotate(bodies[0].orientation,site.local);
        require(norm(subtract(state[site.particle].x,add(bodies[0].center,arm)))<1e-13,"Rigid sample position");
        require(norm(subtract(state[site.particle].v,add(bodies[0].velocity,cross(bodies[0].omega(),arm))))<1e-13,"Rigid sample velocity");
        require(std::abs(dot(arm,arm)-dot(site.local,site.local))<1e-13,"Distance rigidity");
    }
    double qnorm=0;for(auto v:bodies[0].orientation) qnorm+=v*v;
    require(std::abs(qnorm-1)<1e-14,"Quaternion norm");
    // Frictionless Hertz contact is pairwise, central, and energy-consistent.
    Sphere a; a.center={-.95,0,0};a.radius=1.;a.mass=1.;
    Sphere c; c.center={.95,0,0};c.radius=1.;c.mass=1.;
    auto contact=sphere_contacts({a,c},25.);
    require(norm(add(contact.loads[0].force,contact.loads[1].force))<1e-15,"Hertz contact reciprocity");
    require(contact.loads[0].force[0]<0 && contact.energy>0,"Hertz repulsion and potential");
    double delta=.1;
    require(std::abs(contact.energy-.4*25*std::pow(delta,2.5))<1e-13,"Hertz potential");
    // Rigid motion without forces, including nonzero spin, must stay on a sphere.
    Parameters zero=p;zero.sfs=zero.sff=zero.mu=0;zero.sound=10.;
    std::vector<Sample> isolated={{{.19,0,0},{},1,.001,true}};
    std::vector<Sphere> spins={b};
    std::vector<Site> spin_sites={{0,0,{.19,0,0}}};
    sync_sites(isolated,spins,spin_sites);
    auto center0=spins[0].center,vel0=spins[0].velocity,Lspin0=spins[0].angular_momentum;
    (void)advance(isolated,spins,spin_sites,zero,.01,graph);
    require(norm(subtract(spins[0].center,add(center0,scale(vel0,.01))))<1e-15,"Force-free ballistic center");
    require(norm(subtract(spins[0].angular_momentum,Lspin0))<1e-15,"Torque-free angular momentum");
    // Capillary energy must generate the same force and torque on a rigid body.
    Parameters eparam=p;
    eparam.mu=0;eparam.sff=0;eparam.sfs=.01;
    std::vector<Sample> interaction={{{.29,.12,.07},{},1,.001,false},
                                     {{0,0,0},{},1,.001,true}};
    Sphere body;
    body.mass=3.;body.radius=.2;body.center={.015,-.01,.02};
    body.orientation=normalize({.96,.1,.2,.1});
    const std::vector<Site> one_site={{1,0,{.19,0,0}}};
    auto energy=[&](Sphere one) {
        auto ss=interaction;
        sync_sites(ss,{one},one_site);
        return evaluate(ss,graph(ss),eparam).pair_energy;
    };
    auto ss=interaction;
    sync_sites(ss,{body},one_site);
    auto force=evaluate(ss,graph(ss),eparam);
    const auto load=loads_from_fluid(ss,force,{body},one_site)[0];
    const double eps=1e-6;
    Sphere plus=body,minus=body;
    plus.center[0]+=eps;minus.center[0]-=eps;
    const double fx=-(energy(plus)-energy(minus))/(2*eps);
    require(std::abs(fx-load.force[0])<1e-9,"PF capillary body force is an energy gradient");
    plus=body;minus=body;
    plus.orientation=multiply({std::cos(eps/2),0,0,std::sin(eps/2)},body.orientation);
    minus.orientation=multiply({std::cos(eps/2),0,0,-std::sin(eps/2)},body.orientation);
    const double tz=-(energy(plus)-energy(minus))/(2*eps);
    require(std::abs(tz-load.torque[2])<1e-9,"PF capillary rigid torque is an energy gradient");
    // Restart an actual liquid/6DOF trajectory, not just individual fields.
    auto whole=state,part=state;
    auto wb=bodies,pb=bodies;
    for(int k=0;k<8;++k) advance(whole,wb,sites,p,1e-5,graph);
    for(int k=0;k<3;++k) advance(part,pb,sites,p,1e-5,graph);
    pf_restart::Checkpoint cp;cp.kind="immersed";cp.parameters=p;cp.radius=.7;
    cp.nf=2;cp.ns=3;cp.step=3;cp.time=3e-5;cp.spheres=pb;cp.state=part;
    std::stringstream saved;pf_restart::write(saved,cp);
    auto restored=pf_restart::read(saved);
    part=restored.state;pb=restored.spheres;
    for(int k=3;k<8;++k) advance(part,pb,sites,p,1e-5,graph);
    for(size_t i=0;i<whole.size();++i) {
        require(whole[i].x==part[i].x&&whole[i].v==part[i].v&&whole[i].rho==part[i].rho,
                "6DOF resumed fluid+site state exactly matches continuous run");
    }
    require(wb[0].center==pb[0].center && wb[0].velocity==pb[0].velocity &&
            wb[0].angular_momentum==pb[0].angular_momentum &&
            wb[0].orientation==pb[0].orientation,"6DOF trajectory split/restart");
    bool failed=false;
    try {validate_sites(state,bodies,{{2,0,{.19,0,0}},{2,0,{.19,0,0}}});}
    catch(const std::runtime_error&){failed=true;}
    require(failed,"Reject duplicate sites");
    std::cout<<"Rigid reciprocal PF coupling, 6DOF spin, midpoint conservation, and Hertz contact passed\n";
}
