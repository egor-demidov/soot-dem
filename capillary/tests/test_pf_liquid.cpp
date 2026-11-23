#include "../src/pf_liquid.h"
#include <iostream>
using namespace pf_liquid;
void require(bool b,const char *message) {if(!b) throw std::runtime_error(message);}
double norm(Vector x) {return std::sqrt(dot(x,x));}
Pairs graph(const std::vector<Sample> &s) {
    Pairs p;
    for(size_t i=0;i<s.size();++i) for(size_t j=i+1;j<s.size();++j)
        if(!(s[i].solid&&s[j].solid)) p.emplace_back(i,j);
    return p;
}
Vector momentum(const std::vector<Sample> &s,bool angular=false) {
    Vector p{};
    for(auto &a:s) if(!a.solid) p=add(p,angular?cross(a.x,scale(a.v,a.mass)):scale(a.v,a.mass));
    return p;
}
int main() {
    Parameters p; p.dx=.1;p.sff=.0001;p.sfs=.0002;
    std::vector<Sample> s={{{0,0,.1},{.1,.03,0},1.001,.001,false},
                           {{.11,.02,.1},{-.02,.01,.01},.999,.0012,false},
                           {{.03,.1,.15},{0,-.02,.03},1.002,.0009,false},
                           {{.01,.01,-.1},{},1,.001,true}};
    auto r=evaluate(s,graph(s),p);
    Vector F{},T{}; double work=0,internal=0;
    for(size_t i=0;i<s.size();++i) {
        F=add(F,r.force[i]); T=add(T,cross(s[i].x,r.force[i]));
        work+=dot(s[i].v,r.force[i]);
        if(!s[i].solid) internal+=s[i].mass*pressure(s[i].rho,p)/(s[i].rho*s[i].rho)*r.density[i];
    }
    require(norm(F)<1e-14&&norm(T)<1e-14,"Combined operator reciprocity");
    require(r.viscous_power<0,"Viscosity must dissipate");
    Parameters pressure_only=p;pressure_only.sff=pressure_only.sfs=pressure_only.mu=0;
    auto pr=evaluate(s,graph(s),pressure_only); double pw=0;
    for(size_t i=0;i<s.size();++i) pw+=dot(s[i].v,pr.force[i]);
    require(std::abs(pw+internal)<1e-14,"Pressure/continuity power identity");
    (void)work;
    auto P0=momentum(s),L0=momentum(s,true); Vector wallP{},wallL{};
    for(int n=0;n<4;++n) {
        double dt=std::min(1e-5,timestep(s,evaluate(s,graph(s),p),p));
        auto mid=advance(s,p,dt,graph);
        wallP=add(wallP,scale(mid.force[3],dt));
        wallL=add(wallL,scale(cross(s[3].x,mid.force[3]),dt));
    }
    require(norm(subtract(add(momentum(s),wallP),P0))<1e-17,"Integrated linear momentum");
    require(norm(subtract(add(momentum(s,true),wallL),L0))<1e-17,"Integrated angular momentum");
    require(s[3].rho==1&&norm(s[3].v)==0,"Fixed dummy wall state");
    require(s[0].mass==.001&&s[1].mass==.0012,"No mass restoration");
    // Parallel shard reduction must preserve pairwise momentum/torque and rates.
    // This tests the numerical reduction used by the TBB adapter without TBB.
    auto allpairs=graph(s);
    auto serial=evaluate(s,allpairs,p);
    auto pres=pressures(s,p);
    auto chunked=empty_rates(s.size());
    for(size_t shard=0;shard<3;++shard) {
        size_t first=allpairs.size()*shard/3;
        size_t last=allpairs.size()*(shard+1)/3;
        auto partial=evaluate_range(s,allpairs,p,pres,first,last);
        add_rates(chunked,partial);
    }
    for(size_t i=0;i<s.size();++i) {
        require(norm(subtract(serial.force[i],chunked.force[i]))<1e-13,"Shard force consistency");
        require(std::abs(serial.density[i]-chunked.density[i])<1e-13,"Shard continuity consistency");
        require(std::abs(serial.viscous_row[i]-chunked.viscous_row[i])<1e-13,"Shard viscosity consistency");
        require(std::abs(serial.conservative_virial[i]-chunked.conservative_virial[i])<1e-13,"Shard virial consistency");
    }
    require(std::abs(serial.pair_energy-chunked.pair_energy)<1e-13,"Shard energy consistency");
    require(std::abs(serial.viscous_power-chunked.viscous_power)<1e-13,"Shard dissipation consistency");
    // Small 3D liquid-only lattice; same midpoint core and refreshed graph as executable.
    std::vector<Sample> drop;
    for(int i=-2;i<=2;++i) for(int j=-2;j<=2;++j) for(int k=-2;k<=2;++k)
        if(i*i+j*j+k*k<=4) drop.push_back({{i*p.dx,j*p.dx,k*p.dx},{},1,.001,false});
    auto initP=momentum(drop),initL=momentum(drop,true);
    for(int n=0;n<2;++n) {
        auto rr=evaluate(drop,graph(drop),p); advance(drop,p,timestep(drop,rr,p),graph);
    }
    require(norm(subtract(momentum(drop),initP))<1e-17,"Liquid-only smoke momentum");
    require(norm(subtract(momentum(drop,true),initL))<1e-17,"Liquid-only smoke angular momentum");
    // Failure leaves caller's state intact.
    auto before=s; bool rejected=false;
    try {advance(s,p,-1.,graph);} catch(const std::invalid_argument&) {rejected=true;}
    require(rejected&&s[0].rho==before[0].rho,"Invalid step rejection");
    std::cout<<"PF liquid: pressure work, viscosity dissipation, full-step force/torque balance and 3D liquid-only smoke passed\n";
}
