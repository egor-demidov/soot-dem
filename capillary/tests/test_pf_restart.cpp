#include "../src/pf_restart.h"
#include <sstream>
#include <cmath>
#include <iostream>
using namespace pf_restart;
void require(bool b,const char *msg) {if(!b) throw std::runtime_error(msg);}
int main() {
    Checkpoint a;
    a.kind="sessile"; a.parameters.dx=.04; a.radius=.6; a.time=7.0002326908245482;
    a.lastdt=.000299; a.step=23377; a.nf=2;a.ns=1;
    a.wallP={.002,1.34e-17,-.000019}; a.wallL={-1.e-15,2.e-16,3.e-14};
    a.state={{{.1,.2,.3},{.04,.01,0},1.025,.000064,false},
             {{.2,.4,.3},{-.04,0,0},.998,.000064,false},
             {{0,0,-.04},{0,0,0},1.,.000064,true}};
    std::stringstream f; write(f,a);
    auto b=read(f);
    require(b.time==a.time&&b.lastdt==a.lastdt&&b.step==a.step,"Clock roundtrip");
    require(b.wallP==a.wallP&&b.wallL==a.wallL,"Ledger roundtrip");
    require(b.state[0].x==a.state[0].x&&b.state[1].v==a.state[1].v,"State roundtrip");
    check_compatibility(b,a.kind,a.parameters,a.radius,a.state,a.nf);
    b.parameters.sff*=2.;
    bool rejected=false;
    try {check_compatibility(b,a.kind,a.parameters,a.radius,a.state,a.nf);}
    catch(const std::runtime_error&) {rejected=true;}
    require(rejected,"Parameter mismatch rejection");
    std::stringstream truncated("PFSPH_CHECKPOINT 1\nsessile\n");
    rejected=false;
    try {(void)read(truncated);} catch(const std::runtime_error&) {rejected=true;}
    require(rejected,"Truncated file rejection");
    // Splitting a time trajectory across an ASCII checkpoint reproduces the
    // exact integration state with a fixed graph and identical evaluator.
    auto graph=[](const std::vector<Sample>& s) {
        Pairs pairs;
        for(size_t i=0;i<s.size();++i) for(size_t j=i+1;j<s.size();++j)
            if(!(s[i].solid&&s[j].solid)) pairs.emplace_back(i,j);
        return pairs;
    };
    auto continuous=a.state;
    auto broken=a.state;
    const double dt=1e-5;
    for(int k=0;k<8;++k) (void)advance(continuous,a.parameters,dt,graph);
    for(int k=0;k<3;++k) (void)advance(broken,a.parameters,dt,graph);
    Checkpoint cp=a; cp.state=broken; cp.step=a.step+3; cp.time=a.time+3*dt;
    std::stringstream saved; write(saved,cp);
    auto restored=read(saved); broken=restored.state;
    for(int k=3;k<8;++k) (void)advance(broken,a.parameters,dt,graph);
    for(size_t i=0;i<broken.size();++i) {
        require(broken[i].x==continuous[i].x,"Restarted trajectory positions");
        require(broken[i].v==continuous[i].v,"Restarted trajectory velocities");
        require(broken[i].rho==continuous[i].rho,"Restarted trajectory density");
    }
    // Version 2 preserves a moving 6DOF body and the initial conserved baselines.
    Checkpoint moving=a;
    moving.kind="immersed";
    moving.referenceP={.023,.001,0};moving.referenceL={.0,-.001,.012};
    pf_rigid::Sphere sphere;
    sphere.center={.011,-.02,.03};sphere.radius=.18;sphere.mass=3.;
    sphere.velocity={.012,.003,-.004};sphere.angular_momentum={.1,.07,-.05};
    sphere.orientation=pf_rigid::normalize({.9,.2,.1,-.3});
    moving.spheres={sphere};
    // Initial lattice coordinates attached to the body, not its evolved position.
    auto base=a.state;
    moving.state[2].x=add(sphere.center,pf_rigid::rotate(sphere.orientation,base[2].x));
    moving.state[2].v=add(sphere.velocity,cross(sphere.omega(),
                                pf_rigid::rotate(sphere.orientation,base[2].x)));
    std::stringstream stream2;write(stream2,moving);
    auto restored2=read(stream2);
    require(restored2.referenceP==moving.referenceP && restored2.referenceL==moving.referenceL,
            "Version 2 momentum baseline");
    require(restored2.spheres[0].orientation==sphere.orientation,
            "Version 2 quaternion roundtrip");
    require(restored2.spheres[0].angular_momentum==sphere.angular_momentum,
            "Version 2 spin roundtrip");
    std::vector<pf_rigid::Sphere> expected={sphere};
    check_compatibility(restored2,moving.kind,moving.parameters,moving.radius,base,moving.nf,&expected);
    expected[0].mass*=2.;
    rejected=false;
    try {check_compatibility(restored2,moving.kind,moving.parameters,moving.radius,base,moving.nf,&expected);}
    catch(const std::runtime_error&){rejected=true;}
    require(rejected,"Reject changed DEM material mass");
    // Read existing v1 file from the previous stage (no v2 fields).
    std::stringstream legacy;
    legacy<<"PFSPH_CHECKPOINT 1\n"<<a.kind<<"\n";
    legacy<<a.parameters.dx<<" "<<a.parameters.rho0<<" "<<a.parameters.sound<<" "
          <<a.parameters.mu<<" "<<a.parameters.sff<<" "<<a.parameters.sfs<<" "<<a.radius<<"\n";
    legacy<<a.nf<<" "<<a.ns<<" "<<a.step<<" "<<a.time<<" "<<a.lastdt<<"\n";
    for(double v:a.wallP) legacy<<v<<" ";
    for(double v:a.wallL) legacy<<v<<" ";
    legacy<<"\n";
    for(auto &x:a.state) {
        for(double v:x.x) legacy<<v<<" ";
        for(double v:x.v) legacy<<v<<" ";
        legacy<<x.rho<<" "<<x.mass<<" "<<(x.solid?1:0)<<"\n";
    }
    auto old=read(legacy);
    require(old.spheres.empty() && old.step==a.step && old.state.size()==a.state.size(),
            "Backward-compatible v1 checkpoint recovery");
    std::cout<<"PF checkpoint: v1 recovery, v2 moving DEM state, conservation baseline and compatibility verified\n";
}
