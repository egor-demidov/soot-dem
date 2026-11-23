#pragma once
#include "pf_dem.h"
#include <limits>
namespace pf_dem {
struct DensityProbe {
    size_t sample=0;
    double rho_ratio=0,pressure=0,density_rate=0,continuity_divergence=0;
    Vector position{},velocity{};
    long long nearest_body=-1;
    double surface_gap=std::numeric_limits<double>::quiet_NaN();
    size_t liquid_neighbors=0,boundary_neighbors=0;
};
// Read-only probes at reporting times. Neighbor counts are geometric counts
// within kernel support; counts are unweighted; force rates include the union partition.
inline DensityProbe density_probe(const System &s,const Stage &rates,size_t i) {
    if(i>=s.nf)throw std::out_of_range("Density diagnostic sample index");
    const auto &a=s.state[i];DensityProbe p;p.sample=i;p.rho_ratio=a.rho/s.parameters.rho0;
    p.pressure=pf_liquid::pressure(a.rho,s.parameters);p.position=a.x;p.velocity=a.v;
    p.density_rate=rates.fluid.density.at(i);p.continuity_divergence=-p.density_rate/a.rho;
    double best=std::numeric_limits<double>::infinity();
    for(size_t j=0;j<s.bodies.size();++j){auto d=subtract(a.x,s.bodies[j].center);
        double gap=std::sqrt(dot(d,d))-s.bodies[j].radius;
        if(gap<best){best=gap;p.nearest_body=static_cast<long long>(j);p.surface_gap=gap;}}
    const double H=s.parameters.support(),H2=H*H;
    for(size_t j=0;j<s.nf;++j)if(j!=i){auto d=subtract(a.x,s.state[j].x);if(dot(d,d)<H2)++p.liquid_neighbors;}
    for(size_t j=s.nf;j<s.state.size();++j){auto d=subtract(a.x,s.state[j].x);if(dot(d,d)>=H2)continue;
        ++p.boundary_neighbors;
    }
    return p;
}
inline std::pair<size_t,size_t> density_extreme_indices(const System &s) {
    if(s.nf==0)throw std::runtime_error("No liquid density diagnostics");
    size_t lo=0,hi=0;for(size_t i=1;i<s.nf;++i){
        if(s.state[i].rho<s.state[lo].rho)lo=i;
        if(s.state[i].rho>s.state[hi].rho)hi=i;}
    return {lo,hi};
}
inline void write_density_header(std::ostream &f){
    f<<"step,time,last_dt,extreme,sample,rho_over_rho0,pressure_pa,drho_dt,continuity_divergence,x,y,z,vx,vy,vz,nearest_body,surface_gap,liquid_neighbors,boundary_neighbors\n";
}
inline void write_density_probe(std::ostream &f,const System &s,double dt,const char *label,const DensityProbe &p){
    f<<std::setprecision(17)<<s.step<<','<<s.time<<','<<dt<<','<<label<<','<<p.sample<<','<<p.rho_ratio<<','<<p.pressure<<','<<p.density_rate<<','<<p.continuity_divergence;
    for(auto v:{p.position,p.velocity})for(double x:v)f<<','<<x;
    f<<','<<p.nearest_body<<','<<p.surface_gap<<','<<p.liquid_neighbors<<','<<p.boundary_neighbors<<'\n';
}
}
