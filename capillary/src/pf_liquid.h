#pragma once
// Liquid-only WCSPH operators: Whebell et al. (2025), eqs. 2.4-2.17.
// SPHinXsys supplies neighbor lists to the driver; this core is independently testable.
#include "pf_capillary.h"
#include <algorithm>
#include <limits>
#include <vector>
#include <utility>
namespace pf_liquid {
using namespace pf_capillary;
struct Sample {
    Vector x{}, v{};
    double rho=1., mass=1.;
    bool solid=false;
    double quadrature_weight=1.; // transient union-volume weight, rigid samples only
};
struct Parameters {
    double dx=.12, rho0=1., sound=10., mu=.05, sff=.001, sfs=.001;
    double support() const { return 4*dx; }
    void validate() const {
        for(double v : {dx,rho0,sound})
            if(!(v>0 && std::isfinite(v))) throw std::invalid_argument("Invalid positive PF parameter");
        for(double v : {mu,sff,sfs})
            if(!(v>=0 && std::isfinite(v))) throw std::invalid_argument("Invalid PF strength/viscosity");
    }
};
using Pairs=std::vector<std::pair<size_t,size_t>>; // unique, current, undirected, no solid-solid
inline void validate(const Sample &s) {
    if(!finite(s.x)||!finite(s.v)||!(s.rho>0 && std::isfinite(s.rho))||
       !(s.mass>0 && std::isfinite(s.mass))||
       !(s.quadrature_weight>0&&s.quadrature_weight<=1&&std::isfinite(s.quadrature_weight))) throw std::runtime_error("Invalid PF liquid state");
}
inline double pressure(double rho, const Parameters &p) {
    if(!(rho>0 && std::isfinite(rho))) throw std::runtime_error("Invalid PF density");
    double out=p.sound*p.sound*p.rho0/7*(std::pow(rho/p.rho0,7)-1);
    if(!std::isfinite(out)) throw std::runtime_error("PF pressure overflow");
    return out; // No tensile-pressure clipping or density reset.
}
inline double kernel_derivative(double r,double H) {
    double q=r/H;
    if(q>=1) return 0.;
    double z=1-q;
    return (-210./(3.14159265358979323846*H*H*H*H))*q*z*z*z;
}
struct Rates {
    std::vector<Vector> force;
    std::vector<double> density, viscous_row, conservative_virial;
    double pair_energy=0., viscous_power=0.;
};
inline Rates empty_rates(size_t n) {
    return {std::vector<Vector>(n),std::vector<double>(n),
            std::vector<double>(n),std::vector<double>(n)};
}
inline std::vector<double> pressures(const std::vector<Sample> &s,const Parameters &p) {
    p.validate();
    std::vector<double> out; out.reserve(s.size());
    for(const auto &a:s) {validate(a); out.push_back(pressure(a.rho,p));}
    return out;
}
// Evaluate an undirected pair only once. The distance and kernel prefactors are
// shared between all three forces, avoiding a second sqrt/hypot in the PF law.
// The input graph must be unique; no locks or global projections are applied.
inline void accumulate_pair(const std::vector<Sample> &s,size_t i,size_t j,
                            const Parameters &p,const std::vector<double> &P,
                            double H,double dw_prefactor,Rates &out) {
    if(i>=s.size()||j>=s.size()||i==j) throw std::invalid_argument("Invalid PF pair index");
    const auto &a=s[i]; const auto &b=s[j];
    if(a.solid && b.solid) throw std::invalid_argument("Solid-solid pair in PF fluid graph");
    const auto rvec=subtract(a.x,b.x);
    const double r=std::sqrt(dot(rvec,rvec));
    if(r>=H) return;
    if(r==0) throw std::runtime_error("Coincident SPH samples");
    const auto e=scale(rvec,1/r), dv=subtract(a.v,b.v);
    const double q=r/H,z=1-q,dw=dw_prefactor*q*z*z*z;
    const double radial=dot(dv,e);
    const double weight=a.solid?a.quadrature_weight:(b.solid?b.quadrature_weight:1.);
    const double volume_product=a.mass*b.mass/(a.rho*b.rho);
    const double drag=-10*p.mu*volume_product*dw/r;
    const bool fs=a.solid||b.solid;
    const auto kind=fs?PairKind::fluid_solid:PairKind::fluid_fluid;
    const double strength=fs?p.sfs:p.sff;
    const double cap_magnitude=H*strength*profile(q,kind);
    const double cap_energy=strength*H*H*profile_integral_to_cutoff(q,kind);
    const auto conservative=scale(e,cap_magnitude-volume_product*(P[i]+P[j])*dw);
    const auto force=scale(subtract(conservative,scale(e,drag*radial)),weight);
    const double virial=weight*.5*dot(rvec,conservative);
    out.conservative_virial[i]+=virial; out.conservative_virial[j]+=virial;
    out.force[i]=add(out.force[i],force); out.force[j]=subtract(out.force[j],force);
    if(!a.solid) out.density[i]+=weight*a.rho*b.mass/b.rho*radial*dw;
    if(!b.solid) out.density[j]+=weight*b.rho*a.mass/a.rho*radial*dw;
    out.viscous_row[i]+=weight*drag/a.mass; out.viscous_row[j]+=weight*drag/b.mass;
    out.viscous_power-=weight*drag*radial*radial;
    out.pair_energy+=weight*cap_energy;
}
inline Rates evaluate_range(const std::vector<Sample> &s,const Pairs &pairs,
                            const Parameters &p,const std::vector<double> &P,
                            size_t begin,size_t end) {
    if(P.size()!=s.size()||begin>end||end>pairs.size())
        throw std::invalid_argument("Invalid PF pair range/pressure array");
    auto out=empty_rates(s.size());
    const double H=p.support();
    const double dw_prefactor=-210./(3.14159265358979323846*H*H*H*H);
    for(size_t k=begin;k<end;++k)
        accumulate_pair(s,pairs[k].first,pairs[k].second,p,P,H,dw_prefactor,out);
    return out;
}
inline void add_rates(Rates &dst,const Rates &src) {
    if(dst.force.size()!=src.force.size()) throw std::invalid_argument("PF rate sizes differ");
    for(size_t i=0;i<dst.force.size();++i) {
        dst.force[i]=add(dst.force[i],src.force[i]);
        dst.density[i]+=src.density[i];
        dst.viscous_row[i]+=src.viscous_row[i];
        dst.conservative_virial[i]+=src.conservative_virial[i];
    }
    dst.pair_energy+=src.pair_energy;
    dst.viscous_power+=src.viscous_power;
}
inline Rates evaluate(const std::vector<Sample> &s,const Pairs &pairs,const Parameters &p) {
    auto P=pressures(s,p);
    return evaluate_range(s,pairs,p,P,0,pairs.size());
}
struct SerialEvaluator {
    Rates operator()(const std::vector<Sample> &s,const Pairs &pairs,const Parameters &p) const {
        return evaluate(s,pairs,p);
    }
};
inline double timestep(const std::vector<Sample>&s,const Rates &r,const Parameters &p) {
    double amax=0,vmax=0,row=0;
    for(size_t i=0;i<s.size();++i) if(!s[i].solid) {
        amax=std::max(amax,std::sqrt(dot(r.force[i],r.force[i]))/s[i].mass);
        vmax=std::max(vmax,std::sqrt(dot(s[i].v,s[i].v)));
        row=std::max(row,r.viscous_row[i]);
    }
    double dt=.15*p.dx/(p.sound+vmax);
    if(p.mu>0) dt=std::min(dt,p.dx*p.dx*p.rho0/p.mu);
    if(amax>0) dt=std::min(dt,.15*std::sqrt(p.dx/amax));
    if(row>0) dt=std::min(dt,.25/row);
    return dt;
}
// graph(state) MUST refresh neighbors at both old and midpoint positions.
// Evaluator is injectable so the native driver can use a deterministic TBB
// shard reduction while the standalone core has no external dependency.
template<class Graph,class Evaluator=SerialEvaluator>
Rates advance(std::vector<Sample> &s,const Parameters &p,double dt,Graph graph,
              const Rates *initial=nullptr,Evaluator evaluator={}) {
    if(!(dt>0 && std::isfinite(dt))) throw std::invalid_argument("Invalid PF timestep");
    const auto old=s;
    Rates owned_initial;
    if(!initial) {owned_initial=evaluator(old,graph(old),p);initial=&owned_initial;}
    if(initial->force.size()!=s.size()||initial->density.size()!=s.size())
        throw std::invalid_argument("PF initial-rate size mismatch");
    auto mid=old;
    for(size_t i=0;i<s.size();++i) if(!s[i].solid) {
        mid[i].x=add(old[i].x,scale(old[i].v,.5*dt));
        mid[i].v=add(old[i].v,scale(initial->force[i],.5*dt/old[i].mass));
        mid[i].rho=old[i].rho+.5*dt*initial->density[i];
    }
    auto r=evaluator(mid,graph(mid),p);
    auto next=old;
    for(size_t i=0;i<s.size();++i) if(!s[i].solid) {
        next[i].v=add(old[i].v,scale(r.force[i],dt/old[i].mass));
        next[i].x=add(old[i].x,scale(add(old[i].v,next[i].v),.5*dt));
        next[i].rho=old[i].rho+dt*r.density[i];
        validate(next[i]);
    }
    s.swap(next);
    return r;
}
} // namespace pf_liquid
