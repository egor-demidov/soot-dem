#pragma once
// Rigid spherical DEM boundary driven by reciprocal PF-WCSPH sample forces.
// No SPHinXsys dependency; solid sample mass is boundary quadrature weight,
// NOT added to the physical DEM sphere mass. Angular momentum is in lab frame.
#include "pf_liquid.h"
#include <cstddef>
#include <algorithm>
#include <array>
#include <stdexcept>

namespace pf_rigid {
using namespace pf_liquid;
using Quaternion=std::array<double,4>; // scalar-first, body -> world
inline Quaternion normalize(Quaternion q) {
    double n=0.; for(double a:q) n+=a*a;
    if(!(n>0&&std::isfinite(n))) throw std::runtime_error("Invalid DEM quaternion");
    for(double &a:q) a/=std::sqrt(n);
    return q;
}
inline Quaternion multiply(Quaternion a,Quaternion b) {
    return {a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0]};
}
inline Vector rotate(Quaternion q,Vector r) {
    q=normalize(q);
    Vector u{q[1],q[2],q[3]};
    const auto t=scale(cross(u,r),2.);
    return add(r,add(scale(t,q[0]),cross(u,t)));
}
inline Quaternion add_quaternion(Quaternion a,Quaternion b) {
    for(int j=0;j<4;++j) a[j]+=b[j];
    return a;
}
inline Quaternion scale_quaternion(Quaternion q,double factor) {
    for(double &x:q) x*=factor;
    return q;
}
inline Quaternion derivative(Quaternion q,Vector omega) {
    auto out=multiply(Quaternion{0.,omega[0],omega[1],omega[2]},q);
    for(double &x:out) x*=.5;
    return out;
}
struct Sphere {
    Vector center{},velocity{},angular_momentum{};
    Quaternion orientation{1.,0.,0.,0.};
    double radius=1.,mass=1.;
    bool fixed=false;
    double inertia() const {return .4*mass*radius*radius;}
    Vector omega() const {return scale(angular_momentum,1./inertia());}
    void validate() const {
        if(!finite(center)||!finite(velocity)||!finite(angular_momentum)||
           !(radius>0&&std::isfinite(radius))||!(mass>0&&std::isfinite(mass)))
            throw std::runtime_error("Invalid DEM sphere state");
        (void)normalize(orientation);
    }
};
// Each solid sample maps to exactly one sphere. Local position is invariant.
struct Site {size_t particle=0,body=0; Vector local{};};
struct Load {Vector force{},torque{};};
struct Contacts {
    std::vector<Load> loads;
    double energy=0.;
};
inline void validate_sites(const std::vector<Sample> &state,
                           const std::vector<Sphere> &bodies,const std::vector<Site> &sites) {
    for(const auto &body:bodies) body.validate();
    std::vector<int> visited(state.size(),0);
    for(auto &site:sites) {
        if(site.particle>=state.size()||!state[site.particle].solid||
           site.body>=bodies.size()||!finite(site.local)||visited[site.particle]++)
            throw std::runtime_error("Invalid/duplicate DEM boundary site");
    }
    for(size_t i=0;i<state.size();++i)
        if(state[i].solid && !visited[i]) throw std::runtime_error("Unassigned DEM solid sample");
}
inline void sync_sites(std::vector<Sample> &state,const std::vector<Sphere> &bodies,
                       const std::vector<Site> &sites) {
    for(const auto &site:sites) {
        auto &a=state[site.particle];
        const auto &b=bodies[site.body];
        const Vector arm=rotate(b.orientation,site.local);
        a.x=add(b.center,arm);
        a.v=b.fixed?Vector{}:add(b.velocity,cross(b.omega(),arm));
        // Boundary rho/mass are kept at their original quadrature values.
    }
}
inline std::vector<Load> loads_from_fluid(const std::vector<Sample> &state,
                              const Rates &rates,const std::vector<Sphere> &bodies,
                              const std::vector<Site> &sites) {
    if(rates.force.size()!=state.size()) throw std::runtime_error("DEM force array size mismatch");
    std::vector<Load> out(bodies.size());
    for(const auto &site:sites) {
        const auto &f=rates.force[site.particle];
        auto &l=out[site.body];
        l.force=pf_liquid::add(l.force,f);
        l.torque=pf_liquid::add(l.torque,cross(subtract(state[site.particle].x,bodies[site.body].center),f));
    }
    return out;
}
// Frictionless, reversible Hertz-type normal contact. k has force / length^1.5.
// This stage deliberately has no collision damping/tangential friction.
inline Contacts sphere_contacts(const std::vector<Sphere> &bodies,double k) {
    if(!(k>=0&&std::isfinite(k))) throw std::runtime_error("Invalid DEM contact stiffness");
    Contacts result;result.loads.resize(bodies.size());
    if(k==0) return result;
    for(size_t i=0;i<bodies.size();++i) for(size_t j=i+1;j<bodies.size();++j) {
        const Vector d=subtract(bodies[i].center,bodies[j].center);
        double r=std::sqrt(dot(d,d));
        double overlap=bodies[i].radius+bodies[j].radius-r;
        if(overlap<=0) continue;
        if(r==0) throw std::runtime_error("Coincident DEM sphere centers in contact");
        double amp=k*std::pow(overlap,1.5);
        const auto f=scale(d,amp/r);
        result.loads[i].force=add(result.loads[i].force,f);
        result.loads[j].force=subtract(result.loads[j].force,f);
        result.energy+=.4*k*std::pow(overlap,2.5);
    }
    return result;
}
inline std::vector<Load> total_loads(const std::vector<Sample> &s,const Rates &r,
                      const std::vector<Sphere> &b,const std::vector<Site> &sites,
                      double contact_stiffness=0.) {
    auto out=loads_from_fluid(s,r,b,sites);
    auto contacts=sphere_contacts(b,contact_stiffness);
    for(size_t i=0;i<out.size();++i) {
        out[i].force=add(out[i].force,contacts.loads[i].force);
        out[i].torque=add(out[i].torque,contacts.loads[i].torque);
    }
    return out;
}
inline Vector total_momentum(const std::vector<Sample> &s,const std::vector<Sphere> &b,
                             bool angular=false) {
    Vector out{};
    for(const auto &a:s) if(!a.solid) {
        auto p=scale(a.v,a.mass);
        out=add(out,angular?cross(a.x,p):p);
    }
    for(const auto &r:b) if(!r.fixed) {
        auto p=scale(r.velocity,r.mass);
        out=add(out,angular?add(cross(r.center,p),r.angular_momentum):p);
    }
    return out;
}
inline double timestep(const std::vector<Sample> &s,const std::vector<Sphere> &b,
                       const std::vector<Site> &sites,const Rates &rates,
                       const Parameters &p,double stiffness=0.) {
    double dt=pf_liquid::timestep(s,rates,p);
    const auto loads=total_loads(s,rates,b,sites,stiffness);
    // Heuristic explicit-contact stability bound: the Hertz tangent stiffness
    // at a characteristic overlap of dx is (3/2) k sqrt(dx). This is an
    // additional safeguard, not a replacement for timestep convergence tests.
    if(stiffness>0) for(size_t i=0;i<b.size();++i) if(!b[i].fixed) {
        const double k_tangent=1.5*stiffness*std::sqrt(p.dx);
        dt=std::min(dt,.08*std::sqrt(b[i].mass/k_tangent));
    }
    for(size_t i=0;i<b.size();++i) if(!b[i].fixed) {
        const auto &body=b[i];
        double vmax=std::sqrt(dot(body.velocity,body.velocity))+
                    std::sqrt(dot(body.omega(),body.omega()))*body.radius;
        double amag=std::sqrt(dot(loads[i].force,loads[i].force))/body.mass+
                    std::sqrt(dot(loads[i].torque,loads[i].torque))*body.radius/body.inertia();
        dt=std::min(dt,.15*p.dx/(p.sound+vmax));
        if(amag>0) dt=std::min(dt,.15*std::sqrt(p.dx/amag));
        auto w=body.omega();double wmag=std::sqrt(dot(w,w));
        if(wmag>0) dt=std::min(dt,.1/wmag);
    }
    return dt;
}
// Explicit midpoint forces; final translation uses trapezoidal velocity,
// consistent with PF liquid advance. Rigid boundary regenerated at each stage.
// The old and midpoint graph must be rebuilt for moving solid locations.
template<class Graph,class Evaluator=SerialEvaluator>
Rates advance(std::vector<Sample> &s,std::vector<Sphere> &bodies,
              const std::vector<Site> &sites,const Parameters &p,double dt,
              Graph graph,const Rates *initial=nullptr,Evaluator evaluator={},
              double contact_stiffness=0.) {
    if(!(dt>0&&std::isfinite(dt))) throw std::invalid_argument("Invalid coupled timestep");
    validate_sites(s,bodies,sites);
    const auto old=s;const auto old_b=bodies;
    Rates own_initial;
    if(!initial) {own_initial=evaluator(old,graph(old),p); initial=&own_initial;}
    if(initial->force.size()!=s.size()||initial->density.size()!=s.size())
        throw std::runtime_error("Invalid coupled initial rates");
    auto start_load=total_loads(old,*initial,old_b,sites,contact_stiffness);
    auto mid=old;auto mid_b=old_b;
    for(size_t i=0;i<s.size();++i) if(!old[i].solid) {
        mid[i].x=add(old[i].x,scale(old[i].v,.5*dt));
        mid[i].v=add(old[i].v,scale(initial->force[i],.5*dt/old[i].mass));
        mid[i].rho=old[i].rho+.5*dt*initial->density[i];
        validate(mid[i]);
    }
    for(size_t i=0;i<bodies.size();++i) if(!old_b[i].fixed) {
        auto &m=mid_b[i];const auto &o=old_b[i];
        m.center=add(o.center,scale(o.velocity,.5*dt));
        m.velocity=add(o.velocity,scale(start_load[i].force,.5*dt/o.mass));
        m.angular_momentum=add(o.angular_momentum,scale(start_load[i].torque,.5*dt));
        m.orientation=normalize(add_quaternion(o.orientation,scale_quaternion(derivative(o.orientation,o.omega()),.5*dt)));
    }
    sync_sites(mid,mid_b,sites);
    auto rates=evaluator(mid,graph(mid),p);
    auto mid_load=total_loads(mid,rates,mid_b,sites,contact_stiffness);
    auto next=old;auto next_b=old_b;
    for(size_t i=0;i<s.size();++i) if(!old[i].solid) {
        next[i].v=add(old[i].v,scale(rates.force[i],dt/old[i].mass));
        next[i].x=add(old[i].x,scale(add(old[i].v,next[i].v),.5*dt));
        next[i].rho=old[i].rho+dt*rates.density[i];
        validate(next[i]);
    }
    for(size_t i=0;i<bodies.size();++i) if(!old_b[i].fixed) {
        auto &n=next_b[i];const auto &o=old_b[i];
        n.velocity=add(o.velocity,scale(mid_load[i].force,dt/o.mass));
        n.angular_momentum=add(o.angular_momentum,scale(mid_load[i].torque,dt));
        n.center=add(o.center,scale(add(o.velocity,n.velocity),.5*dt));
        n.orientation=normalize(add_quaternion(o.orientation,scale_quaternion(derivative(mid_b[i].orientation,mid_b[i].omega()),dt)));
        n.validate();
    }
    sync_sites(next,next_b,sites);
    s.swap(next);bodies.swap(next_b);
    return rates;
}
} // namespace pf_rigid
