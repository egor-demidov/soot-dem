#pragma once
// Two equal, non-overlapping DEM spheres with a cylindrical initial volume of
// free-surface liquid in the gap. Pure geometry and force decomposition: no
// SPHinXsys dependency, so the exact native setup can be unit-tested.
#include "pf_rigid.h"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace pf_bridge {
using namespace pf_liquid;
struct Geometry {
    double sphere_radius=.32;
    double gap=.24;                 // nearest surface-to-surface separation
    double liquid_radius=.4;       // transverse radius of initial liquid cylinder
    double penetration_factor=.65; // cylinder extent behind each sphere's tip / R
    void validate() const {
        if (!(sphere_radius>0 && gap>0 && liquid_radius>0 &&
              penetration_factor>0 && penetration_factor<=1 &&
              std::isfinite(sphere_radius) && std::isfinite(gap) &&
              std::isfinite(liquid_radius) && std::isfinite(penetration_factor)))
            throw std::invalid_argument("Invalid two-sphere bridge geometry");
    }
    double half_separation() const {return sphere_radius + .5*gap;}
    double half_cylinder_length() const {return .5*gap+penetration_factor*sphere_radius;}
    double solid_x_extent() const {return half_separation()+sphere_radius;}
    Vector center(size_t body) const {
        if(body>1) throw std::out_of_range("Two-sphere body index");
        return {body==0?-half_separation():half_separation(),0.,0.};
    }
    // Signed radial clearance outside sphere surface; negative is penetration.
    double clearance(Vector x,size_t body) const {
        return std::sqrt(dot(subtract(x,center(body)),subtract(x,center(body))))-sphere_radius;
    }
    bool in_solid(Vector x) const {return clearance(x,0)<0 || clearance(x,1)<0;}
    bool in_liquid(Vector x) const {
        return std::abs(x[0])<half_cylinder_length() &&
               x[1]*x[1]+x[2]*x[2]<liquid_radius*liquid_radius &&
               !in_solid(x);
    }
    size_t body_of(Vector x) const {
        if(!in_solid(x)) throw std::invalid_argument("Not a two-sphere solid site");
        return x[0]<0?0:1; // gap > 0 means surfaces cannot cross the midpoint
    }
};
struct BridgeForces {
    std::vector<pf_rigid::Load> pf_solid;
    double liquid_min_clearance=std::numeric_limits<double>::infinity();
    size_t neck_particles=0;
    size_t wetting_left=0,wetting_right=0;
    bool liquid_connected=false;
};
// Only the *fluid-solid* PF pair contribution, not total capillary bridge force.
// Total solid loads also include the WCSPH pressure and viscous forces.
inline BridgeForces diagnostics(const std::vector<Sample> &s,const Pairs &pairs,
                                const Parameters &p,const Geometry &g,
                                const std::vector<pf_rigid::Sphere> &b,
                                const std::vector<pf_rigid::Site> &sites) {
    if(b.size()!=2||b[0].radius!=g.sphere_radius||b[1].radius!=g.sphere_radius)
        throw std::runtime_error("Bridge diagnostic requires two equal spheres");
    std::vector<int> body_for_particle(s.size(),-1);
    for(const auto &site:sites) {
        if(site.particle>=s.size()||site.body>=2||!s[site.particle].solid||
           body_for_particle[site.particle]>=0)
            throw std::runtime_error("Invalid bridge site map");
        body_for_particle[site.particle]=static_cast<int>(site.body);
    }
    BridgeForces out;out.pf_solid.resize(2);
    // A short-range liquid graph avoids falsely treating gaps smaller than the
    // full four-dx force support as a continuous liquid bridge.
    std::vector<size_t> parent(s.size());
    for(size_t i=0;i<s.size();++i) parent[i]=i;
    auto root=[&](size_t i) {
        while(parent[i]!=i) {parent[i]=parent[parent[i]];i=parent[i];}
        return i;
    };
    for(const auto &[i,j]:pairs) {
        if(i>=s.size()||j>=s.size()||i==j) throw std::runtime_error("Invalid bridge pair");
        if(!s[i].solid&&!s[j].solid) {
            auto d=subtract(s[i].x,s[j].x);
            if(dot(d,d)<2.56*p.dx*p.dx) parent[root(i)]=root(j);
        }
    }
    std::vector<unsigned> wetting(s.size(),0);
    for(size_t i=0;i<s.size();++i) if(!s[i].solid) {
        const auto &a=s[i];
        const double d0=std::sqrt(dot(subtract(a.x,b[0].center),subtract(a.x,b[0].center)))-b[0].radius;
        const double d1=std::sqrt(dot(subtract(a.x,b[1].center),subtract(a.x,b[1].center)))-b[1].radius;
        out.liquid_min_clearance=std::min({out.liquid_min_clearance,d0,d1});
        if(std::abs(a.x[0])<p.dx) ++out.neck_particles;
        if(d0<=1.6*p.dx) {++out.wetting_left;wetting[root(i)]|=1;}
        if(d1<=1.6*p.dx) {++out.wetting_right;wetting[root(i)]|=2;}
    }
    for(size_t i=0;i<s.size();++i) if(parent[i]==i && wetting[i]==3)
        out.liquid_connected=true;
    for(const auto &[i,j]:pairs) {
        if(i>=s.size()||j>=s.size()||i==j) throw std::runtime_error("Invalid bridge pair");
        if(s[i].solid==s[j].solid) continue;
        const size_t solid=s[i].solid?i:j,fluid=s[i].solid?j:i;
        const int body=body_for_particle[solid];
        if(body<0) throw std::runtime_error("Bridge solid lacks body mapping");
        auto interaction=pf_capillary::evaluate(s[solid].x,s[fluid].x,p.support(),p.sfs,
                                                pf_capillary::PairKind::fluid_solid);
        auto &f=out.pf_solid[body];
        f.force=add(f.force,interaction.force_on_i);
        f.torque=add(f.torque,cross(subtract(s[solid].x,b[body].center),interaction.force_on_i));
    }
    return out;
}
} // namespace pf_bridge
