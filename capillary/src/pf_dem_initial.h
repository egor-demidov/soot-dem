#pragma once
#include "pf_dem.h"
#include <optional>
namespace pf_dem {
struct Droplet {Vector center{};double radius=1.;};
struct FlatWall {
    double z=0.;
    double half_width=0.;
};
struct InitialGeometry {
    std::vector<Sphere> bodies;
    std::vector<Droplet> droplets;
    std::optional<FlatWall> wall;
    Vector lower{-1,-1,-1},upper{1,1,1};
    bool levelset=false;
    double alpha=.5,threshold=.5,rm=0.; // rm=0 uses each particle radius
    size_t max_samples=2000000;
    bool in_solid(Vector x) const {
        for(const auto &b:bodies){auto d=subtract(x,b.center);if(dot(d,d)<b.radius*b.radius)return true;}
        return false;
    }
    bool in_liquid(Vector x) const {
        if((wall && x[2]<=wall->z) || in_solid(x))return false;
        if(levelset){
            double psi=0.;
            for(const auto &b:bodies){auto d=subtract(x,b.center);double r=std::sqrt(dot(d,d));
                psi+=.5*std::erfc((r/(rm>0?rm:b.radius)-1)/alpha);}
            return psi>=threshold;
        }
        for(const auto &d:droplets){auto v=subtract(x,d.center);if(dot(v,v)<d.radius*d.radius)return true;}
        return false;
    }
    void validate() const {
        if(!finite(lower)||!finite(upper)||max_samples<1||
           !(alpha>0&&std::isfinite(alpha))||!(threshold>0&&std::isfinite(threshold))||
           !(rm>=0&&std::isfinite(rm)))throw std::invalid_argument("Invalid aggregate/liquid initialization");
        for(int k=0;k<3;++k)if(!(upper[k]>lower[k]))throw std::invalid_argument("Invalid sampling bounds");
        for(const auto &b:bodies)b.validate();
        if(!levelset&&droplets.empty())throw std::invalid_argument("No liquid droplets specified");
        for(const auto &d:droplets)if(!finite(d.center)||!(d.radius>0&&std::isfinite(d.radius)))throw std::invalid_argument("Invalid droplet");
        if(bodies.empty() && (levelset || droplets.size()!=1))
            throw std::invalid_argument("Calibration requires exactly one liquid droplet");
        if(wall) {
            if(!bodies.empty() || !std::isfinite(wall->z) ||
               !(wall->half_width>0 && std::isfinite(wall->half_width)))
                throw std::invalid_argument("Flat wall requires a sessile calibration without DEM particles");
            const auto &d=droplets.front();
            if(!(std::abs(d.center[2]-wall->z)<d.radius))
                throw std::invalid_argument("Initial droplet must intersect the flat wall plane");
        }
        // Explicit bounds must contain the complete solid geometry.
        for(const auto &b:bodies)for(int k=0;k<3;++k)
            if(b.center[k]-b.radius<lower[k]||b.center[k]+b.radius>upper[k])throw std::invalid_argument("Bounds truncate DEM sphere");
    }
};
inline System initialize(const Parameters &p,const InitialGeometry &g,ContactModel cm={}) {
    p.validate();g.validate();std::vector<Sample>s;std::vector<Site>sites;
    const double mass=p.rho0*p.dx*p.dx*p.dx;
    auto insert=[&](Vector x,bool solid){if(s.size()>=g.max_samples)throw std::runtime_error("Initial sample limit exceeded; increase max_samples deliberately");
        s.push_back({x,{},p.rho0,mass,solid});};
    // Guard iteration count as well as stored samples (a large empty box can be expensive).
    long long count[3];double first[3];long double work=1;
    const Vector origin=g.bodies.empty()?g.droplets.front().center:Vector{};
    for(int k=0;k<3;++k) {
        first[k]=(std::floor((g.lower[k]-origin[k])/p.dx)+.5)*p.dx+origin[k];
        if(g.wall && k==2)first[k]=g.wall->z+.5*p.dx;
        const double cells=std::ceil((g.upper[k]-first[k])/p.dx);
        if(!(cells>0&&cells<10000000))throw std::runtime_error("Invalid/oversized fluid sampling dimension");
        count[k]=static_cast<long long>(cells);
        work*=count[k];
    }
    if(work>50.L*g.max_samples)throw std::runtime_error("Sampling box too large for max_samples budget");
    for(long long i=0;i<count[0];++i)for(long long j=0;j<count[1];++j)for(long long k=0;k<count[2];++k){
        Vector x{first[0]+i*p.dx,first[1]+j*p.dx,first[2]+k*p.dx};if(g.in_liquid(x))insert(x,false);}
    size_t nf=s.size();if(!nf)throw std::runtime_error("Initial liquid is unresolved/empty; reduce dx or change geometry");
    for(size_t j=0;j<g.bodies.size();++j){const auto &b=g.bodies[j];double cells=std::ceil(b.radius/p.dx);
        if(cells>100000)throw std::runtime_error("Oversized DEM sample extent");
        auto n=static_cast<long long>(cells);size_t begin=s.size();
        if(std::pow(2.L*n+1,3)>50.L*g.max_samples)throw std::runtime_error("DEM sampling exceeds work budget");
        for(long long ix=-n;ix<=n;++ix)for(long long iy=-n;iy<=n;++iy)for(long long iz=-n;iz<=n;++iz){
            Vector local{ix*p.dx,iy*p.dx,iz*p.dx};if(dot(local,local)>=b.radius*b.radius)continue;
            sites.push_back({s.size(),j,local});insert(add(b.center,rotate(b.orientation,local)),true);}
        if(begin==s.size())throw std::runtime_error("Unresolved DEM boundary");
    }
    if(g.wall) {
        const auto &wall=*g.wall;
        const auto &center=g.droplets.front().center;
        const double half_cells=std::ceil(wall.half_width/p.dx);
        if(wall.half_width<g.droplets.front().radius+p.support() ||
           !(half_cells>0 && 4.L*std::ceil(p.support()/p.dx)*half_cells*half_cells<=g.max_samples-nf))
            throw std::runtime_error("Flat wall extent is too small for the droplet or exceeds max_samples");
        const auto n=static_cast<long long>(half_cells);
        const auto layers=static_cast<int>(std::ceil(p.support()/p.dx));
        for(long long i=-n;i<n;++i)
            for(long long j=-n;j<n;++j)
                for(int k=0;k<layers;++k)
                    insert({center[0]+(i+.5)*p.dx,center[1]+(j+.5)*p.dx,
                            wall.z-(k+.5)*p.dx},true);
    }
    return System(p,std::move(s),nf,g.bodies,std::move(sites),cm);
}
} // namespace pf_dem
