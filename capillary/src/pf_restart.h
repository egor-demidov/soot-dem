#pragma once
// Versioned restart including exact position/velocity/density, clock and wall
// impulse ledger. No SPHinXsys dependency; does not serialize neighbor lists.
#include "pf_rigid.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <cstdio>
#include <string>
namespace pf_restart {
using namespace pf_liquid;
struct Checkpoint {
    std::string kind;
    Parameters parameters;
    double radius=0, time=0, lastdt=0;
    double bridge_gap=0, contact_stiffness=0; // v3; zeros for v1/v2
    size_t nf=0, ns=0;
    long long step=0;
    Vector wallP{},wallL{};
    Vector referenceP{},referenceL{}; // Conserved initial total momentum about origin.
    std::vector<pf_rigid::Sphere> spheres;
    std::vector<Sample> state;
};
inline void validate_checkpoint(const Checkpoint &c) {
    c.parameters.validate();
    if((c.kind!="sessile" && c.kind!="droplet" && c.kind!="immersed" && c.kind!="bridge") ||
       !(c.radius>c.parameters.dx && std::isfinite(c.radius)) ||
       !(c.time>=0&&std::isfinite(c.time)) ||
       !(c.lastdt>=0&&std::isfinite(c.lastdt)) ||
       !(c.bridge_gap>=0&&std::isfinite(c.bridge_gap)) ||
       !(c.contact_stiffness>=0&&std::isfinite(c.contact_stiffness)) ||
       c.step<0 || c.nf==0 || c.nf+c.ns!=c.state.size() ||
       !finite(c.wallP) || !finite(c.wallL) || !finite(c.referenceP) || !finite(c.referenceL) ||
       (c.kind=="immersed" && c.spheres.size()!=1) ||
       (c.kind=="bridge" && (c.spheres.size()!=2 || c.bridge_gap<=0)) ||
       (c.kind!="immersed" && c.kind!="bridge" && !c.spheres.empty()))
        throw std::runtime_error("Invalid PF checkpoint metadata");
    for(const auto &b:c.spheres) b.validate();
    for(size_t i=0;i<c.state.size();++i) {
        validate(c.state[i]);
        if(c.state[i].solid!=(i>=c.nf))
            throw std::runtime_error("PF checkpoint particle ordering mismatch");
    }
}
inline void write(std::ostream &f,const Checkpoint &c) {
    validate_checkpoint(c);
    f<<std::setprecision(std::numeric_limits<double>::max_digits10);
    f<<"PFSPH_CHECKPOINT 3\n"<<c.kind<<'\n';
    const auto &p=c.parameters;
    f<<p.dx<<' '<<p.rho0<<' '<<p.sound<<' '<<p.mu<<' '<<p.sff<<' '<<p.sfs<<' '<<c.radius<<' '<<c.bridge_gap<<' '<<c.contact_stiffness<<'\n';
    f<<c.nf<<' '<<c.ns<<' '<<c.step<<' '<<c.time<<' '<<c.lastdt<<'\n';
    for(double v:c.wallP) f<<v<<' ';
    for(double v:c.wallL) f<<v<<' ';
    f<<'\n';
    for(double v:c.referenceP) f<<v<<' ';
    for(double v:c.referenceL) f<<v<<' ';
    f<<'\n'<<c.spheres.size()<<'\n';
    for(const auto &b:c.spheres) {
        for(double v:b.center) f<<v<<' ';
        for(double v:b.velocity) f<<v<<' ';
        for(double v:b.angular_momentum) f<<v<<' ';
        for(double v:b.orientation) f<<v<<' ';
        f<<b.radius<<' '<<b.mass<<' '<<(b.fixed?1:0)<<'\n';
    }
    for(const auto &a:c.state) {
        for(double x:a.x) f<<x<<' ';
        for(double v:a.v) f<<v<<' ';
        f<<a.rho<<' '<<a.mass<<' '<<(a.solid?1:0)<<'\n';
    }
    if(!f) throw std::runtime_error("PF checkpoint write failed");
}
inline Checkpoint read(std::istream &f) {
    std::string magic; int version=0;
    if(!(f>>magic>>version) || magic!="PFSPH_CHECKPOINT" || (version!=1&&version!=2&&version!=3))
        throw std::runtime_error("Invalid/unsupported PF checkpoint format");
    Checkpoint c;
    auto &p=c.parameters;
    if(!(f>>c.kind>>p.dx>>p.rho0>>p.sound>>p.mu>>p.sff>>p.sfs>>c.radius))
        throw std::runtime_error("Truncated PF checkpoint parameters");
    if(version>=3 && !(f>>c.bridge_gap>>c.contact_stiffness))
        throw std::runtime_error("Truncated PF bridge options");
    if(!(f>>c.nf>>c.ns>>c.step>>c.time>>c.lastdt))
        throw std::runtime_error("Truncated PF checkpoint header");
    for(double &v:c.wallP) if(!(f>>v)) throw std::runtime_error("Truncated PF wall momentum");
    for(double &v:c.wallL) if(!(f>>v)) throw std::runtime_error("Truncated PF wall angular momentum");
    if(version>=2) {
        for(double &v:c.referenceP) if(!(f>>v)) throw std::runtime_error("Truncated conserved PF momentum");
        for(double &v:c.referenceL) if(!(f>>v)) throw std::runtime_error("Truncated conserved PF angular momentum");
        size_t n=0;
        if(!(f>>n)||n>10000) throw std::runtime_error("Invalid DEM sphere count");
        c.spheres.resize(n);
        for(auto &b:c.spheres) {
            for(double &v:b.center) if(!(f>>v)) throw std::runtime_error("Truncated DEM center");
            for(double &v:b.velocity) if(!(f>>v)) throw std::runtime_error("Truncated DEM velocity");
            for(double &v:b.angular_momentum) if(!(f>>v)) throw std::runtime_error("Truncated DEM spin");
            for(double &v:b.orientation) if(!(f>>v)) throw std::runtime_error("Truncated DEM quaternion");
            int fixed=-1;
            if(!(f>>b.radius>>b.mass>>fixed)||(fixed!=0&&fixed!=1))
                throw std::runtime_error("Invalid/truncated DEM material or flags");
            b.fixed=fixed!=0;
        }
    }
    if(c.nf>100000000 || c.ns>100000000 || c.nf+c.ns>100000000)
        throw std::runtime_error("PF checkpoint particle count out of range");
    c.state.resize(c.nf+c.ns);
    for(auto &a:c.state) {
        int solid=-1;
        for(double &x:a.x) if(!(f>>x)) throw std::runtime_error("Truncated PF positions");
        for(double &v:a.v) if(!(f>>v)) throw std::runtime_error("Truncated PF velocities");
        if(!(f>>a.rho>>a.mass>>solid) || (solid!=0&&solid!=1))
            throw std::runtime_error("Invalid/truncated PF sample");
        a.solid=solid==1;
    }
    std::string extra;
    if(f>>extra) throw std::runtime_error("Unexpected data after PF checkpoint");
    validate_checkpoint(c);
    return c;
}
inline void write_atomic(const std::string &filename,const Checkpoint &c) {
    const std::string tmp=filename+".tmp";
    {
        std::ofstream f(tmp,std::ios::trunc);
        if(!f) throw std::runtime_error("Cannot create PF checkpoint: "+tmp);
        write(f,c);
        f.close();
        if(!f) throw std::runtime_error("Cannot close PF checkpoint: "+tmp);
    }
    if(std::rename(tmp.c_str(),filename.c_str())!=0)
        throw std::runtime_error("Cannot atomically replace PF checkpoint: "+filename);
}
inline Checkpoint read_file(const std::string &filename) {
    std::ifstream f(filename);
    if(!f) throw std::runtime_error("Cannot open PF checkpoint: "+filename);
    return read(f);
}
inline void check_compatibility(const Checkpoint &c,const std::string &kind,
                                const Parameters &p,double radius,
                                const std::vector<Sample> &initial,size_t nf,
                                const std::vector<pf_rigid::Sphere> *expected_spheres=nullptr,
                                const std::vector<pf_rigid::Site> *sites=nullptr,
                                double bridge_gap=0.,double contact_stiffness=0.) {
    auto equal=[](double a,double b) {return std::abs(a-b)<=1e-13*std::max({1.,std::abs(a),std::abs(b)});};
    const auto &q=c.parameters;
    if(c.kind!=kind || !equal(c.radius,radius) ||
       !equal(p.dx,q.dx)||!equal(p.rho0,q.rho0)||!equal(p.sound,q.sound)||
       !equal(p.mu,q.mu)||!equal(p.sff,q.sff)||!equal(p.sfs,q.sfs)||
       !equal(c.bridge_gap,bridge_gap)||!equal(c.contact_stiffness,contact_stiffness)||
       c.nf!=nf||c.state.size()!=initial.size())
        throw std::runtime_error("PF checkpoint incompatible with run parameters or sample counts");
    if(expected_spheres) {
        if(c.spheres.size()!=expected_spheres->size())
            throw std::runtime_error("PF checkpoint DEM body count mismatch");
        for(size_t j=0;j<c.spheres.size();++j) {
            const auto &x=c.spheres[j],&y=(*expected_spheres)[j];
            if(!equal(x.radius,y.radius)||!equal(x.mass,y.mass)||x.fixed!=y.fixed)
                throw std::runtime_error("PF checkpoint DEM radius/mass/fixity mismatch");
        }
    }
    std::vector<int> site_for(initial.size(),-1);
    if(sites) for(size_t k=0;k<sites->size();++k) {
        if((*sites)[k].particle>=initial.size()||site_for[(*sites)[k].particle]>=0)
            throw std::runtime_error("Invalid restart DEM site mapping");
        site_for[(*sites)[k].particle]=static_cast<int>(k);
    }
    if(c.kind=="bridge" && !sites) throw std::runtime_error("Bridge restart requires site mapping");
    for(size_t i=0;i<initial.size();++i) {
        if(!equal(c.state[i].mass,initial[i].mass) || c.state[i].solid!=initial[i].solid)
            throw std::runtime_error("PF checkpoint particle mass/type mismatch");
        if(initial[i].solid) {
            auto expected=initial[i].x;
            if(c.kind=="immersed" || c.kind=="bridge") {
                if(sites) {
                    if(site_for[i]<0) throw std::runtime_error("Missing DEM restart site");
                    const auto &site=(*sites)[static_cast<size_t>(site_for[i])];
                    const auto &b=c.spheres.at(site.body);
                    expected=add(b.center,pf_rigid::rotate(b.orientation,site.local));
                } else {
                    const auto &b=c.spheres.at(0);
                    expected=add(b.center,pf_rigid::rotate(b.orientation,initial[i].x));
                }
            }
            for(size_t k=0;k<3;++k)
                if(!equal(c.state[i].x[k],expected[k]))
                    throw std::runtime_error("PF checkpoint boundary geometry/order mismatch");
        }
    }
}
} // namespace pf_restart
