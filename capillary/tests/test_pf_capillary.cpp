#include "../src/pf_capillary.h"
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
using namespace pf_capillary;
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
bool close(double a,double b,double tolerance=1e-10) {
    return std::abs(a-b)<=tolerance*std::max({1.,std::abs(a),std::abs(b)});
}
double norm(Vector a) { return std::sqrt(dot(a,a)); }
Vector rotate(Vector v,double angle) {
    return {std::cos(angle)*v[0]-std::sin(angle)*v[1],
            std::sin(angle)*v[0]+std::cos(angle)*v[1],v[2]};
}
struct State { double energy=0; std::vector<Vector> force; };
State configuration(const std::vector<Vector> &x) {
    // Last two samples belong to one rigid solid. No solid-solid capillary force.
    State out;out.force.resize(x.size(),Vector{0,0,0});
    for (size_t i=0;i<x.size();i++) for (size_t j=i+1;j<x.size();j++) {
        if (i>=3 && j>=3) continue;
        const auto kind=j>=3?PairKind::fluid_solid:PairKind::fluid_fluid;
        auto p=evaluate(x[i],x[j],1.,j>=3?.006:.01,kind);
        out.energy+=p.energy;accumulate(p,out.force[i],out.force[j]);
    }
    return out;
}
int main() {
    try {
        for (auto kind : {PairKind::fluid_fluid,PairKind::fluid_solid}) {
            require(profile(0,kind)==1,"profile origin");
            require(profile(1,kind)==0 && profile(1.5,kind)==0,"compact support");
            require(profile(.05,kind)>0 && profile(.8,kind)<0,"repulsive/attractive signs");
            require(std::abs(profile(1-1e-7,kind))<1e-11,"smooth force cutoff");
            for (double r : {.11,.32,.71,.999}) {
                const auto pair=evaluate({r,0,0},{0,0,0},1.,.03,kind);
                const double eps=1e-6;
                const double ep=evaluate({r+eps,0,0},{0,0,0},1.,.03,kind).energy;
                const double em=evaluate({r-eps,0,0},{0,0,0},1.,.03,kind).energy;
                require(close(pair.force_on_i[0],-(ep-em)/(2*eps)),"force/energy gradient");
            }
        }
        require(std::abs(profile(.5,PairKind::fluid_solid))<1e-14,"published fs zero");
        const double root=(2+std::sqrt(86.))/41.;
        require(std::abs(profile(root,PairKind::fluid_fluid))<1e-14,"published ff zero");
        std::vector<Vector> x{{.02,.04,.1},{.28,-.1,.2},{-.15,.33,.25},{.1,.1,-.15},{-.3,.15,-.1}};
        const auto original=configuration(x);
        Vector total{0,0,0},moment{0,0,0};
        for (size_t i=0;i<x.size();i++) {
            total=add(total,original.force[i]);moment=add(moment,cross(x[i],original.force[i]));
            for (int a=0;a<3;a++) {
                auto xp=x,xm=x;xp[i][a]+=1e-6;xm[i][a]-=1e-6;
                require(close(original.force[i][a],-(configuration(xp).energy-configuration(xm).energy)/2e-6),"network gradient");
            }
        }
        require(norm(total)<1e-16 && norm(moment)<1e-16,"network force/torque reciprocity");
        auto moved=x;
        for (auto &v:moved) v=add(rotate(v,.61),Vector{2.,-1.,.3});
        const auto transformed=configuration(moved);
        require(close(transformed.energy,original.energy,1e-14),"energy frame invariance");
        for(size_t i=0;i<x.size();i++)
            require(norm(subtract(transformed.force[i],rotate(original.force[i],.61)))<1e-15,"force covariance");
        const Vector center{.2,-.2,.1};RigidLoad body;
        for(size_t i=3;i<x.size();i++) body.add_sample(x[i],original.force[i],center);
        require(norm(body.torque)>1e-7,"nontrivial rigid torque test");
        const double eps=1e-6;
        for(int axis=0;axis<3;axis++) {
            auto xp=x,xm=x;Vector direction{0,0,0};direction[axis]=1;
            for(size_t i=3;i<x.size();i++) {xp[i]=add(x[i],scale(direction,eps));xm[i]=subtract(x[i],scale(direction,eps));}
            require(close(body.force[axis],-(configuration(xp).energy-configuration(xm).energy)/(2*eps)),"rigid force virtual work");
            xp=x;xm=x;
            for(size_t i=3;i<x.size();i++) {
                auto displacement=cross(direction,subtract(x[i],center));
                xp[i]=add(x[i],scale(displacement,eps));xm[i]=subtract(x[i],scale(displacement,eps));
            }
            require(close(body.torque[axis],-(configuration(xp).energy-configuration(xm).energy)/(2*eps)),"rigid torque virtual work");
        }
        // Pair impulses preserve total linear and angular momentum with unequal masses.
        const auto pair=evaluate(x[0],x[1],1.,.01,PairKind::fluid_fluid);
        const Vector p_i=scale(pair.force_on_i,.002),p_j=scale(pair.force_on_i,-.002);
        require(norm(add(p_i,p_j))==0,"pair impulse momentum");
        require(norm(add(cross(x[0],p_i),cross(x[1],p_j)))<1e-18,"pair impulse angular momentum");
        // Explicit mass division belongs in the caller, never in an action/reaction force.
        require(close(norm(scale(p_i,1./.001))/norm(scale(p_j,1./1.)),1000),"unequal mass velocity kicks");
        const auto scaled=evaluate(scale(x[0],3),scale(x[1],3),3.,.01,PairKind::fluid_fluid);
        require(norm(subtract(scaled.force_on_i,scale(pair.force_on_i,3)))<1e-15,"force dimensional scaling");
        require(close(scaled.energy,9*pair.energy,1e-14),"energy dimensional scaling");
        bool rejected=false;
        try {evaluate({0,0,0},{0,0,0},1.,.01,PairKind::fluid_fluid);}
        catch(const std::domain_error &) {rejected=true;}
        require(rejected,"coincident samples must fail explicitly");
        rejected=false;
        try {evaluate(x[0],x[1],0.,.01,PairKind::fluid_fluid);}
        catch(const std::invalid_argument &) {rejected=true;}
        require(rejected,"invalid support");
        std::cout<<"PF-SPH pair law: profile, energy, force, torque, rigid virtual work and scaling checks passed.\n";
        return 0;
    } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
