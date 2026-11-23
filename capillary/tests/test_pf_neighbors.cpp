#include "../src/pf_neighbors.h"
#include "../src/pf_rigid.h"
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace pf_liquid;
using namespace pf_neighbors;
static void require(bool cond,const char* msg) {
    if(!cond) throw std::runtime_error(msg);
}
static Pairs active(const std::vector<Sample> &s,const Pairs &candidates,double H) {
    Pairs out;
    for(const auto &pair:candidates) {
        const auto d=subtract(s[pair.first].x,s[pair.second].x);
        if(dot(d,d)<H*H) out.push_back(pair);
    }
    return out;
}
static void equivalent(const std::vector<Sample>&s,const Pairs&candidates,Parameters p) {
    const auto expected=brute_force(s,std::count_if(s.begin(),s.end(),[](const Sample&x){return !x.solid;}),p.support());
    require(active(s,candidates,p.support())==expected,"Neighbor candidate list missed or duplicated a live pair");
    const auto oracle=evaluate(s,expected,p);
    const auto got=evaluate(s,candidates,p);
    require(oracle.force.size()==got.force.size(),"Wrong rates size");
    auto similar=[](double a,double b) {return std::abs(a-b)<1.e-11*(1+std::abs(a));};
    for(size_t i=0;i<s.size();++i) {
        for(int k=0;k<3;++k)
            require(similar(oracle.force[i][k],got.force[i][k]),"Neighbor force differs from reference");
        require(similar(oracle.density[i],got.density[i]),"Neighbor continuity differs");
        require(similar(oracle.viscous_row[i],got.viscous_row[i]),"Neighbor viscous rate differs");
        require(similar(oracle.conservative_virial[i],got.conservative_virial[i]),"Neighbor virial differs");
    }
    require(similar(oracle.pair_energy,got.pair_energy),"Neighbor capillary energy differs");
    require(similar(oracle.viscous_power,got.viscous_power),"Neighbor dissipation differs");
}
int main() {
    Parameters p;
    p.dx=.08;p.mu=.05;p.sff=.001;p.sfs=.002;
    const double H=p.support();
    std::mt19937_64 rng(84219);
    std::uniform_real_distribution<double> xyz(-.8,.8);
    std::uniform_real_distribution<double> jig(-1.,1.);
    const size_t nf=100,ns=75;
    std::vector<Sample> s;
    for(size_t i=0;i<nf+ns;++i) {
        Sample a;
        a.x={xyz(rng),xyz(rng),xyz(rng)};
        a.v={.005*jig(rng),.005*jig(rng),.005*jig(rng)};
        a.rho=1.+.03*jig(rng);a.mass=.0005;a.solid=i>=nf;
        s.push_back(a);
    }
    // Fixed-time cell linked list and brute-force oracle must agree.
    const auto cellpairs=build_cells(s,nf,H);
    require(cellpairs==brute_force(s,nf,H),"Cell-linked graph differs from exhaustive reference");
    equivalent(s,cellpairs,p);
    Verlet cache(nf,H,.25*H);
    equivalent(s,cache.graph(s),p);
    require(cache.builds()==1&&cache.reuses()==0,"Initial cache build count");
    // Identical state must reuse the graph.
    for(int i=0;i<10;++i) equivalent(s,cache.graph(s),p);
    require(cache.builds()==1&&cache.reuses()==10,"Static graph was rebuilt needlessly");
    // Stress moving fluid and moving solid sites, including negative cells.
    for(int step=0;step<120;++step) {
        for(auto &a:s) {
            for(int axis=0;axis<3;++axis) a.x[axis]+=0.0015*jig(rng);
        }
        equivalent(s,cache.graph(s),p);
    }
    // Force a displacement beyond half the skin to trigger a rebuild.
    for(int axis=0;axis<3;++axis) s[0].x[axis]+=.08;
    equivalent(s,cache.graph(s),p);
    require(cache.builds()>1&&cache.reuses()>cache.builds(),"Verlet cache did not both rebuild and reuse");
    cache.invalidate();
    const auto count=cache.builds();
    equivalent(s,cache.graph(s),p);
    require(cache.builds()==count+1,"Explicit invalidation did not rebuild");
    // Every call with zero skin must rebuild a short-range cell graph.
    Verlet no_skin(nf,H,0.);
    for(int step=0;step<4;++step) equivalent(s,no_skin.graph(s),p);
    require(no_skin.builds()==4 && no_skin.reuses()==0,"Zero skin must always rebuild");
    // Full midpoint time integration must agree, including moving samples.
    auto state_a=s,state_b=s;
    Verlet advancing(nf,H,.25*H);
    auto reference_graph=[&](const auto &x){return brute_force(x,nf,H);};
    auto cached_graph=[&](const auto &x)->const Pairs&{return advancing.graph(x);};
    for(int step=0;step<4;++step) {
        auto r_a=evaluate(state_a,reference_graph(state_a),p);
        auto r_b=evaluate(state_b,cached_graph(state_b),p);
        auto dt=std::min({1.e-4,timestep(state_a,r_a,p),timestep(state_b,r_b,p)});
        advance(state_a,p,dt,reference_graph,&r_a);
        advance(state_b,p,dt,cached_graph,&r_b);
        for(size_t i=0;i<state_a.size();++i) {
            for(int k=0;k<3;++k) {
                require(std::abs(state_a[i].x[k]-state_b[i].x[k])<1.e-12,"Integrated positions differ");
                require(std::abs(state_a[i].v[k]-state_b[i].v[k])<1.e-12,"Integrated velocities differ");
            }
            require(std::abs(state_a[i].rho-state_b[i].rho)<1.e-12,"Integrated density differs");
        }
    }
    // Ordering or NaN errors must be detected, not silently skipped.
    bool caught=false;
    try {auto bad=s;bad[0].solid=true;build_cells(bad,nf,H);}catch(const std::invalid_argument&) {caught=true;}
    require(caught,"Fluid-solid ordering violation was missed");
    std::cout<<"PF neighbor cell/Verlet equivalence, motion, cache reuse and integrator checks passed. "
             <<"builds="<<cache.builds()<<" reuses="<<cache.reuses()<<'\n';
}
