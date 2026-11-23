#pragma once
#include <array>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <istream>
#include <ostream>
namespace pf_dem {
struct SootSettings {
    bool enabled=false;
    // Contact: normal; sliding(k,gamma,mu,phi); rolling; torsion.
    // Bond: normal, sliding, rolling, torsion (k,gamma each).
    // Then d_crit,A,h0,max_dt,normal_yield,shear_yield (zero yield disables).
    std::array<double,28> v{100,.02,20,0,.5,.8,5,0,.1,.8,5,0,.1,.8,
                            100,.02,20,0,5,0,5,0,.02,0,.01,.0001,0,0};
    int neck_mode=1; // 0 none, 1 upstream automatic initial graph, 2 explicit graph
    std::vector<std::pair<size_t,size_t>> edges;
    void validate(size_t n) const {
        for(double x:v)if(!(x>=0&&std::isfinite(x)))throw std::runtime_error("Invalid soot-dem coefficient");
        if(!(v[24]>0&&v[25]>0)||neck_mode<0||neck_mode>2||n<1||n>2000)
            throw std::runtime_error("Invalid soot-dem h0/max_dt/neck mode/body count (limit 2000)");
        for(size_t i:{size_t(2),size_t(6),size_t(10)})
            if(v[i]==0&&(v[i+1]>0||v[i+2]>0))throw std::runtime_error("Nonzero friction/damping requires positive shear stiffness");
        auto ordered=edges;std::sort(ordered.begin(),ordered.end());
        if(std::adjacent_find(ordered.begin(),ordered.end())!=ordered.end())throw std::runtime_error("Duplicate initial neck");
        for(auto e:edges)if(e.first>=e.second||e.second>=n)throw std::runtime_error("Invalid initial neck endpoints");
        if(neck_mode!=2&&!edges.empty())throw std::runtime_error("Neck list requires explicit mode");
    }
    void write(std::ostream &f) const {
        for(double x:v)f<<x<<' ';
        f<<'\n'<<neck_mode<<' '<<edges.size()<<'\n';
        for(auto e:edges)f<<e.first<<' '<<e.second<<'\n';
    }
    void read(std::istream &f,size_t n) {
        enabled=true;for(double &x:v)f>>x;size_t count;
        if(!(f>>neck_mode>>count)||count>n*n)throw std::runtime_error("Invalid soot-dem settings record");
        edges.resize(count);for(auto &e:edges)f>>e.first>>e.second;
        if(!f)throw std::runtime_error("Truncated soot-dem settings");
        validate(n);
    }
};
}
