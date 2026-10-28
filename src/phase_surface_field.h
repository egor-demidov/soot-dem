#pragma once
#include "interface_plane_field.h"
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>

// Cell-centered, wall-extended fill. Positive value points toward gas.
// The callback reads the geometry-stage field, not solid bookkeeping fill.
struct PhaseSurfaceField {
    std::array<int,3> size;
    std::function<double(int,int,int)> fill;
};
inline std::optional<InterfaceFieldSample> sampleInterfaceField(
    const Eigen::Vector3d& point,const PhaseSurfaceField& field,double) {
    if(!point.allFinite() || !field.fill)throw std::runtime_error("Invalid surface phase query");
    const Eigen::Vector3d local=point-Eigen::Vector3d::Constant(.5);
    std::array<int,3> base;Eigen::Vector3d u;
    for(int k=0;k<3;++k){if(local[k]<0 || local[k]>=field.size[k]-1)return std::nullopt;
        base[k]=static_cast<int>(std::floor(local[k]));u[k]=local[k]-base[k];
        if(base[k]<0 || base[k]+1>=field.size[k])return std::nullopt;}
    double phi=0;Eigen::Vector3d gradient=Eigen::Vector3d::Zero();
    for(int z=0;z<2;++z)for(int y=0;y<2;++y)for(int x=0;x<2;++x){
        const int bit[3]={x,y,z};double w[3];
        for(int k=0;k<3;++k)w[k]=bit[k]?u[k]:1-u[k];
        const double value=field.fill(base[0]+x,base[1]+y,base[2]+z);
        if(!std::isfinite(value))throw std::runtime_error("Nonfinite surface phase sample");
        phi+=w[0]*w[1]*w[2]*value;
        for(int k=0;k<3;++k)gradient[k]-=(bit[k]?1.:-1.)*w[(k+1)%3]*w[(k+2)%3]*value;
    }
    return InterfaceFieldSample{.5-phi,gradient,8};
}
