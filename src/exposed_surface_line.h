#pragma once
#include "surface_contact_line.h"
#include <algorithm>
#include <cmath>
#include <utility>

// Validate the complete interpolated sphere contour first. Then intersect each
// segment with the complement of other analytic spheres. New endpoints are
// allowed only where the segment crosses an occluding solid boundary.
inline SurfaceLineResult exposeSurfaceLine(const SurfaceLineResult& full,
 const std::vector<Eigen::Vector3d>& centers,std::size_t owner,double radius,double sigma) {
    if(full.badEndpoints || full.invalidSegments || full.unsupportedTriangles || full.degenerateTriangles)
        throw std::runtime_error("Incomplete phase contour before solid-union clipping");
    if(owner>=centers.size() || full.segments.size()!=full.solidNormals.size() ||
       full.segments.size()!=full.interfaceNormals.size())throw std::runtime_error("Invalid exposed line metadata");
    SurfaceLineResult out;
    for(std::size_t i=0;i<full.segments.size();++i){
        const auto a=full.segments[i][0];const Eigen::Vector3d d=full.segments[i][1]-a;
        const double aa=d.squaredNorm();
        std::vector<std::array<double,2>> visible{{0.,1.}};
        for(std::size_t j=0;j<centers.size() && !visible.empty();++j){
            if(j==owner)continue;
            const Eigen::Vector3d q=a-centers[j];const double b=q.dot(d),c=q.squaredNorm()-radius*radius;
            const double disc=b*b-aa*c;if(disc<=0)continue;
            const double lo=std::max(0.,(-b-std::sqrt(disc))/aa),hi=std::min(1.,(-b+std::sqrt(disc))/aa);
            if(hi<=lo)continue;
            std::vector<std::array<double,2>> next;
            for(const auto& v:visible){
                if(hi<=v[0] || lo>=v[1]){next.push_back(v);continue;}
                if(lo>v[0]){next.push_back({v[0],lo});++out.occlusionCuts;}
                if(hi<v[1]){next.push_back({hi,v[1]});++out.occlusionCuts;}
            }
            visible=std::move(next);
        }
        for(const auto& interval:visible){
            const Eigen::Vector3d p=a+interval[0]*d,q=a+interval[1]*d;
            if((q-p).norm()<1e-14)continue;
            auto load=contactLineSegmentLoad(p,q,full.solidNormals[i],full.interfaceNormals[i],centers[owner],sigma);
            out.load.force+=load.force;out.load.torque+=load.torque;out.length+=(q-p).norm();
            out.segments.push_back({p,q});out.solidNormals.push_back(full.solidNormals[i]);out.interfaceNormals.push_back(full.interfaceNormals[i]);
        }
    }
    return out;
}
