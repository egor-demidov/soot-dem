#pragma once
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

// Same distance and nonbranching test as the all-pairs validator. Sort along
// the widest axis and examine only endpoints within the matching tolerance.
// No snapping, quantization, synthetic links, or relaxed closure criterion.
inline std::size_t countBadSurfaceEndpoints(
    const std::vector<std::array<Eigen::Vector3d,2>>& segments, double tolerance) {
    if (!std::isfinite(tolerance) || tolerance<=0)
        throw std::invalid_argument("Invalid contact-line endpoint tolerance");
    if(segments.empty())return 0;
    Eigen::Vector3d lower=segments[0][0],upper=lower;
    for(const auto& segment:segments)for(const auto& p:segment){
        if(!p.allFinite())throw std::invalid_argument("Nonfinite contact-line endpoint");
        lower=lower.cwiseMin(p);upper=upper.cwiseMax(p);
    }
    Eigen::Index axis;
    (upper-lower).maxCoeff(&axis);
    std::vector<std::size_t> order(2*segments.size()),matches(order.size(),0);
    std::iota(order.begin(),order.end(),std::size_t(0));
    auto point=[&](std::size_t i)->const Eigen::Vector3d&{return segments[i/2][i%2];};
    std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){
        return point(a)[axis]<point(b)[axis] || (point(a)[axis]==point(b)[axis] && a<b);
    });
    const double squared=tolerance*tolerance;
    for(std::size_t i=0;i<order.size();++i)for(std::size_t j=i+1;j<order.size();++j){
        const auto a=order[i],b=order[j];
        if(point(b)[axis]-point(a)[axis]>tolerance)break;
        if(a/2!=b/2 && (point(a)-point(b)).squaredNorm()<squared){++matches[a];++matches[b];}
    }
    return std::count_if(matches.begin(),matches.end(),[](std::size_t n){return n!=1;});
}
