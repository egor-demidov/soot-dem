#pragma once
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

// Inward normal at the first entry of a fluid-to-solid segment into the union.
// A first entry from outside the union is exposed, including overlapping spheres.
inline Eigen::Vector3d sphereUnionWallNormal(const Eigen::Vector3d& fluid,
 const Eigen::Vector3d& solid,const std::vector<Eigen::Vector3d>& centers,double radius) {
 if(!fluid.allFinite() || !solid.allFinite() || !std::isfinite(radius) || radius<=0)
  throw std::runtime_error("Invalid sphere-union wall query");
 const Eigen::Vector3d d=solid-fluid;
 const double a=d.squaredNorm();
 if(a==0)throw std::runtime_error("Zero wall query segment");
 double first=std::numeric_limits<double>::infinity();std::size_t owner=centers.size();
 for(std::size_t i=0;i<centers.size();++i){
  const Eigen::Vector3d offset=fluid-centers[i];
  const double c=offset.squaredNorm()-radius*radius;
  if(c < -1e-10*radius*radius)throw std::runtime_error("Fluid cell lies inside analytic solid union");
  const double b=offset.dot(d),disc=b*b-a*c;
  if(disc<=0)continue; // A tangent touch does not cross into solid.
  const double t=(-b-std::sqrt(disc))/a;
  if(t>=-1e-12 && t<=1+1e-12 && t<first){first=t;owner=i;}
 }
 if(owner==centers.size())throw std::runtime_error("No exposed sphere crossing on fluid-solid link");
 return (centers[owner]-(fluid+std::max(0.0,first)*d)).normalized();
}
