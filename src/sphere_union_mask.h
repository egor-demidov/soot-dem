#pragma once
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

// Exact sphere-union membership at owned lattice cell centers. Rasterize only
// each sphere's conservative cell bounding box; overlaps are an ordinary union.
// This is rebuilt for every proposed geometry, never reused across remaps.
class SphereUnionMask {
    std::array<int,3> size_;
    std::vector<unsigned char> occupied_;
    std::size_t index(int x, int y, int z) const {
        return (std::size_t(z)*size_[1]+y)*size_[0]+x;
    }
public:
    SphereUnionMask(std::array<int,3> size, const Eigen::Vector3d& firstCenter,
                    const std::vector<Eigen::Vector3d>& centers, double radius) : size_(size) {
        if (!firstCenter.allFinite() || !std::isfinite(radius) || radius <= 0 ||
            size[0]<=0 || size[1]<=0 || size[2]<=0)
            throw std::invalid_argument("Invalid sphere-union lattice");
        occupied_.assign(std::size_t(size[0])*size[1]*size[2], 0);
        const double radiusSquared=radius*radius;
        for (const auto& center:centers) {
            if (!center.allFinite()) throw std::invalid_argument("Invalid sphere-union center");
            std::array<int,3> lo,hi;
            for(int a=0;a<3;++a) {
                // Clamp before integer conversion, including far-away spheres.
                lo[a]=int(std::clamp(std::floor(center[a]-firstCenter[a]-radius)-1., 0., double(size[a])));
                hi[a]=int(std::clamp(std::ceil(center[a]-firstCenter[a]+radius)+1., -1., double(size[a]-1)));
            }
            for(int z=lo[2];z<=hi[2];++z)for(int y=lo[1];y<=hi[1];++y)for(int x=lo[0];x<=hi[0];++x) {
                const Eigen::Vector3d position=firstCenter+Eigen::Vector3d(x,y,z);
                if ((position-center).squaredNorm()<=radiusSquared) occupied_[index(x,y,z)]=1;
            }
        }
    }
    bool operator()(int x,int y,int z) const { return occupied_[index(x,y,z)]!=0; }
};
