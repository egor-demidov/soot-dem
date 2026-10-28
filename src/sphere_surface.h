//
// Created by egor on 10/27/26.
//

#ifndef SOOT_AFM_SPHERE_SURFACE_H
#define SOOT_AFM_SPHERE_SURFACE_H

#pragma once

#include <Eigen/Dense>
#include <array>
#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>

using SphereTriangle = std::array<Eigen::Vector3d, 3>;

inline std::vector<SphereTriangle> makeSphereSurface(
    const Eigen::Vector3d& center,
    double radius,
    int levels
) {
    using Vec = Eigen::Vector3d;

    const Vec X{1, 0, 0};
    const Vec Y{0, 1, 0};
    const Vec Z{0, 0, 1};

    // Outward-facing triangles of a unit octahedron.
    std::vector<SphereTriangle> mesh = {
        { Z,  X,  Y},
        { Z,  Y, -X},
        { Z, -X, -Y},
        { Z, -Y,  X},
        {-Z,  Y,  X},
        {-Z, -X,  Y},
        {-Z, -Y, -X},
        {-Z,  X, -Y}
    };

    for (int level = 0; level < levels; ++level)
    {
        std::vector<SphereTriangle> refined;
        refined.reserve(4 * mesh.size());

        for (const auto& triangle : mesh)
        {
            const Vec ab =
                (triangle[0] + triangle[1]).normalized();

            const Vec bc =
                (triangle[1] + triangle[2]).normalized();

            const Vec ca =
                (triangle[2] + triangle[0]).normalized();

            refined.push_back({triangle[0], ab, ca});
            refined.push_back({ab, triangle[1], bc});
            refined.push_back({ca, bc, triangle[2]});
            refined.push_back({ab, bc, ca});
        }

        mesh = std::move(refined);
    }

    for (auto& triangle : mesh)
        for (auto& vertex : triangle)
            vertex = center + radius * vertex;

    return mesh;
}

// values contains a signed liquid/gas indicator at each vertex.
//
// Returns:
//   1: a segment was found
//   0: no segment, or a single-point tangency
//  -1: ambiguous intersection along an entire edge
inline int cutSurfaceTriangle(
    const SphereTriangle& triangle,
    const std::array<double, 3>& values,
    double valueTolerance,
    Eigen::Vector3d& a,
    Eigen::Vector3d& b
) {
    using Vec = Eigen::Vector3d;

    const double edgeScale = std::max({
        (triangle[1] - triangle[0]).norm(),
        (triangle[2] - triangle[1]).norm(),
        (triangle[0] - triangle[2]).norm()
    });

    std::vector<Vec> hits;

    auto addPoint = [&](const Vec& point)
    {
        for (const auto& existing : hits)
        {
            if ((point - existing).norm() < 1e-10 * edgeScale)
                return;
        }

        hits.push_back(point);
    };

    for (int i = 0; i < 3; ++i)
    {
        const int j = (i + 1) % 3;

        const double fi = values[i];
        const double fj = values[j];

        if (std::abs(fi) <= valueTolerance &&
            std::abs(fj) <= valueTolerance)
        {
            return -1;
        }

        if (std::abs(fi) <= valueTolerance)
            addPoint(triangle[i]);

        if (std::abs(fj) <= valueTolerance)
            addPoint(triangle[j]);

        if ((fi < -valueTolerance && fj > valueTolerance) ||
            (fi > valueTolerance && fj < -valueTolerance))
        {
            const double fraction = fi / (fi - fj);

            addPoint(
                triangle[i]
                + fraction * (triangle[j] - triangle[i])
            );
        }
    }

    if (hits.size() < 2)
        return 0;

    if (hits.size() != 2)
        return -1;

    a = hits[0];
    b = hits[1];

    return 1;
}

#endif //SOOT_AFM_SPHERE_SURFACE_H
