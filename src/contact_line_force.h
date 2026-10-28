//
// Created by egor on 10/27/26.
//

#ifndef SOOT_AFM_CONTACT_LINE_FORCE_H
#define SOOT_AFM_CONTACT_LINE_FORCE_H

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <array>
#include <vector>

struct ContactLineLoad {
    Eigen::Vector3d force = Eigen::Vector3d::Zero();
    Eigen::Vector3d torque = Eigen::Vector3d::Zero();
};

inline ContactLineLoad contactLineSegmentLoad(
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& solidNormal,
    const Eigen::Vector3d& interfaceNormal,
    const Eigen::Vector3d& particleCenter,
    double sigma
) {
    if (!a.allFinite() ||
        !b.allFinite() ||
        !particleCenter.allFinite() ||
        !solidNormal.allFinite() ||
        !interfaceNormal.allFinite() ||
        !std::isfinite(sigma) ||
        sigma < 0.0 ||
        solidNormal.norm() == 0.0 ||
        interfaceNormal.norm() == 0.0)
    {
        throw std::invalid_argument(
            "Invalid contact-line segment"
        );
    }

    const double length = (b - a).norm();

    if (length == 0.0)
        return {};

    const Eigen::Vector3d ns = solidNormal.normalized();
    const Eigen::Vector3d ni = interfaceNormal.normalized();

    // Tangent to the liquid-gas interface, perpendicular to
    // the contact line, pointing away from the solid.
    Eigen::Vector3d direction =
        ns - ns.dot(ni) * ni;

    if (direction.norm() < 1e-10) {
        throw std::runtime_error(
            "Degenerate contact-line geometry"
        );
    }

    direction.normalize();

    ContactLineLoad load;

    load.force = sigma * length * direction;

    const Eigen::Vector3d midpoint = 0.5 * (a + b);

    load.torque =
        (midpoint - particleCenter).cross(load.force);

    return load;
}

// Coordinates are local to a unit cell: [0,1]^3.
//
// axis: 0=x, 1=y, 2=z
// side: -1 selects coordinate 0, +1 selects coordinate 1
//
// Returns:
//   1: one segment
//   0: no segment, or a single-point tangency
//  -1: degenerate intersection along an edge or entire face
inline int interfacePlaneFaceSegment(
    const Eigen::Vector3d& normal,
    const Eigen::Vector3d& planePoint,
    int axis,
    int side,
    Eigen::Vector3d& a,
    Eigen::Vector3d& b
) {
    using Vec = Eigen::Vector3d;

    constexpr double eps = 1e-10;

    const int u = (axis + 1) % 3;
    const int v = (axis + 2) % 3;

    const int uv[4][2] = {
        {0, 0}, {1, 0}, {1, 1}, {0, 1}
    };

    std::array<Vec, 4> corners;

    for (int i = 0; i < 4; ++i) {
        corners[i] = Vec::Zero();
        corners[i][axis] = side > 0 ? 1.0 : 0.0;
        corners[i][u] = uv[i][0];
        corners[i][v] = uv[i][1];
    }

    std::vector<Vec> hits;

    auto addPoint = [&](const Vec& point) {
        for (const auto& existing : hits) {
            if ((point - existing).norm() < eps)
                return;
        }

        hits.push_back(point);
    };

    for (int i = 0; i < 4; ++i) {
        const Vec& p = corners[i];
        const Vec& q = corners[(i + 1) % 4];

        const double dp = normal.dot(p - planePoint);
        const double dq = normal.dot(q - planePoint);

        if (std::abs(dp) <= eps && std::abs(dq) <= eps)
            return -1;

        if (std::abs(dp) <= eps)
            addPoint(p);

        if (std::abs(dq) <= eps)
            addPoint(q);

        if ((dp < -eps && dq > eps) ||
            (dp > eps && dq < -eps))
        {
            addPoint(p + (dp / (dp - dq)) * (q - p));
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

#endif //SOOT_AFM_CONTACT_LINE_FORCE_H
