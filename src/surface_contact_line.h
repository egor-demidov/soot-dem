//
// Created by egor on 10/27/26.
//

#ifndef SOOT_AFM_SURFACE_CONTACT_LINE_H
#define SOOT_AFM_SURFACE_CONTACT_LINE_H

#include <optional>
#include <stdexcept>

#include "contact_line_force.h"
#include "sphere_surface.h"
#include "interface_plane_field.h"

struct SurfaceLineResult {
    ContactLineLoad load;
    double length = 0.0;

    std::vector<std::array<Eigen::Vector3d, 2>> segments;

    std::size_t unsupportedTriangles = 0;
    std::size_t degenerateTriangles = 0;
    std::size_t invalidSegments = 0;
    std::size_t badEndpoints = 0;
};

template<class PlaneContainer>
SurfaceLineResult reconstructSurfaceLine(
    const std::vector<SphereTriangle>& mesh,
    const PlaneContainer& planes,
    const Eigen::Vector3d& particleCenter,
    double sigma,
    double supportRadius,
    std::optional<double> wallContactAngleDegrees = std::nullopt
) {
    using Vec = Eigen::Vector3d;

    SurfaceLineResult result;

    const double tolerance = 1e-10 * supportRadius;

    for (const auto& triangle : mesh)
    {
        std::array<double, 3> values;
        bool supported = true;

        for (int i = 0; i < 3; ++i)
        {
            const auto sample = sampleInterfaceField(
                triangle[i], planes, supportRadius
            );

            if (!sample) {
                supported = false;
                break;
            }

            values[i] = sample->value;
        }

        if (!supported) {
            ++result.unsupportedTriangles;
            continue;
        }

        Vec a, b;

        const int hit = cutSurfaceTriangle(
            triangle, values, tolerance, a, b
        );

        if (hit < 0) {
            ++result.degenerateTriangles;
            continue;
        }

        if (hit == 0)
            continue;

        const Vec midpoint = 0.5 * (a + b);

        const auto sample = sampleInterfaceField(
            midpoint, planes, supportRadius
        );

        const Vec faceNormal =
            (triangle[1] - triangle[0])
            .cross(triangle[2] - triangle[0]);

        const double length = (b - a).norm();

        if (!sample ||
            length == 0.0 ||
            faceNormal.squaredNorm() == 0.0)
        {
            ++result.invalidSegments;
            continue;
        }

        const Vec tangent = (b - a) / length;

        Vec ni =
            sample->gradient
            - sample->gradient.dot(tangent) * tangent;

        if (ni.norm() < 1e-10) {
            ++result.invalidSegments;
            continue;
        }

        ni.normalize();

        const Vec ns = faceNormal.normalized();

        if (wallContactAngleDegrees)
        {
            const double angle = *wallContactAngleDegrees;

            if (!std::isfinite(angle) ||
                angle <= 0.0 || angle >= 180.0)
            {
                throw std::invalid_argument(
                    "Wall contact angle must be between 0 and 180 degrees"
                );
            }

            // Tangent to the solid, pointing toward the dry side.
            Vec dryDirection = ni - ni.dot(ns) * ns;

            if (dryDirection.norm() < 1e-10)
            {
                ++result.invalidSegments;
                continue;
            }

            dryDirection.normalize();

            const double theta =
                angle * std::acos(-1.0) / 180.0;

            // Liquid-to-gas interface normal at the wall.
            ni = std::cos(theta) * ns
               + std::sin(theta) * dryDirection;
        }

        if ((ns - ns.dot(ni) * ni).norm() < 1e-10) {
            ++result.invalidSegments;
            continue;
        }

        const auto load = contactLineSegmentLoad(
            a, b, ns, ni, particleCenter, sigma
        );

        result.load.force += load.force;
        result.load.torque += load.torque;

        result.length += length;
        result.segments.push_back({a, b});
    }

    // A closed, non-branching contour has exactly one other
    // segment endpoint matching each endpoint.
    const double endpointToleranceSquared =
        std::pow(1e-8 * supportRadius, 2);

    for (std::size_t i = 0; i < result.segments.size(); ++i)
    {
        for (int end = 0; end < 2; ++end)
        {
            std::size_t matches = 0;

            for (std::size_t j = 0;
                 j < result.segments.size(); ++j)
            {
                if (i == j)
                    continue;

                for (int other = 0; other < 2; ++other)
                {
                    const double distanceSquared =
                        (result.segments[i][end]
                         - result.segments[j][other])
                        .squaredNorm();

                    if (distanceSquared <
                        endpointToleranceSquared)
                    {
                        ++matches;
                    }
                }
            }

            if (matches != 1)
                ++result.badEndpoints;
        }
    }

    return result;
}

#endif //SOOT_AFM_SURFACE_CONTACT_LINE_H
