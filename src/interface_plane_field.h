//
// Created by egor on 10/27/26.
//

#ifndef SOOT_AFM_INTERFACE_PLANE_FIELD_H
#define SOOT_AFM_INTERFACE_PLANE_FIELD_H

#include <Eigen/Dense>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

struct InterfaceFieldSample {
    double value;
    Eigen::Vector3d gradient;
    std::size_t sampleCount;
};

template<class PlaneContainer>
std::optional<InterfaceFieldSample> sampleInterfaceField(
    const Eigen::Vector3d& x,
    const PlaneContainer& planes,
    double supportRadius
) {
    using Vec = Eigen::Vector3d;

    if (!x.allFinite() ||
        !std::isfinite(supportRadius) ||
        supportRadius <= 0.0)
    {
        throw std::invalid_argument(
            "Invalid interface-field query"
        );
    }

    const double h2 = supportRadius * supportRadius;

    double weightSum = 0.0;
    double numerator = 0.0;

    Vec weightGradient = Vec::Zero();
    Vec numeratorGradient = Vec::Zero();

    std::size_t count = 0;

    for (const auto& plane : planes)
    {
        const Vec offset = x - plane.point;
        const double distanceSquared = offset.squaredNorm();

        if (distanceSquared >= h2)
            continue;

        const double s =
            std::sqrt(distanceSquared) / supportRadius;

        const double oneMinus = 1.0 - s;

        // Smooth weight with compact support.
        const double w =
            std::pow(oneMinus, 4) * (4.0 * s + 1.0);

        const Vec dw =
            (-20.0 * std::pow(oneMinus, 3) / h2) * offset;

        const double d = plane.normal.dot(offset);

        weightSum += w;
        numerator += w * d;

        weightGradient += dw;
        numeratorGradient += dw * d + w * plane.normal;

        ++count;
    }

    // No phase is invented where support is insufficient.
    if (count < 3 || weightSum < 1e-12)
        return std::nullopt;

    const double value = numerator / weightSum;

    const Vec gradient =
        (numeratorGradient - value * weightGradient)
        / weightSum;

    if (!std::isfinite(value) || !gradient.allFinite())
        return std::nullopt;

    return InterfaceFieldSample{
        value,
        gradient,
        count
    };
}

#endif //SOOT_AFM_INTERFACE_PLANE_FIELD_H
