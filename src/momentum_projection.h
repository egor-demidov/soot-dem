#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// Numerical conservation projection for the single-particle prototype.
// This is not a replacement for a consistent, converged capillary-stress model.
namespace capillary_coupling {

struct MomentumState {
    Eigen::Vector3d linear = Eigen::Vector3d::Zero();
    Eigen::Vector3d angular = Eigen::Vector3d::Zero();
};

struct ProjectionCell {
    Eigen::Vector3d position;
    double mass; // Resolved liquid mass phi*rho, in lattice units.
};

struct MomentumProjection {
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    Eigen::Vector3d translation = Eigen::Vector3d::Zero();
    Eigen::Vector3d rotation = Eigen::Vector3d::Zero();

    Eigen::Vector3d velocityIncrement(const Eigen::Vector3d &position) const {
        return translation + rotation.cross(position - centroid);
    }
};

// Minimize sum(m*|du|^2)/2 subject to the requested total linear/angular
// impulse. All positions and angular impulses use the SAME reference origin.
// A rigid velocity increment about the liquid centroid is the minimizer.
inline MomentumProjection makeMomentumProjection(const std::vector<ProjectionCell> &cells,
                                                 const Eigen::Vector3d &linearImpulse,
                                                 const Eigen::Vector3d &angularImpulse,
                                                 const Eigen::Vector3d &origin) {
    if (!linearImpulse.allFinite() || !angularImpulse.allFinite() || !origin.allFinite())
        throw std::invalid_argument("Non-finite projection impulse or origin");

    MomentumProjection result;
    if (linearImpulse.isZero(0.0) && angularImpulse.isZero(0.0))
        return result;

    double mass = 0.0;
    for (const auto &cell : cells) {
        if (!cell.position.allFinite() || !std::isfinite(cell.mass) || cell.mass <= 0.0)
            throw std::invalid_argument("Projection cells require finite positions and positive mass");
        mass += cell.mass;
        // Work relative to origin to avoid loss of accuracy at large coordinates.
        result.centroid += cell.mass * (cell.position - origin);
    }
    if (!std::isfinite(mass) || mass <= 0.0)
        throw std::runtime_error("No liquid cells available for momentum projection");

    const Eigen::Vector3d relativeCentroid = result.centroid / mass;
    result.centroid = origin + relativeCentroid;
    result.translation = linearImpulse / mass;

    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();
    for (const auto &cell : cells) {
        const Eigen::Vector3d r = (cell.position - origin) - relativeCentroid;
        inertia += cell.mass * (r.squaredNorm() * Eigen::Matrix3d::Identity() - r * r.transpose());
    }
    const Eigen::Vector3d spinImpulse = angularImpulse - relativeCentroid.cross(linearImpulse);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(inertia);
    if (solver.info() != Eigen::Success || !solver.eigenvalues().allFinite() ||
        solver.eigenvalues().minCoeff() <= 1e-12 * solver.eigenvalues().maxCoeff() ||
        solver.eigenvalues().minCoeff() <= 0.0)
        throw std::runtime_error("Insufficient spatial support for angular-momentum projection");

    result.rotation = solver.eigenvectors() * (solver.eigenvalues().cwiseInverse().asDiagonal() *
                                               (solver.eigenvectors().transpose() * spinImpulse));
    if (!result.translation.allFinite() || !result.rotation.allFinite())
        throw std::runtime_error("Non-finite momentum projection");
    return result;
}

// For D3Q19 (cs^2=1/3), this odd PDF increment changes j by rho*du,
// preserves density, and leaves the raw second moment unchanged.
inline double momentumPdfIncrement(double weight, double density, const Eigen::Vector3d &direction,
                                   const Eigen::Vector3d &velocityIncrement) {
    return 3.0 * weight * density * direction.dot(velocityIncrement);
}

struct InterfaceMassPlan {
    std::vector<double> delta;
    double availablePositiveMass = 0.0;
    std::size_t eligibleCount = 0;
};

// Existing interface states retain the native signed fill; only freshly
// initialized cells must start with a physical fill fraction.
inline bool validRemapFill(double fill, bool initializePdfs) {
    return std::isfinite(fill) && (!initializePdfs || (fill >= 0.0 && fill <= 1.0));
}

// Allocate the signed covered-minus-refill mass without clipping native states.
// Negative/depleted interface cells may receive mass but cannot donate mass.
inline InterfaceMassPlan makeInterfaceMassPlan(const std::vector<double> &signedMass,
                                               const std::vector<bool> &eligible, double requiredDelta) {
    if (signedMass.size() != eligible.size() || !std::isfinite(requiredDelta))
        throw std::invalid_argument("Invalid interface mass-plan inputs");
    InterfaceMassPlan plan;
    plan.delta.assign(signedMass.size(), 0.0);
    std::size_t correctionIndex = signedMass.size();
    double largestCapacity = -1.0;
    for (std::size_t i = 0; i < signedMass.size(); ++i) {
        if (!std::isfinite(signedMass[i]))
            throw std::invalid_argument("Non-finite interface mass in remap plan");
        if (!eligible[i])
            continue;
        ++plan.eligibleCount;
        const double capacity = std::max(0.0, signedMass[i]);
        plan.availablePositiveMass += capacity;
        if (capacity > largestCapacity) {
            largestCapacity = capacity;
            correctionIndex = i;
        }
    }
    if (!std::isfinite(plan.availablePositiveMass))
        throw std::runtime_error("Non-finite total donor capacity");
    if (requiredDelta == 0.0)
        return plan;
    if (requiredDelta > 0.0 && plan.eligibleCount == 0)
        throw std::runtime_error("Remap has displaced mass but no surviving interface recipients");
    if (requiredDelta < 0.0 &&
        (plan.availablePositiveMass <= 0.0 || -requiredDelta > plan.availablePositiveMass))
        throw std::runtime_error("Insufficient remap donor mass: required=" + std::to_string(-requiredDelta) +
                                 " available=" + std::to_string(plan.availablePositiveMass));

    double allocated = 0.0;
    for (std::size_t i = 0; i < signedMass.size(); ++i) {
        if (!eligible[i])
            continue;
        plan.delta[i] = requiredDelta > 0.0
                            ? requiredDelta / double(plan.eligibleCount)
                            : requiredDelta * std::max(0.0, signedMass[i]) / plan.availablePositiveMass;
        allocated += plan.delta[i];
    }
    // Keep all of the requested mass, including summation roundoff.
    plan.delta.at(correctionIndex) += requiredDelta - allocated;
    return plan;
}

} // namespace capillary_coupling
