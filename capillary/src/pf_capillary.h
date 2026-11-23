#pragma once
// Gas-free PF-SPH capillary pair law, independent of SPHinXsys.
// Whebell et al., JFM 1019 A37 (2025), DOI:10.1017/jfm.2025.10601,
// equations (2.19)-(2.21). Analytic potential derived by integrating that law.
// s_ff and s_fs require calibration; neither is automatically physical gamma.
#include <array>
#include <cmath>
#include <stdexcept>

namespace pf_capillary {
using Vector = std::array<double, 3>;
inline Vector add(const Vector &a, const Vector &b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}
inline Vector subtract(const Vector &a, const Vector &b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
inline Vector scale(const Vector &a, double s) { return {s*a[0], s*a[1], s*a[2]}; }
inline double dot(const Vector &a, const Vector &b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
inline Vector cross(const Vector &a, const Vector &b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
inline bool finite(const Vector &v) {
    return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}
enum class PairKind { fluid_fluid, fluid_solid };
inline double alpha(PairKind kind) {
    switch (kind) {
    case PairKind::fluid_fluid: return -23.5;
    case PairKind::fluid_solid: return -11.;
    }
    throw std::invalid_argument("Unknown PF-SPH pair kind");
}
inline double profile(double q, PairKind kind) {
    if (!(q >= 0 && std::isfinite(q))) throw std::invalid_argument("Invalid PF-SPH distance");
    const double a = alpha(kind);
    if (q >= 1) return 0.;
    // Factored form of 1+a*q^2-2*(2+a)*q^3+(3+a)*q^4.
    return (1-q)*(1-q)*(1+2*q+(3+a)*q*q);
}
inline double profile_integral_to_cutoff(double q, PairKind kind) {
    if (!(q >= 0 && std::isfinite(q))) throw std::invalid_argument("Invalid PF-SPH distance");
    const double a = alpha(kind);
    if (q >= 1) return 0.;
    const double z = 1-q;
    // Integral_q^1 f(t)dt; this avoids cancellation close to the cutoff.
    return z*z*z*((a+6)/3 - (a+4)*z/2 + (a+3)*z*z/5);
}
struct Interaction {
    Vector force_on_i{0.,0.,0.};
    double energy = 0.;
    double scalar_virial = 0.; // (x_i-x_j).F_ij, not a pressure by itself.
};
inline Interaction evaluate(const Vector &x_i, const Vector &x_j,
                             double support, double strength, PairKind kind) {
    if (!(support > 0 && strength >= 0 && std::isfinite(support) && std::isfinite(strength)) ||
        !finite(x_i) || !finite(x_j))
        throw std::invalid_argument("Invalid PF-SPH coordinates or parameters");
    alpha(kind); // Validate even for an inactive pair.
    const Vector displacement = subtract(x_i,x_j);
    const double distance = std::hypot(displacement[0],displacement[1],displacement[2]);
    if (!std::isfinite(distance)) throw std::invalid_argument("PF-SPH distance overflow");
    if (distance >= support || strength == 0) return {};
    // The published law has nonzero magnitude at r=0 and no defined direction.
    // Self pairs must be excluded by the neighbor iterator. Do not hide overlaps
    // of distinct samples with an arbitrary vector or a force clip.
    if (distance == 0) throw std::domain_error("Coincident active PF-SPH samples");
    const double q = distance/support;
    const double magnitude = support*strength*profile(q,kind);
    Interaction out;
    out.force_on_i = scale(displacement,magnitude/distance);
    out.energy = strength*support*support*profile_integral_to_cutoff(q,kind);
    out.scalar_virial = distance*magnitude;
    if (!finite(out.force_on_i) || !std::isfinite(out.energy) || !std::isfinite(out.scalar_virial))
        throw std::overflow_error("PF-SPH interaction overflow");
    return out;
}
// Process a UNIQUE undirected pair once. The caller owns neighbor traversal and
// chooses ff/fs symmetrically. There is no gas state or solid-solid capillary pair.
inline void accumulate(const Interaction &pair, Vector &force_i, Vector &force_j) {
    if (&force_i == &force_j) throw std::invalid_argument("PF-SPH pair accumulators alias");
    force_i = add(force_i,pair.force_on_i);
    force_j = subtract(force_j,pair.force_on_i);
}
struct RigidLoad {
    Vector force{0.,0.,0.}, torque{0.,0.,0.};
    void add_sample(const Vector &sample_position, const Vector &sample_force,
                    const Vector &body_center) {
        force = add(force,sample_force);
        torque = add(torque,cross(subtract(sample_position,body_center),sample_force));
    }
};
} // namespace pf_capillary
