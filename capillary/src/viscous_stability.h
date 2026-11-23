#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// The pinned scalar SPH viscosity has dv_i/dt=sum_j b_ij(v_j-v_i), b_ij>=0.
// Inner links use factor 2; phase and mirrored-wall links use factor 4.
// The numerator is cached only while the library's dW*V neighbor data is fixed.
inline double viscous_rate_numerator(double mu, double dW_volume, double distance,
                                     double smoothing_length, double factor) {
    if (!(mu > 0 && smoothing_length > 0 && distance >= 0 &&
          dW_volume <= 0 && std::isfinite(mu) && std::isfinite(smoothing_length) &&
          std::isfinite(distance) && std::isfinite(dW_volume)))
        throw std::runtime_error("Invalid scalar viscous link");
    return -factor * mu * dW_volume / (distance + .01 * smoothing_length);
}
inline double harmonic_viscosity(double mu_i, double mu_j) {
    if (!(mu_i > 0 && mu_j > 0 && std::isfinite(mu_i) && std::isfinite(mu_j)))
        throw std::runtime_error("Invalid viscosity");
    return 2 * std::min(mu_i, mu_j) / (1 + std::min(mu_i, mu_j)/std::max(mu_i, mu_j));
}
inline double viscous_row_rate(double numerator, double density) {
    if (!(numerator >= 0 && std::isfinite(numerator) && density > 0 && std::isfinite(density)))
        throw std::runtime_error("Invalid viscous row rate");
    return numerator / density;
}
inline double viscous_step_limit(double max_rate) {
    if (!(max_rate >= 0 && std::isfinite(max_rate)))
        throw std::runtime_error("Invalid maximum viscous rate");
    // dt*sum b <= .5: a conservative convex-combination bound for the isolated
    // forward-Euler viscosity update. Not a stability proof for all SPH terms.
    return max_rate > 0 ? .5 / max_rate : std::numeric_limits<double>::infinity();
}
