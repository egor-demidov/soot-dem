#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

// One clock for nested advection/acoustic stepping. Separate elapsed counters
// can disagree at the final endpoint and request a zero acoustic step.
inline double advection_endpoint(double time, double duration, double end) {
    if (!(duration > 0 && std::isfinite(duration)))
        throw std::runtime_error("Invalid advection time step");
    const double endpoint = std::min(end, time + duration);
    if (!(endpoint > time))
        throw std::runtime_error("Advection clock cannot advance");
    return endpoint;
}
inline double checked_acoustic_step(double time, double endpoint, double limit) {
    const double dt = std::min(limit, endpoint - time);
    if (!(limit > 0 && std::isfinite(limit) && dt > 0 && time + dt > time))
        throw std::runtime_error("Invalid or non-advancing acoustic time step");
    return dt;
}
