//
// Created by egor on 2/24/24.
//

#include "energy.h"

double compute_ke(std::vector<Eigen::Vector3d> const & v, std::vector<Eigen::Vector3d> const & omega, double mass, double inertia) {
    double ke = 0.0;
    for (auto const & pt : v) {
        ke += 1.0 / 2.0 * mass * pt.dot(pt);
    }
    for (auto const & pt : omega) {
        ke += 1.0 / 2.0 * inertia * pt.dot(pt);
    }
    return ke;
}
