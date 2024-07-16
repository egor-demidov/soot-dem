/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include "energy.h"

double compute_ke_trs(std::vector<Eigen::Vector3d> const & v, double mass) {
    double ke = 0.0;
    for (auto const & pt : v) {
        ke += 1.0 / 2.0 * mass * pt.dot(pt);
    }
    return ke;
}

double compute_ke_rot(std::vector<Eigen::Vector3d> const & omega, double inertia) {
    double ke = 0.0;
    for (auto const & pt : omega) {
        ke += 1.0 / 2.0 * inertia * pt.dot(pt);
    }
    return ke;
}

double compute_ke(std::vector<Eigen::Vector3d> const & v, std::vector<Eigen::Vector3d> const & omega, double mass, double inertia) {
    double ke = 0.0;
    ke += compute_ke_rot(omega, inertia);
    ke += compute_ke_trs(v, mass);
    return ke;
}

double compute_linear_momentum(std::vector<Eigen::Vector3d> const & v, double mass) {
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();

    for (auto const & pt : v) {
        acc += pt * mass;
    }

    return acc.norm();
}

double compute_angular_momentum(std::vector<Eigen::Vector3d> const & x, std::vector<Eigen::Vector3d> const & v, std::vector<Eigen::Vector3d> const & omega,
                                double mass, double inertia) {

    Eigen::Vector3d acc = Eigen::Vector3d::Zero();

    for (long i = 0; i < v.size(); i ++) {
        acc += omega[i] * inertia + x[i].cross(v[i] * mass);
    }

    return acc.norm();
}
