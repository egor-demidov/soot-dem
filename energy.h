//
// Created by egor on 2/24/24.
//

#ifndef SOOT_AFM_ENERGY_H
#define SOOT_AFM_ENERGY_H

#include <vector>

#include <Eigen/Eigen>

double compute_ke(std::vector<Eigen::Vector3d> const & v, std::vector<Eigen::Vector3d> const & omega, double mass, double inertia);
double compute_linear_momentum(std::vector<Eigen::Vector3d> const & v, double mass);
double compute_angular_momentum(std::vector<Eigen::Vector3d> const & x, std::vector<Eigen::Vector3d> const & v,
                                std::vector<Eigen::Vector3d> const & omega, double mass, double inertia);

#endif //SOOT_AFM_ENERGY_H
