/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_ENERGY_H
#define SOOT_AFM_ENERGY_H

#include <vector>

#include <Eigen/Eigen>

double compute_ke_trs(std::vector<Eigen::Vector3d> const & v, double mass);
double compute_ke_rot(std::vector<Eigen::Vector3d> const & omega, double inertia);
double compute_ke(std::vector<Eigen::Vector3d> const & v, std::vector<Eigen::Vector3d> const & omega, double mass, double inertia);
double compute_linear_momentum(std::vector<Eigen::Vector3d> const & v, double mass);
double compute_angular_momentum(std::vector<Eigen::Vector3d> const & x, std::vector<Eigen::Vector3d> const & v,
                                std::vector<Eigen::Vector3d> const & omega, double mass, double inertia);

#endif //SOOT_AFM_ENERGY_H
