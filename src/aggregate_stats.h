//
// Created by egor on 3/4/24.
//

#ifndef SOOT_AFM_AGGREGATE_STATS_H
#define SOOT_AFM_AGGREGATE_STATS_H

#include <vector>

#include <Eigen/Eigen>

Eigen::Vector3d center_of_mass(std::vector<Eigen::Vector3d> const & x);

double r_gyration(std::vector<Eigen::Vector3d> const & x);

#endif //SOOT_AFM_AGGREGATE_STATS_H
