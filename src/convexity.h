//
// Created by egor on 6/12/25.
//

#ifndef CONVEXITY_H
#define CONVEXITY_H

#include <vector>
#include <Eigen/Eigen>

double compute_convexity();

double compute_z_convexity(std::vector<Eigen::Vector3d> const & pos, double r_part);

#endif //CONVEXITY_H
