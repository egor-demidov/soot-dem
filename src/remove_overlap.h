//
// Created by egor on 3/3/24.
//

#ifndef SOOT_AFM_REMOVE_OVERLAP_H
#define SOOT_AFM_REMOVE_OVERLAP_H

#include <vector>

#include <Eigen/Eigen>

void remove_overlap(std::vector<Eigen::Vector3d> & x, double r_part, double d_crit, size_t n_iter);

#endif //SOOT_AFM_REMOVE_OVERLAP_H
