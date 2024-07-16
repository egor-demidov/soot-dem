/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_REMOVE_OVERLAP_H
#define SOOT_AFM_REMOVE_OVERLAP_H

#include <vector>

#include <Eigen/Eigen>

void remove_overlap(std::vector<Eigen::Vector3d> & x, double r_part, double d_crit, size_t n_iter);

#endif //SOOT_AFM_REMOVE_OVERLAP_H
