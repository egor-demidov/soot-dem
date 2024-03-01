//
// Created by egor on 3/1/24.
//

#ifndef SOOT_AFM_MAIN_LOOP_H
#define SOOT_AFM_MAIN_LOOP_H

#include <vector>
#include <chrono>

#include <Eigen/Eigen>

template <typename granular_system_t>
void execute_main_loop(granular_system_t & granular_system,
                       std::vector<Eigen::Vector3d> const & x,
                       std::vector<Eigen::Vector3d> const & v,
                       std::vector<Eigen::Vector3d> const & theta,
                       std::vector<Eigen::Vector3d> const & omega,
                       size_t n_steps, size_t dump_period, double dt) {

    // Need to add geometry dumps
    // Need to add code to compute AFM forces
}

#endif //SOOT_AFM_MAIN_LOOP_H
