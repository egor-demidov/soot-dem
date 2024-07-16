/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_IO_COMMON_H
#define SOOT_AFM_IO_COMMON_H

#include <vector>
#include <deque>
#include <chrono>

#include <Eigen/Eigen>

#include "parameter_loader.h"

std::vector<Eigen::Vector3d> load_aggregate(parameter_store_t const & parameter_store);

void print_header(parameter_store_t const & parameter_store, std::string const & name);

struct state_printer_t {
    state_printer_t(std::vector<Eigen::Vector3d> const & x,
                    std::vector<Eigen::Vector3d> const & v,
                    std::vector<Eigen::Vector3d> const & theta,
                    std::vector<Eigen::Vector3d> const & omega,
                    double mass,
                    double inertia,
                    long total_n_iter) :
                    x{x}, v{v}, theta{theta}, omega{omega}, mass{mass},
                    inertia{inertia}, total_n_iter{total_n_iter}, n_iter{0} {}

private:
    static constexpr long n_avg_points = 20l;
    std::vector<Eigen::Vector3d> const & x, & v, & theta, & omega;
    const double mass, inertia;
    const long total_n_iter;
    long n_iter;
    std::deque<long> durations;
    std::chrono::high_resolution_clock::time_point prev_time;
    friend std::ostream & operator << (std::ostream & os, state_printer_t & state_printer);
};

#endif //SOOT_AFM_IO_COMMON_H
