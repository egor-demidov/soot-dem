//
// Created by egor on 3/6/24.
//

#include <iomanip>

#include "reader.h"
#include "energy.h"
#include "aggregate_stats.h"
#include "git.h"

#include "io_common.h"

void print_header(parameter_store_t const & parameter_store, std::string const & name) {
    std::cout << name << " (part of soot-dem project) by Egor Demidov" << std::endl;
    std::cout << "Compiler version: " << COMPILER_STRING << std::endl;
    std::cout << "Compiled on: " << SYSTEM_STRING << std::endl;
    std::cout << "Built from git commit " << git::CommitSHA1() << " (" << git::Branch() << ")" << std::endl;

    const char * omp_num_threads = std::getenv("OMP_NUM_THREADS");
    if (omp_num_threads != nullptr)
        std::cout << "OMP_NUM_THREADS: " << omp_num_threads << std::endl;
    else
        std::cout << "OMP_NUM_THREADS was not set" << std::endl;

    if (git::AnyUncommittedChanges()) {
        std::cout << "WARNING: There were uncommitted changes at compilation" << std::endl;
    }

    std::cout << std::endl;
    std::cout << parameter_store << std::endl;
}

std::vector<Eigen::Vector3d> load_aggregate(parameter_store_t const & parameter_store) {
    const std::string aggregate_type = get_string_parameter(parameter_store, "aggregate_type");
    const std::string aggregate_path = get_string_parameter(parameter_store, "aggregate_path");
    const double r_part = get_real_parameter(parameter_store, "r_part");

    std::vector<Eigen::Vector3d> x0;

    if (aggregate_type == "mackowski") {
        x0 = load_mackowski_aggregate(aggregate_path, r_part);
    } else if (aggregate_type == "flage") {
        x0 = load_flage_aggregate(aggregate_path, r_part);
    } else if (aggregate_type == "vtk") {
        x0 = load_vtk_aggregate(aggregate_path, r_part);
    } else {
        std::cerr << "Unrecognized aggregate type: `" << aggregate_type << "`" << std::endl;
        exit(EXIT_FAILURE);
    }

    return x0;
}

std::ostream & operator << (std::ostream & os, state_printer_t & state_printer) {
    if (state_printer.n_iter == 0) {
        state_printer.prev_time = std::chrono::high_resolution_clock::now();
    } else {
        auto current_time = std::chrono::high_resolution_clock::now();
        long duration = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - state_printer.prev_time).count();
        state_printer.durations.emplace_back(duration);
        if (state_printer.durations.size() > state_printer_t::n_avg_points)
            state_printer.durations.pop_front();
        state_printer.prev_time = current_time;
    }

    double ke = compute_ke(state_printer.v, state_printer.theta, state_printer.mass, state_printer.inertia);
    double lin_momentum = compute_linear_momentum(state_printer.v, state_printer.mass);
    double ang_momentum = compute_angular_momentum(state_printer.x, state_printer.v, state_printer.omega, state_printer.mass, state_printer.inertia);

    os << "Dump: " << state_printer.n_iter << "/" << state_printer.total_n_iter << "\tKE: " << ke << "\tP: " << lin_momentum << "\tL: " << ang_momentum;

    double r_g = r_gyration(state_printer.x);

    os << "\tR_g: " << r_g;

    if (state_printer.n_iter < state_printer_t::n_avg_points) {
        os << "\testimating remaining time ...";
    } else {
        long total = std::accumulate(std::begin(state_printer.durations), std::end(state_printer.durations), 0l);
        long duration_per_step = total / state_printer_t::n_avg_points;

        size_t remaining_steps = state_printer.total_n_iter - state_printer.n_iter;
        auto remaining_time_minutes = remaining_steps * duration_per_step / 1000l / 60l;
        auto remaining_time_seconds = remaining_steps * duration_per_step / 1000l % 60l;

        os << "\tremaining: " << remaining_time_minutes << ":" << std::setw(2) << std::setfill('0') << remaining_time_seconds;
    }

    state_printer.n_iter ++;

    return os;
}
