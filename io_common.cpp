//
// Created by egor on 3/6/24.
//

#include "reader.h"
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
