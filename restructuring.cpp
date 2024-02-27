//
// Created by egor on 2/26/24.
//


#include <iostream>
#include <fstream>
#include <vector>


#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>

#include "aggregate.h"
#include "writer.h"
#include "energy.h"
#include "reader.h"
#include "central_force.h"

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;

using binary_force_container_t =
        binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t>;

using central_force_model_t =
        central_force_functor<Eigen::Vector3d, double>;

using unary_force_container_t =
        unary_force_functor_container<Eigen::Vector3d, double, central_force_model_t>;

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main() {
    // General simulation parameters
    const double dt = 1e-14;
    const double t_tot = 3.0e-7;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 300;
    const size_t dump_period = n_steps / n_dumps;
    const size_t n_thermo_dumps = 10000;
    const size_t thermo_dump_period = n_steps / n_thermo_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 1.4e-8;
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);

    // Parameters for the contact model
    const double k = 10000.0;
    const double gamma_n = 5.0e-9;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the bonded contact model
    const double k_bond = 100000.0;
    const double gamma_n_bond = 2.0*sqrt(2.0*mass*k_bond);
    const double gamma_t_bond = 0.2 * gamma_n_bond;
    const double gamma_r_bond = 0.05 * gamma_n_bond;
    const double gamma_o_bond = 0.05 * gamma_n_bond;
    const double d_crit = 1.0e-9; // Critical separation

    // Parameters for the Van der Waals model
    const double A = 1.0e-19;
    const double h0 = 1.0e-9;

    // Parameter for the central force model
    const double k_central = 8.0e-6;

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0 = load_mackowski_aggregate("../mackowski_aggregates/aggregate_3.txt", r_part);

    // Fill the remaining buffers with zeros
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    aggregate_model_t aggregate_model(k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
                                      mu_o, phi, k_bond, gamma_n_bond, k_bond, gamma_t_bond, k_bond, gamma_r_bond, k_bond, gamma_o_bond,
                                      d_crit, A, h0, x0, x0.size(), r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    central_force_model_t central_model(x0, k, mass, Eigen::Vector3d::Zero());

    binary_force_container_t binary_force_functors {aggregate_model};

    unary_force_container_t unary_force_functors {central_model};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    granular_system_t system(x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump " << n / dump_period << "/" << n_dumps << std::endl;
        }
        system.do_step(dt);
        central_model.update_center_of_mass(system.get_x());
    }

    return 0;
}