//
// Created by egor on 12/23/2025.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

#include <Eigen/Eigen>

#include <libgran/granular_system/granular_system_neighbor_list.h>

#include "../src/aggregate.h"
#include "../src/writer.h"
#include "../src/io_common.h"
#include "../src/random_engine.h"

#include "mechanical_testing_step_handler.h"
#include "force_writer.h"

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;

using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>;

using granular_system_t = granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        mechanical_testing_step_handler, binary_force_container_t, unary_force_container_t>;

int main() {

    // General simulation parameters
    const double dt = 5e-15;
    const double t_tot = 5e-9;
    const auto n_steps = long(t_tot / dt);
    const long n_dumps = 100;
    const long dump_period = n_steps / n_dumps;
    const long neighbor_update_period = 20;
    const long rng_seed = 0;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 14e-9;
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);
    const double r_verlet = 70.0e-9;

    // Parameters for the contact model
    const double k_n = 10000.0;
    const double gamma_n = 5e-9;
    const double k_t = 10000.0;
    const double gamma_t = 1e-9;
    const double mu_t = 1.0;
    const double phi_t = 1.0;
    const double k_r = 10000.0;
    const double gamma_r = 2.5e-10;
    const double mu_r = 0.1;
    const double phi_r = 1.0;
    const double k_o = 10000.0;
    const double gamma_o = 2.5e-10;
    const double mu_o = 0.1;
    const double phi_o = 1.0;

    // Parameters for the bonded contact model
    const double k_n_bond = 1000000.0;
    const double gamma_n_bond = 1.25e-6;
    const double k_t_bond = 10000000.0;
    const double gamma_t_bond = 2.5e-7;
    const double k_r_bond = 10000000.0;
    const double gamma_r_bond = 6.25e-8;
    const double k_o_bond = 10000000.0;
    const double gamma_o_bond = 6.25e-8;
    const double d_crit = 1e-9; // Critical separation

    // Parameters for the Van der Waals model
    const double A = 1.0e-19;
    const double h0 = 1.0e-9;

    // Deformation rate
    const double v_deform = 1.0;

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0 = {
        {0.0, 0.0, -r_part},
        {0.0, 0.0, r_part},
    };

    v0 = {
        Eigen::Vector3d::Zero(),
        {0.0, 0.0, v_deform}
    };

    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    aggregate_model_t aggregate_model {
        k_n, gamma_n,
        k_t, gamma_t, mu_t, phi_t,
        k_r, gamma_r, mu_r, phi_r,
        k_o, gamma_o, mu_o, phi_o,
        k_n_bond, gamma_n_bond,
        k_t_bond, gamma_t_bond,
        k_r_bond, gamma_r_bond,
        k_o_bond, gamma_o_bond,
        d_crit, A, h0, x0, x0.size(),
        r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0};

    binary_force_container_t binary_force_functors {aggregate_model};

    unary_force_container_t unary_force_functors;

    mechanical_testing_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    granular_system_t system(x0.size(), r_verlet, x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    seed_random_engine(rng_seed);

    state_printer_t state_printer(system.get_x(), system.get_v(), system.get_theta(), system.get_omega(), mass, inertia, n_dumps);

    std::filesystem::create_directory("run");

    ForceWriter force_writer("tensile_forces.csv");

    for (long n = 0; n < n_steps; n ++) {
        if (n % neighbor_update_period == 0) {
            system.update_neighbor_list();
        }
        if (n % dump_period == 0) {
            std::cout << state_printer << std::endl;

            dump_particles("run", n / dump_period, system.get_x(),
                           system.get_v(), system.get_a(),
                           system.get_omega(), system.get_alpha(), r_part);
            dump_necks("run", n / dump_period, system.get_x(), aggregate_model.get_bonded_contacts(), r_part);

            force_writer.add_line(n / dump_period,
                    system.get_x()[1] - system.get_x()[0],
                    (system.get_x()[1] - system.get_x()[0]).norm() - 2.0 * r_part,
                    std::get<0>(aggregate_model.get_sinter_model().get_contact_springs()[0]),
                    std::get<1>(aggregate_model.get_sinter_model().get_contact_springs()[0]),
                    std::get<2>(aggregate_model.get_sinter_model().get_contact_springs()[0]),
                    system.get_a()[0] * mass,
                    system.get_alpha()[0] * mass);
        }

        system.do_step(dt);
    }

    return 0;
}
