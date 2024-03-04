//
// Created by egor on 2/28/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>

#include "coating_force.h"
#include "aggregate.h"

#include "writer.h"
#include "energy.h"
#include "writer.h"
#include "reader.h"
#include "break_neck.h"
#include "remove_overlap.h"
#include "aggregate_stats.h"

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;
using coating_model_t = binary_coating_functor<Eigen::Vector3d, double>;

using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t, coating_model_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>;

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main() {
    // General simulation parameters
    const double dt = 5e-14;
    const double t_tot = 5.0e-7;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 500;
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

    // Parameters for the coating model
    const double f_coat_mag = 1e-11;
    const double f_coat_cutoff = 6.0 * r_part;
    const double f_coat_drop_rate = 10.0 / f_coat_cutoff;

    // Necking fraction
    const double frac_necks = 0.80;

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
  
    x0 = load_mackowski_aggregate("../mackowski_aggregates/aggregate_4.txt", r_part);
    remove_overlap(x0, r_part, d_crit, 1000);

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

    coating_model_t coating_model(f_coat_cutoff, f_coat_mag, f_coat_drop_rate, mass, Eigen::Vector3d::Zero());

    binary_force_container_t binary_force_functors {aggregate_model, coating_model};

    unary_force_container_t unary_force_functors;

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    granular_system_t system(x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    // Count the number of necks
    size_t n_necks = std::count(aggregate_model.get_bonded_contacts().begin(),
                                aggregate_model.get_bonded_contacts().end(), true) / 2;

    auto target_n_necks = size_t(double(n_necks) * frac_necks);

    std::cout << "Breaking " << n_necks - target_n_necks << " necks ..." << std::endl;

    for (size_t i = n_necks; i > target_n_necks; i --) {
        break_random_neck(aggregate_model.get_bonded_contacts(), x0.size());
    }

    auto prev_time = std::chrono::high_resolution_clock::now();
    long mean_time_per_step = 0;

    // Initial radius of gyration
    double rg_0 = r_gyration(x0);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            auto curr_time = std::chrono::high_resolution_clock::now();

            auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - prev_time).count();

            mean_time_per_step = (mean_time_per_step * n / dump_period + time_elapsed) * dump_period / (n + dump_period);

            size_t remaining_steps = n_dumps - n / dump_period;
            auto remaining_time_minutes = remaining_steps * mean_time_per_step / 1000 / 60;
            auto remaining_time_seconds = remaining_steps * mean_time_per_step / 1000 % 60;
            prev_time = curr_time;

            std::cout << "Dump " << n / dump_period << "/" << n_dumps
                      << " E: " << compute_ke(system.get_v(), system.get_omega(), mass, inertia) << " remaining: " <<
                      remaining_time_minutes << ":" << std::setfill('0') << std::setw(2) << remaining_time_seconds << std::endl;

            dump_particles("run", n / dump_period, system.get_x(), r_part);
            dump_necks("run", n / dump_period, system.get_x(), aggregate_model.get_bonded_contacts(), r_part);
        }

        system.do_step(dt);
    }

    // Final radius of gyration
    double rg_f = r_gyration(system.get_x());

    std::cout << "Rg_0: " << rg_0 << " Rg_f: " << rg_f << " Rg_f/Rg_0: " << rg_f / rg_0 << std::endl;

    return 0;
}
