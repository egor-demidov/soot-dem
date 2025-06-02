//
// Created by egor on 7/18/24.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system_neighbor_list.h>

#include "coating_force.h"
#include "aggregate.h"

#include "writer.h"
#include "energy.h"
#include "writer.h"
#include "reader.h"
#include "break_neck.h"
#include "aggregate_stats.h"
#include "io_common.h"
#include "parameter_loader.h"
#include "random_engine.h"

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;
using coating_model_t = sphere_coating_functor<Eigen::Vector3d, double>;

using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double, coating_model_t>;

using granular_system_t = granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main(int argc, const char ** argv) {

    if (argc < 2) {
        std::cerr << "Path to the input file must be provided as an argument" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto parameter_store = load_parameters(argv[1]);

    print_header(parameter_store, "restructuring_breaking_necks_sphere_force");

    if (parameter_store.simulation_type != "restructuring_breaking_necks_sphere_force") {
        std::cerr << "Parameter file simulation type must be `restructuring_breaking_necks_sphere_force`" << std::endl;
        exit(EXIT_FAILURE);
    }

    // General simulation parameters
    const double dt = get_real_parameter(parameter_store, "dt");
    const double t_tot = get_real_parameter(parameter_store, "t_tot");
    const auto n_steps = long(t_tot / dt);
    const long n_dumps = get_integer_parameter(parameter_store, "n_dumps");
    const long dump_period = n_steps / n_dumps;
    const long neighbor_update_period = get_integer_parameter(parameter_store, "neighbor_update_period");
    const long rng_seed = get_integer_parameter(parameter_store, "rng_seed");

    // General parameters
    const double rho = get_real_parameter(parameter_store, "rho");
    const double r_part = get_real_parameter(parameter_store, "r_part");
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);
    const double r_verlet = get_real_parameter(parameter_store, "r_verlet");

    // Parameters for the contact model
    const double k_n = get_real_parameter(parameter_store, "k_n");
    const double gamma_n = get_real_parameter(parameter_store, "gamma_n");
    const double k_t = get_real_parameter(parameter_store, "k_t");
    const double gamma_t = get_real_parameter(parameter_store, "gamma_t");
    const double mu_t = get_real_parameter(parameter_store, "mu_t");
    const double phi_t = get_real_parameter(parameter_store, "phi_t");
    const double k_r = get_real_parameter(parameter_store, "k_r");
    const double gamma_r = get_real_parameter(parameter_store, "gamma_r");
    const double mu_r = get_real_parameter(parameter_store, "mu_r");
    const double phi_r = get_real_parameter(parameter_store, "phi_r");
    const double k_o = get_real_parameter(parameter_store, "k_o");
    const double gamma_o = get_real_parameter(parameter_store, "gamma_o");
    const double mu_o = get_real_parameter(parameter_store, "mu_o");
    const double phi_o = get_real_parameter(parameter_store, "phi_o");

    // Parameters for the bonded contact model
    const double k_n_bond = get_real_parameter(parameter_store, "k_n_bond");
    const double gamma_n_bond = get_real_parameter(parameter_store, "gamma_n_bond");
    const double k_t_bond = get_real_parameter(parameter_store, "k_t_bond");
    const double gamma_t_bond = get_real_parameter(parameter_store, "gamma_t_bond");
    const double k_r_bond = get_real_parameter(parameter_store, "k_r_bond");
    const double gamma_r_bond = get_real_parameter(parameter_store, "gamma_r_bond");
    const double k_o_bond = get_real_parameter(parameter_store, "k_o_bond");
    const double gamma_o_bond = get_real_parameter(parameter_store, "gamma_o_bond");
    const double d_crit = get_real_parameter(parameter_store, "d_crit"); // Critical separation
    const double e_crit = get_real_parameter(parameter_store, "e_crit");

    // Parameters for the Van der Waals model
    const double A = get_real_parameter(parameter_store, "A");
    const double h0 = get_real_parameter(parameter_store, "h0");

    // Parameters for the coating model
    const double f_coat_max = get_real_parameter(parameter_store, "f_coat_max");
    // const double f_coat_cutoff = get_real_parameter(parameter_store, "f_coat_cutoff");
    // const double f_coat_drop_rate = get_real_parameter(parameter_store, "f_coat_drop_rate");

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0 = load_aggregate(parameter_store);

    // Fill the remaining buffers with zeros
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());
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

    coating_model_t coating_model(f_coat_max, mass, Eigen::Vector3d::Zero(), r_part, x0, t_tot);

    binary_force_container_t binary_force_functors {aggregate_model};

    unary_force_container_t unary_force_functors {coating_model};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    granular_system_t system(x0.size(), r_verlet, x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    seed_random_engine(rng_seed);

    state_printer_t state_printer(system.get_x(), system.get_v(), system.get_theta(), system.get_omega(), mass, inertia, n_dumps);

    std::filesystem::create_directory("run");

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
            dump_sphere("run", n / dump_period, coating_model.get_COM(), coating_model.get_radius(), r_part);
        }

        coating_model.set_time(n * dt);
        coating_model.updateCOM(system.get_x());
        coating_model.update_interface_particles(system.get_x());
        system.do_step(dt);
        break_strained_necks(aggregate_model, system.get_x(), k_n_bond, k_t_bond, k_r_bond, k_o_bond, e_crit, r_part);
    }

    return 0;
}
