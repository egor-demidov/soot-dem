/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include <filesystem>

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
#include "rect_substrate.h"

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;
using coating_model_t = binary_coating_functor<Eigen::Vector3d, double>;
using rect_substrate_model_t = rect_substrate_with_coating<Eigen::Vector3d, double>;

using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t, coating_model_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double, rect_substrate_model_t>;

using granular_system_t = granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main(int argc, const char ** argv) {
    if (argc < 2) {
        std::cerr << "Path to the input file must be provided as an argument" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto parameter_store = load_parameters(argv[1]);

    print_header(parameter_store, "anchored_restructuring");

    if (parameter_store.simulation_type != "anchored_restructuring") {
        std::cerr << "Parameter file simulation type must be `anchored_restructuring`" << std::endl;
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
    const double substrate_size = get_real_parameter(parameter_store, "substrate_size");

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

    // Parameters for the Van der Waals model
    const double A = get_real_parameter(parameter_store, "A");
    const double h0 = get_real_parameter(parameter_store, "h0");

    // Parameters for the coating model
    const double f_coat_mag = get_real_parameter(parameter_store, "f_coat_max");
    const double f_coat_cutoff = get_real_parameter(parameter_store, "f_coat_cutoff");
    const double f_coat_drop_rate = get_real_parameter(parameter_store, "f_coat_drop_rate");

    // Necking fraction
    const double frac_necks = get_real_parameter(parameter_store, "frac_necks");

    // Parameter for the substrate friction model
    const double k_n_substrate = get_real_parameter(parameter_store, "k_n_substrate");
    const double gamma_n_substrate = get_real_parameter(parameter_store, "gamma_n_substrate");
    const double k_t_substrate = get_real_parameter(parameter_store, "k_t_substrate");
    const double gamma_t_substrate = get_real_parameter(parameter_store, "gamma_t_substrate");
    const double mu_t_substrate = get_real_parameter(parameter_store, "mu_t_substrate");
    const double phi_t_substrate = get_real_parameter(parameter_store, "phi_t_substrate");
    const double k_r_substrate = get_real_parameter(parameter_store, "k_r_substrate");
    const double gamma_r_substrate = get_real_parameter(parameter_store, "gamma_r_substrate");
    const double mu_r_substrate = get_real_parameter(parameter_store, "mu_r_substrate");
    const double phi_r_substrate = get_real_parameter(parameter_store, "phi_r_substrate");
    const double k_o_substrate = get_real_parameter(parameter_store, "k_o_substrate");
    const double gamma_o_substrate = get_real_parameter(parameter_store, "gamma_o_substrate");
    const double mu_o_substrate = get_real_parameter(parameter_store, "mu_o_substrate");
    const double phi_o_substrate = get_real_parameter(parameter_store, "phi_o_substrate");

    // Parameters for the substrate Van der Waals model
    const double A_substrate = get_real_parameter(parameter_store, "A_substrate");
    const double h0_substrate = get_real_parameter(parameter_store, "h0_substrate");

    // Substrate vertices
    const std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> substrate_vertices {
            {-substrate_size / 2.0, -substrate_size / 2.0, 0.0},
            {substrate_size / 2.0, -substrate_size / 2.0 , 0.0},
            {substrate_size / 2.0 , substrate_size / 2.0, 0.0},
            {-substrate_size / 2.0, substrate_size / 2.0, 0.0}
    };

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0 = load_aggregate(parameter_store);

    if (x0.size() == 0) {
        std::cerr << "Loaded an empty aggregate" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Loaded an aggregate of size " << x0.size() << std::endl;

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

    coating_model_t coating_model(f_coat_cutoff, f_coat_mag, f_coat_drop_rate, mass, Eigen::Vector3d::Zero());

    // Create an instance of rectangular substrate model
    rect_substrate_model_t substrate_model {substrate_vertices, x0.size(), k_n_substrate, gamma_n_substrate, k_t_substrate,
                                            gamma_t_substrate, mu_t_substrate, phi_t_substrate, k_r_substrate, gamma_r_substrate, mu_r_substrate,
                                            phi_r_substrate, k_o_substrate, gamma_o_substrate,
                                            mu_o_substrate, phi_o_substrate, A_substrate, h0_substrate, r_part, mass, inertia, dt, f_coat_cutoff - r_part, f_coat_mag, f_coat_drop_rate, Eigen::Vector3d::Zero(), 0.0};

    binary_force_container_t binary_force_functors {aggregate_model, coating_model};

    unary_force_container_t unary_force_functors {substrate_model};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    granular_system_t system(x0.size(), r_verlet, x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    // Count the number of necks
    size_t n_necks = std::count(aggregate_model.get_bonded_contacts().begin(),
                                aggregate_model.get_bonded_contacts().end(), true) / 2;

    auto target_n_necks = size_t(double(n_necks) * frac_necks);

    std::cout << "Breaking " << n_necks - target_n_necks << " necks ..." << std::endl;

    seed_random_engine(rng_seed);

    for (size_t i = n_necks; i > target_n_necks; i --) {
        break_random_neck(aggregate_model.get_bonded_contacts(), x0.size());
    }

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
            substrate_model.dump_mesh("run", n / dump_period);
        }

        system.do_step(dt);
    }

    return 0;
}
