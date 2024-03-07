//
// Created by egor on 2/24/24.
//

#include <iostream>
#include <fstream>
#include <vector>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system_neighbor_list.h>

#include "parameter_loader.h"
#include "rect_substrate.h"
#include "aggregate.h"
#include "writer.h"
#include "energy.h"
#include "remove_overlap.h"
#include "io_common.h"

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;
using rect_substrate_model_t = rect_substrate<Eigen::Vector3d, double>;

using binary_force_container_t =
        binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t>;

using unary_force_container_t =
        unary_force_functor_container<Eigen::Vector3d, double, rect_substrate_model_t>;

using granular_system_t = granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main(int argc, const char ** argv) {
    if (argc < 2) {
        std::cerr << "Path to the input file must be provided as an argument" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto parameter_store = load_parameters(argv[1]);

    print_header(parameter_store, "aggregate_deposition");

    if (parameter_store.simulation_type != "deposition") {
        std::cerr << "Parameter file simulation type must be `deposition`" << std::endl;
        exit(EXIT_FAILURE);
    }

    // General simulation parameters
    const double dt = get_real_parameter(parameter_store, "dt");
    const double t_tot = get_real_parameter(parameter_store, "t_tot");
    const auto n_steps = long(t_tot / dt);
    const long n_dumps = get_integer_parameter(parameter_store, "n_dumps");
    const long dump_period = n_steps / n_dumps;
    const long neighbor_update_period = get_integer_parameter(parameter_store, "neighbor_update_period");
    const double substrate_size = get_real_parameter(parameter_store, "substrate_size");
    const double vz0 = get_real_parameter(parameter_store, "vz0");
    const double rot_x = get_real_parameter(parameter_store, "rot_x");
    const double rot_y = get_real_parameter(parameter_store, "rot_y");
    const double rot_z = get_real_parameter(parameter_store, "rot_z");
    const long n_overlap_iter = get_integer_parameter(parameter_store, "n_overlap_iter");

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

    remove_overlap(x0, r_part, d_crit, n_overlap_iter);

    Eigen::Vector3d center_of_mass = Eigen::Vector3d::Zero();
    for (auto const & x : x0) {
        center_of_mass += x;
    }
    center_of_mass /= double(x0.size());
    Eigen::Matrix3d rot = (Eigen::AngleAxis(rot_x / 180.0 * M_PI, Eigen::Vector3d::UnitX())
                          * Eigen::AngleAxis(rot_y / 180.0 * M_PI, Eigen::Vector3d::UnitY())
                          * Eigen::AngleAxis(rot_z / 180.0 * M_PI, Eigen::Vector3d::UnitZ())).toRotationMatrix();
    std::transform(x0.begin(), x0.end(), x0.begin(), [&center_of_mass, &rot] (auto const & x) -> Eigen::Vector3d {
        return rot * (x - center_of_mass) + center_of_mass;
    });

    double z_min = (*std::min_element(x0.begin(), x0.end(), [](auto const & x1, auto const & x2) -> bool {
        return x1[2] < x2[2];
    }))[2];

    std::transform(x0.begin(), x0.end(), x0.begin(), [r_part, z_min] (auto const & x) {
        return x + Eigen::Vector3d::UnitZ() * (-z_min + 1.1 * r_part);
    });

    // Fill the remaining buffers with zeros
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d {0, 0, -vz0});
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

    // Create an instance of rectangular substrate model
    rect_substrate_model_t substrate_model {substrate_vertices, x0.size(), k_n_substrate, gamma_n_substrate, k_t_substrate,
                                            gamma_t_substrate, mu_t_substrate, phi_t_substrate, k_r_substrate, gamma_r_substrate, mu_r_substrate,
                                            phi_r_substrate, k_o_substrate, gamma_o_substrate,
                                            mu_o_substrate, phi_o_substrate, A_substrate, h0_substrate, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0};

    binary_force_container_t binary_force_functors {aggregate_model};

    unary_force_container_t unary_force_functors {substrate_model};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    granular_system_t system(x0.size(), r_verlet, x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    state_printer_t state_printer(system.get_x(), system.get_v(), system.get_theta(), system.get_omega(), mass, inertia, n_dumps);

    for (long n = 0; n < n_steps; n ++) {
        if (n & neighbor_update_period) {
            system.update_neighbor_list();
        }

        if (n % dump_period == 0) {
            std::cout << state_printer << std::endl;

            dump_particles("run", n / dump_period, system.get_x(), r_part);
            dump_necks("run", n / dump_period, system.get_x(), aggregate_model.get_bonded_contacts(), r_part);
            substrate_model.dump_mesh("run", n / dump_period);
        }

        system.do_step(dt);
    }

    return 0;
}
