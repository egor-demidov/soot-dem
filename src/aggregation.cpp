//
// Created by egor on 7/8/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include <filesystem>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/contact_force/contact_force.h>
#include <libgran/granular_system/granular_system_neighbor_list.h>

#include "coating_force.h"
#include "aggregate.h"

#include "writer.h"
#include "energy.h"
#include "writer.h"
#include "reader.h"
#include "break_neck.h"
#include "remove_overlap.h"
#include "aggregate_stats.h"
#include "io_common.h"
#include "parameter_loader.h"
#include "random_engine.h"

using contact_force_model_t = contact_force_functor<Eigen::Vector3d, double>;
using hamaker_force_model_t = hamaker_functor<Eigen::Vector3d, double>;
using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, contact_force_model_t, hamaker_force_model_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>;

using granular_system_t = granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        rotational_step_handler, binary_force_container_t, unary_force_container_t>;

class granular_system_neighbor_list_mutable_velocity : public granular_system_t {
public:
    using granular_system_t::granular_system_neighbor_list;

    std::vector<Eigen::Vector3d> & get_v() {
        return this->v;
    }
};

bool particle_overlaps(Eigen::Vector3d const & particle,
                       std::vector<Eigen::Vector3d> const & xs,
                       double r_part) {
    return std::any_of(xs.begin(), xs.end(), [&particle, r_part] (Eigen::Vector3d const & x) {
        return (x - particle).norm() < 2.0 * r_part;
    });
}

Eigen::Vector3d get_random_unit_vector() {
    Eigen::Vector3d vec;
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    do {
        vec = {
                dist(get_random_engine()),
                dist(get_random_engine()),
                dist(get_random_engine())
        };
    } while (vec.norm() == 0);
    return vec.normalized();
}

const std::array<Eigen::Vector3d, 6> box_faces {
        -Eigen::Vector3d::UnitX(),
        Eigen::Vector3d::UnitX(),
        -Eigen::Vector3d::UnitY(),
        Eigen::Vector3d::UnitY(),
        -Eigen::Vector3d::UnitZ(),
        Eigen::Vector3d::UnitZ(),
};

void bounce_off_walls(std::vector<Eigen::Vector3d> const & particles,
                      std::vector<Eigen::Vector3d> & velocities,
                      double r_part, double box_size) {
    for (size_t n = 0; n < box_faces.size(); n ++) {
        for (size_t i = 0; i < particles.size(); i ++) {
            if (velocities[i].dot(box_faces[n]) > 0.0 && box_size / 2.0 - particles[i].dot(box_faces[n]) < r_part) {
                velocities[i] -= 2.0 * velocities[i].dot(box_faces[n]) * box_faces[n];
//                std::cout << "Condition met\n";
            }
        }
    }
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        std::cerr << "Path to the input file must be provided as an argument" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto parameter_store = load_parameters(argv[1]);

    print_header(parameter_store, "aggregation");

    if (parameter_store.simulation_type != "aggregation") {
        std::cerr << "Parameter file simulation type must be `aggregation`" << std::endl;
        exit(EXIT_FAILURE);
    }

    // General simulation parameters
    const double dt = get_real_parameter(parameter_store, "dt");
    const double t_tot = get_real_parameter(parameter_store, "t_tot");
    const auto n_steps = long(t_tot / dt);
    const long n_dumps = get_integer_parameter(parameter_store, "n_dumps");
    const long dump_period = n_steps / n_dumps;
    const long neighbor_update_period = get_integer_parameter(parameter_store, "neighbor_update_period");
    const long n_overlap_iter = get_integer_parameter(parameter_store, "n_overlap_iter");
    const long rng_seed = get_integer_parameter(parameter_store, "rng_seed");
    const long n_part = get_integer_parameter(parameter_store, "n_part");
    const double box_size = get_real_parameter(parameter_store, "box_size");
    const double v0_part = get_real_parameter(parameter_store, "v0_part");

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

    // Parameters for the Van der Waals model
    const double A = get_real_parameter(parameter_store, "A");
    const double h0 = get_real_parameter(parameter_store, "h0");

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    seed_random_engine(rng_seed);
    x0.resize(n_part);
    std::uniform_real_distribution<double> x0_dist(-(box_size / 2.0 - r_part), box_size / 2.0 - r_part);
    for (long n = 0; n < n_part; n ++) {
        Eigen::Vector3d particle;
        do {
            particle = {
                    x0_dist(get_random_engine()),
                    x0_dist(get_random_engine()),
                    x0_dist(get_random_engine())
            };
        } while (particle_overlaps(particle, x0, r_part));
        x0[n] = particle;
    }

    // Generate initial velocities
    v0.resize(x0.size());
    std::transform(v0.begin(), v0.end(), v0.begin(), [v0_part](auto const & v [[maybe_unused]]) {
        return get_random_unit_vector() * v0_part;
    });

    // Fill the remaining buffers with zeros
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    contact_force_model_t contact_model {
            x0.size(),
            k_n, gamma_n,
            k_t, gamma_t, mu_t, phi_t,
            k_r, gamma_r, mu_r, phi_r,
            k_o, gamma_o, mu_o, phi_o,
            r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0};

    hamaker_force_model_t hamaker_model {
        A, h0, r_part, mass, Eigen::Vector3d::Zero(), 0.0
    };

    binary_force_container_t binary_force_functors {contact_model, hamaker_model};

    unary_force_container_t unary_force_functors;

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    granular_system_neighbor_list_mutable_velocity system(x0.size(), r_verlet, x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

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
        }

        system.do_step(dt);

        bounce_off_walls(system.get_x(), system.get_v(), r_part, box_size);
    }

    return 0;
}
