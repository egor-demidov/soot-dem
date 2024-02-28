//
// Created by egor on 2/24/24.
//

#include <iostream>
#include <fstream>
#include <vector>


#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/contact_force/contact_force.h>
#include <libgran/granular_system/granular_system.h>

#include "rect_substrate.h"
#include "writer.h"
#include "energy.h"
#include "reader.h"

using vdw_model_t = hamaker_functor<Eigen::Vector3d, double>;
using contact_model_t = contact_force_functor<Eigen::Vector3d, double>;
using rect_substrate_model_t = rect_substrate<Eigen::Vector3d, double>;

using binary_force_container_t =
        binary_force_functor_container<Eigen::Vector3d, double, vdw_model_t, contact_model_t>;

using unary_force_container_t =
        unary_force_functor_container<Eigen::Vector3d, double, rect_substrate_model_t>;

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main() {
    // General simulation parameters
    const double dt = 1e-13;
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

    // Parameters for the Van der Waals model
    const double A = 1.0e-19;
    const double h0 = 1.0e-9;

    // Substrate vertices
    const std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> substrate_vertices {
            {-25.0 * r_part, -25.0 * r_part, 0.0},
            {25.0 * r_part, -25.0 * r_part, 0.0},
            {25.0 * r_part, 25.0 * r_part, 0.0},
            {-25.0 * r_part, 25.0 * r_part, 0.0}
    };

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0 = load_mackowski_aggregate("../mackowski_aggregates/aggregate_4.txt", r_part);

    Eigen::Vector3d center_of_mass = Eigen::Vector3d::Zero();
    for (auto const & x : x0) {
        center_of_mass += x;
    }
    center_of_mass /= double(x0.size());
    Eigen::Matrix3d rot = Eigen::AngleAxis(-M_PI / 2.0, Eigen::Vector3d::UnitX()).toRotationMatrix();
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
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d {0, 0, -1.0});
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    contact_model_t contact_model(x0.size(),
                                                                       k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                                                       r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of Hamaker model
    vdw_model_t vdw_model(A, h0,
                                                           r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of rectangular substrate model
    rect_substrate_model_t substrate_model {substrate_vertices, x0.size(), k, gamma_n, k,
                                            gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
                                            mu_o, phi, A, h0, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0};

    binary_force_container_t binary_force_functors {vdw_model, contact_model};

    unary_force_container_t unary_force_functors {substrate_model};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    granular_system_t system(x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump " << n / dump_period << "/" << n_dumps
                << " E: " << compute_ke(system.get_v(), system.get_omega(), mass, inertia) << std::endl;

            dump_particles("run", n / dump_period, system.get_x(), r_part);
            substrate_model.dump_mesh("run", n / dump_period);
        }

        system.do_step(dt);
    }

    return 0;
}
