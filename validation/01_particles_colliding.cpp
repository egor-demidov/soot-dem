//
// Created by egor on 6/27/24.
//

#include <filesystem>
#include <fstream>
#include <iostream>

#include <Eigen/Eigen>

#include <libgran/granular_system/granular_system.h>
#include <libgran/contact_force/contact_force.h>

#include "../src/energy.h"

static const std::filesystem::path OUTPUT_FILE = "/dev/stdout";

using contact_force_functor_t = contact_force_functor<Eigen::Vector3d, double>;
using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, contact_force_functor_t>;

using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>;

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main() {
    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 1.0e-10;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 1000;
    const size_t dump_period = n_steps / n_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 1.4e-8;
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);

    // Parameters for the contact model
    const double k = 10000.0;
    const double gamma_n = 0.0;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Initialize the two particles
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
    x0.emplace_back(Eigen::Vector3d::Zero());
    x0.emplace_back(0.0, 2.00001*r_part, 0.0);

    v0.emplace_back(0.0, 1.0, 0.0);
    v0.emplace_back(Eigen::Vector3d::Zero());

    // Initialize the remaining buffers
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    // Create an instance of step_handler
    // Using field type Eigen::Vector3d with container std::vector
    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;

    // Create an instance of contact force model
    // Using field type Eigen::Vector3d with real type double
    contact_force_functor_t contact_force_model(x0.size(),
                                                k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                                r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    binary_force_container_t binary_force_functors{contact_force_model};

    unary_force_container_t unary_force_functors;

    // Create an instance of granular_system using the contact force model
    // Using velocity Verlet integrator for rotational systems and a default
    // step handler for rotational systems
    granular_system_t system(x0, v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(),
     0.0, step_handler_instance, binary_force_functors, unary_force_functors);

    std::vector<double> t_span, ke_trs_span, ke_rot_span, ke_tot_span, lm_span, am_span;


    for (size_t n = 0; n < n_steps; n ++) {
        if (double(n) * dt * 1e6 <= 4.0e-6) {
            t_span.emplace_back(double(n) * dt * 1e6);
            ke_trs_span.emplace_back(compute_ke_trs(system.get_v(), mass));
            ke_rot_span.emplace_back(compute_ke_rot(system.get_omega(), inertia));
            ke_tot_span.emplace_back(ke_trs_span.back() + ke_rot_span.back());
            lm_span.emplace_back(compute_linear_momentum(system.get_v(), mass));
            am_span.emplace_back(compute_angular_momentum(system.get_x(), system.get_v(), system.get_omega(), mass, inertia));
        }
        system.do_step(dt);
    }

    double ke_trs_max = *std::max_element(ke_trs_span.begin(), ke_trs_span.end());
    double ke_rot_max = *std::max_element(ke_rot_span.begin(), ke_rot_span.end());
    double ke_tot_max = *std::max_element(ke_tot_span.begin(), ke_tot_span.end());
    double lm_trs_max = *std::max_element(lm_span.begin(), lm_span.end());
    double am_rot_max = *std::max_element(am_span.begin(), am_span.end());

    std::transform(ke_trs_span.begin(), ke_trs_span.end(), ke_trs_span.begin(), [ke_trs_max, ke_tot_max] (auto ke) {
        return ke / ke_tot_max;
    });
    std::transform(ke_rot_span.begin(), ke_rot_span.end(), ke_rot_span.begin(), [ke_rot_max, ke_tot_max] (auto ke) {
        return ke / ke_tot_max;
    });
    std::transform(lm_span.begin(), lm_span.end(), lm_span.begin(), [lm_trs_max] (auto ke) {
        return ke / lm_trs_max;
    });
    std::transform(am_span.begin(), am_span.end(), am_span.begin(), [am_rot_max] (auto ke) {
        return ke / am_rot_max;
    });

    std::ofstream ofs(OUTPUT_FILE);

    if (!ofs.good()) {
        std::cerr << "Unable to create a data file" << std::endl;
        return EXIT_FAILURE;
    }

    ofs << "t\tE_trs\tE_rot\tp\tl\n";
    for (size_t i = 0; i < t_span.size(); i ++) {
        ofs << t_span[i] << "\t" << ke_trs_span[i] << "\t" << ke_rot_span[i] << "\t" << lm_span[i] << "\t" << am_span[i] << "\n";
    }

    return 0;
}