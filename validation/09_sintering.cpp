/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include <filesystem>
#include <fstream>
#include <iostream>

#include <Eigen/Eigen>

#include <libgran/granular_system/granular_system.h>
#include <libgran/sinter_bridge/alt_sinter_bridge.h>
#include <libgran/contact_force/contact_force.h>

#include "../src/energy.h"
#include "../src/writer.h"

static const std::filesystem::path OUTPUT_FILE = "/dev/stdout";

using contact_functor_t = contact_force_functor<Eigen::Vector3d, double>;
using sinter_functor_t = alt_sinter_functor<Eigen::Vector3d, double>; // Bonded contact force
using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, sinter_functor_t>;

using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>;

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main() {
    // General simulation parameters
    const double dt = 1e-14;
    const double t_tot = 1.0e-7;
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

    // Parameters for the non-bonded contact model
    const double k = 10000.0;
    const double gamma_n = 5.0e-9;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the bonded contact model
    const double k_bond = 1000.0;
    const double gamma_n_bond = 1.25e-8;
    const double gamma_t_bond = 0.2 * gamma_n_bond;
    const double gamma_r_bond = 0.05 * gamma_n_bond;
    const double gamma_o_bond = 0.05 * gamma_n_bond;
    const double d_crit = 1.0e-9;

    // Initialize the particles
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
    x0.emplace_back(0, 0, 0);
    x0.emplace_back(0, 2.0*r_part, 0);
    x0.emplace_back(2.0*r_part, 2.0*r_part, 0);

    v0.emplace_back(1, 0, 0);
    v0.emplace_back(0, 0, 0);
    v0.emplace_back(0, 1, 0);

    // Initialize the remaining buffers
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    contact_functor_t contact_model(x0.size(),
                                    k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                    r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    sinter_functor_t sinter_model(x0.size(), x0,
                                  k_bond, gamma_n_bond, k_bond, gamma_t_bond, k_bond, gamma_r_bond, k_bond, gamma_o_bond,
                                  r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0, d_crit, contact_model);

    binary_force_container_t
            binary_force_functors{sinter_model};

    unary_force_container_t unary_force_functors;

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;

    // Create an instance of granular_system using the contact force model
    // Using velocity Verlet integrator for rotational systems and a default
    // step handler for rotational systems
    granular_system_t system(x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(),
                             0.0, step_handler_instance,
                             binary_force_functors, unary_force_functors);

    // Buffers for thermo data
    std::vector<double> t_span, ke_trs_span, ke_rot_span, ke_tot_span, lm_span, am_span, lm_span_norm, am_span_norm;

    for (size_t n = 0; n < n_steps; n ++) {
//        if (n % dump_period == 0) {
//            std::cout << "Dump #" << n / dump_period << std::endl;
//            dump_particles("run", n / dump_period, system.get_x(), system.get_theta(),
//                           system.get_v(), system.get_a(),
//                           system.get_omega(), system.get_alpha(), r_part);
//            dump_necks("run", n / dump_period, system.get_x(), std::vector<bool>{false, true, false, true, false, true, false, false,
//                                                                                 true}, r_part);
//        }
        if (n % thermo_dump_period == 0) {
            auto ke_trs = compute_ke_trs(system.get_v(), mass);
            auto ke_rot = compute_ke_rot(system.get_omega(), inertia);
            auto ke_tot = ke_trs + ke_rot;
            auto lm = compute_linear_momentum(system.get_v(), mass);
            auto am = compute_angular_momentum(system.get_x(), system.get_v(), system.get_omega(), mass, inertia);
            t_span.emplace_back(double(n) * dt);
            ke_trs_span.emplace_back(ke_trs);
            ke_rot_span.emplace_back(ke_rot);
            ke_tot_span.emplace_back(ke_tot);
            lm_span.emplace_back(lm);
            am_span.emplace_back(am);
        }
        system.do_step(dt);
    }

    std::ofstream ofs(OUTPUT_FILE);

    if (!ofs.good()) {
        std::cerr << "Unable to create a data file" << std::endl;
        return EXIT_FAILURE;
    }

    double ke_trs_max = *std::max_element(ke_trs_span.begin(), ke_trs_span.end());
    double ke_rot_max = *std::max_element(ke_rot_span.begin(), ke_rot_span.end());
    double ke_tot_max = *std::max_element(ke_tot_span.begin(), ke_tot_span.end());
    double lm_max = *std::max_element(lm_span.begin(), lm_span.end());
    double am_max = *std::max_element(am_span.begin(), am_span.end());

    lm_span_norm.resize(lm_span.size());
    am_span_norm.resize(am_span.size());

    // Normalize the buffers
    std::transform(ke_trs_span.begin(), ke_trs_span.end(), ke_trs_span.begin(), [ke_trs_max, ke_tot_max] (auto ke) {
        return ke / ke_tot_max;
    });
    std::transform(ke_rot_span.begin(), ke_rot_span.end(), ke_rot_span.begin(), [ke_rot_max, ke_tot_max] (auto ke) {
        return ke / ke_tot_max;
    });
    std::transform(lm_span.begin(), lm_span.end(), lm_span_norm.begin(), [lm_max] (auto lm) {
        return lm / lm_max;
    });
    std::transform(am_span.begin(), am_span.end(), am_span_norm.begin(), [am_max] (auto am) {
        return am / am_max;
    });

    // Write the data
    ofs << "t\tE_trs\tE_rot\tE_tot\tP\tL\tPnorm\tLnorm\n";
    for (size_t i = 0; i < ke_trs_span.size(); i ++) {
        ofs << t_span[i] << "\t"
            << ke_trs_span[i] << "\t"
            << ke_rot_span[i] << "\t"
            << ke_tot_span[i] << "\t"
            << lm_span[i] << "\t"
            << am_span[i] << "\t"
            << lm_span_norm[i] << "\t"
            << am_span_norm[i] << "\n";
    }

    return 0;
}