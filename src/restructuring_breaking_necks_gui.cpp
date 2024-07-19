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
#include "remove_overlap.h"
#include "aggregate_stats.h"
#include "io_common.h"
#include "parameter_loader.h"
#include "random_engine.h"

#include "preview_renderer.h"

static const std::filesystem::path aggregate_file_path = "aggregate_example.vtk";
static constexpr double A = 1e-19;
static constexpr double d_crit = 1e-09;
static constexpr double dt = 5e-15;
static constexpr double f_coat_cutoff = 5.6e-08;
static constexpr double f_coat_drop_rate = 1.78571e+08;
static constexpr double f_coat_max = 1e-09;
static constexpr double gamma_n = 5e-09;
static constexpr double gamma_n_bond = 1.25e-06;
static constexpr double gamma_o = 2.5e-10;
static constexpr double gamma_o_bond = 6.25e-08;
static constexpr double gamma_r = 2.5e-10;
static constexpr double gamma_r_bond = 6.25e-08;
static constexpr double gamma_t = 1e-09;
static constexpr double gamma_t_bond = 2.5e-07;
static constexpr double h0 = 1e-09;
static constexpr double k_n = 10000;
static constexpr double k_n_bond = 1e+06;
static constexpr double k_o = 10000;
static constexpr double k_o_bond = 1e+07;
static constexpr double k_r = 10000;
static constexpr double k_r_bond = 1e+07;
static constexpr double k_t = 10000;
static constexpr double k_t_bond = 1e+07;
static constexpr double mu_o = 0.1;
static constexpr double mu_r = 0.1;
static constexpr double mu_t = 1;
static constexpr double phi_o = 1;
static constexpr double phi_r = 1;
static constexpr double phi_t = 1;
static constexpr double r_part = 1.4e-08;
static constexpr double r_verlet = 7e-08;
static constexpr double rho = 1700;
static constexpr double t_tot = 5e-08;
static constexpr long n_dumps = 500;
static constexpr long n_overlap_iter = 10000;
static constexpr long neighbor_update_period = 20;
static constexpr long rng_seed = 0;

static constexpr auto n_steps = long(t_tot / dt);
static constexpr long dump_period = n_steps / n_dumps;
static const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
static const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;
using coating_model_t = binary_coating_functor<Eigen::Vector3d, double>;

using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t, coating_model_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>;

using granular_system_t = granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main() {


    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0 = load_vtk_aggregate(aggregate_file_path, r_part);

    Renderer renderer(x0, r_part);

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

    coating_model_t coating_model(f_coat_cutoff, f_coat_max, f_coat_drop_rate, mass, Eigen::Vector3d::Zero());

    binary_force_container_t binary_force_functors {aggregate_model, coating_model};

    unary_force_container_t unary_force_functors;

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
            renderer.update_preview(system.get_x(), n / dump_period);
        }

        system.do_step(dt);
        break_strained_necks(aggregate_model, system.get_x(), k_n_bond, k_t_bond, k_r_bond, k_o_bond, r_part);
    }

    return 0;
}
