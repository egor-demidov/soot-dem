//
// Created by egor on 7/18/24.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numbers>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system_neighbor_list.h>

#include "applied_force.h"
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


// A custom step handler that prevents velocity updates of select primary particles
template <typename field_container_t, typename field_value_t>
struct mechanical_testing_step_handler : public rotational_step_handler<field_container_t, field_value_t> {
public:
    mechanical_testing_step_handler(size_t particle_i, size_t particle_j)
        : particle_i_(particle_i)
        , particle_j_(particle_j)
    {}

    // This method increments the specified value in the x buffer
    void increment_v(long n,
                     field_value_t const & dx [[maybe_unused]],
                     typename field_container_t::iterator x_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator a_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator omega_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator alpha_begin_itr [[maybe_unused]]) const {
        if (n == particle_i_ || n == particle_j_)
            return;

        rotational_step_handler<field_container_t, field_value_t>::increment_v(n, dx, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

private:
    const size_t particle_i_, particle_j_;
};

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;

using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>;

using granular_system_t = granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        mechanical_testing_step_handler, binary_force_container_t, unary_force_container_t>;

int main(int argc, const char ** argv) {

    if (argc < 2) {
        std::cerr << "Path to the input file must be provided as an argument" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto parameter_store = load_parameters(argv[1]);

    print_header(parameter_store, "mechanics_testing_breaking_necks");

    if (parameter_store.simulation_type != "mechanics_testing_breaking_necks") {
        std::cerr << "Parameter file simulation type must be `mechanics_testing_breaking_necks`" << std::endl;
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

    // Parameters for the Van der Waals model
    const double A = get_real_parameter(parameter_store, "A");
    const double h0 = get_real_parameter(parameter_store, "h0");

    // Parameters for the applied force model
    const long applied_particle_i = get_integer_parameter(parameter_store, "applied_particle_i");
    const long applied_particle_j = get_integer_parameter(parameter_store, "applied_particle_j");
    const double deformation_velocity_magnitude = get_real_parameter(parameter_store, "deformation_velocity_magnitude");
    const double deformation_velocity_ramp_time = get_real_parameter(parameter_store, "deformation_velocity_ramp_time");

    // Parameters for neck breakage
    const double neck_radius = get_real_parameter(parameter_store, "neck_radius");
    const double normal_strength = get_real_parameter(parameter_store, "normal_strength");
    const double shear_strength = get_real_parameter(parameter_store, "shear_strength");

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

    // Deformation velocity parameters
    auto velocity_dir = (x0[applied_particle_j] - x0[applied_particle_i]).normalized();
    double velocity_ramp_slope = deformation_velocity_magnitude / deformation_velocity_ramp_time;

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

    unary_force_container_t unary_force_functors {};

    mechanical_testing_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance(applied_particle_i, applied_particle_j);

    granular_system_t system(x0.size(), r_verlet, x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    seed_random_engine(rng_seed);

    state_printer_t state_printer(system.get_x(), system.get_v(), system.get_theta(), system.get_omega(), mass, inertia, n_dumps);

    std::filesystem::create_directory("run");

    // assigning neck strengths with a random normal distribution
    std::vector<double> normal_yield_forces, shear_yield_forces;
    
    std::vector<bool> bonded_contacts = aggregate_model.get_bonded_contacts();

    int num_necks = std::count(bonded_contacts.begin(), bonded_contacts.end(), true)/2;

    normal_yield_forces.resize(num_necks);
    shear_yield_forces.resize(num_necks);
    for (long i = 0; i < num_necks; i++) {
        normal_yield_forces[i] = std::numbers::pi * neck_radius * neck_radius * normal_strength;
        shear_yield_forces[i] = std::numbers::pi * neck_radius * neck_radius * shear_strength;
    }

    // TODO: REMOVE BEFORE PRODUCTION
    // BREAKS ALL NECKS
    // std::fill(aggregate_model.get_bonded_contacts().begin(), aggregate_model.get_bonded_contacts().end(), false);

    std::filesystem::create_directory("run");

    std::ofstream applied_force_log("applied_force_log.csv");

    applied_force_log << "t, applied_i_x, applied_i_y, applied_i_z, applied_j_x, applied_j_y, applied_j_z, net_i_x, net_i_y, net_i_z, net_j_x, net_j_y, net_j_z\n";

    for (long n = 0; n < n_steps; n ++) {
        double t = static_cast<double>(n) * dt;

        if (n % neighbor_update_period == 0) {
            system.update_neighbor_list();
        }
        if (n % dump_period == 0) {
            std::cout << state_printer << std::endl;

            std::vector<Eigen::Vector3d> xs_pulled {
                system.get_x()[applied_particle_i],
                system.get_x()[applied_particle_j]
            };

            dump_particles("run/pulled_particles_" + std::to_string(n / dump_period), xs_pulled, r_part);

            dump_particles("run", n / dump_period, system.get_x(),
                           system.get_v(), system.get_a(),
                           system.get_omega(), system.get_alpha(), r_part);
            dump_necks("run", n / dump_period, system.get_x(), aggregate_model.get_bonded_contacts(), r_part,
                       k_n_bond, k_t_bond, k_r_bond, k_o_bond, aggregate_model.get_sinter_model().get_contact_springs(), x0.size());

            // auto [applied_a_i, applied_alpha_i] = applied_force_model(applied_particle_i, system.get_x(), system.get_v(), system.get_theta(), system.get_omega(), t);
            // auto [applied_a_j, applied_alpha_j] = applied_force_model(applied_particle_j, system.get_x(), system.get_v(), system.get_theta(), system.get_omega(), t);
            auto applied_force_i = Eigen::Vector3d::Zero() * mass;
            auto applied_force_j = Eigen::Vector3d::Zero() * mass;
            auto net_force_i = system.get_a()[applied_particle_i] * mass;
            auto net_force_j = system.get_a()[applied_particle_j] * mass;

            applied_force_log << t << ", "
                                << applied_force_i[0] << ", "
                                << applied_force_i[1] << ", "
                                << applied_force_i[2] << ", "
                                << applied_force_j[0] << ", "
                                << applied_force_j[1] << ", "
                                << applied_force_j[2] << ", "
                                << net_force_i[0] << ", "
                                << net_force_i[1] << ", "
                                << net_force_i[2] << ", "
                                << net_force_j[0] << ", "
                                << net_force_j[1] << ", "
                                << net_force_j[2] << std::endl;
        }

        // Update deformation velocity
        double vel = std::min(t * velocity_ramp_slope, deformation_velocity_magnitude);
        system.get_v()[applied_particle_i] = -velocity_dir * vel;
        system.get_v()[applied_particle_j] = velocity_dir * vel;

        system.do_step(dt);

        break_strained_necks(aggregate_model, system.get_x(), k_n_bond, k_t_bond, k_r_bond, k_o_bond, shear_yield_forces, normal_yield_forces, r_part);
    }

    return 0;
}
