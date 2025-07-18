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
#include "rect_substrate.h"

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;
using coating_model_t = hemisphere_coating_model_t<Eigen::Vector3d, double>;
using rect_substrate_model_t = rect_substrate_with_coating<Eigen::Vector3d, double>;

using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double, rect_substrate_model_t, coating_model_t>;

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
    const double e_crit = get_real_parameter(parameter_store, "e_crit");

    // Parameters for the Van der Waals model
    const double A = get_real_parameter(parameter_store, "A");
    const double h0 = get_real_parameter(parameter_store, "h0");

    // Parameters for the coating model
    const double f_coat_max = get_real_parameter(parameter_store, "f_coat_max");
    const double f_coat_cutoff = get_real_parameter(parameter_store, "f_coat_cutoff");
    const double f_coat_drop_rate = get_real_parameter(parameter_store, "f_coat_drop_rate");
    const double drag_coefficient = get_real_parameter(parameter_store, "drag_coefficient");

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

    // Neck Path
    bool has_necks_path = has_path_parameter(parameter_store, "necks_path");

    // Neck strength
    bool has_neck_strength_parameters = has_real_parameter(parameter_store, "mean_neck_radius") && has_real_parameter(parameter_store, "std_neck_radius") && has_real_parameter(parameter_store, "neck_strength_constant");

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

    // Load necks if path to necks has been passes or send errors is no neck parameters passes or both
    if(has_necks_path && has_neck_strength_parameters){
        std::cerr << "Both `necks_path` and `neck strength parameters` were provided, but are mutually exclusive" << std::endl;
        exit(EXIT_FAILURE);
    }else if(!has_necks_path && !has_neck_strength_parameters){
        std::cerr << "Either `necks_path` or `neck strength parameters` must be provided for this simulation type" << std::endl;
        exit(EXIT_FAILURE);
    }else if (has_necks_path) {
        const std::filesystem::path necks_path = get_path_parameter(parameter_store, "necks_path");

        auto bonded_contacts = load_necks(necks_path, x0.size());

        aggregate_model.get_bonded_contacts() = bonded_contacts;

        size_t n_necks = std::count(aggregate_model.get_bonded_contacts().begin(),
                                    aggregate_model.get_bonded_contacts().end(), true) / 2;
        std::cout << "Loaded necks " << n_necks << " from file" << std::endl;
    }

    coating_model_t coating_model(f_coat_max, mass, Eigen::Vector3d::Zero(), r_part, x0, drag_coefficient, t_tot);
    // Create an instance of rectangular substrate model
    rect_substrate_model_t substrate_model {substrate_vertices, x0.size(), k_n_substrate, gamma_n_substrate, k_t_substrate,
                                            gamma_t_substrate, mu_t_substrate, phi_t_substrate, k_r_substrate, gamma_r_substrate, mu_r_substrate,
                                            phi_r_substrate, k_o_substrate, gamma_o_substrate,
                                            mu_o_substrate, phi_o_substrate, A_substrate, h0_substrate, r_part, mass, inertia, dt, f_coat_cutoff - r_part, f_coat_max, f_coat_drop_rate, Eigen::Vector3d::Zero(), 0.0};

    binary_force_container_t binary_force_functors {aggregate_model};

    unary_force_container_t unary_force_functors {substrate_model, coating_model};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    granular_system_t system(x0.size(), r_verlet, x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    seed_random_engine(rng_seed);

    state_printer_t state_printer(system.get_x(), system.get_v(), system.get_theta(), system.get_omega(), mass, inertia, n_dumps);

    std::filesystem::create_directory("run");

    // assigning neck strengths with a random normal distribution
    std::vector<double> neck_strengths;
    
    if(has_necks_path){
        const std::filesystem::path necks_path = get_path_parameter(parameter_store, "necks_path");
        neck_strengths = load_neck_strengths(necks_path);

        for(int i = 0; i < neck_strengths.size(); i++){
            std::cout << neck_strengths[i] << " ";
        }
        std::cout << std::endl;
    }else{
        // Parameters for random neck distribution
        const double mean_neck_width = get_real_parameter(parameter_store, "mean_neck_radius");
        const double std_neck_width = get_real_parameter(parameter_store, "std_neck_radius");
        const double neck_strength_constant = get_real_parameter(parameter_store, "neck_strength_constant");

        std::normal_distribution neck_width_dist{mean_neck_width, std_neck_width};
        std::vector<bool> bonded_contacts = aggregate_model.get_bonded_contacts();

        int num_necks = std::count(bonded_contacts.begin(), bonded_contacts.end(), true)/2;

        neck_strengths.resize(num_necks);
        for (long i = 0; i < num_necks; i++) {
            double neck_width = neck_width_dist(get_random_engine());
            neck_strengths[i] = neck_strength_constant * neck_width * neck_width;
        }
    }

    for (long n = 0; n < n_steps; n ++) {
        if (n % neighbor_update_period == 0) {
            system.update_neighbor_list();
        }
        if (n % dump_period == 0) {
            std::cout << state_printer << std::endl;

            dump_particles("run", n / dump_period, system.get_x(),
                           system.get_v(), system.get_a(),
                           system.get_omega(), system.get_alpha(), r_part);
            dump_necks("run", n / dump_period, system.get_x(), aggregate_model.get_bonded_contacts(), r_part, neck_strengths);
            dump_sphere("run", n / dump_period, coating_model.get_COM(), coating_model.get_radius(), r_part);

            // For debugging: write positions of interface particles
            std::vector<Eigen::Vector3d> interface_particle_xs(coating_model.get_interface_particles().size());
            for (long i = 0; i < interface_particle_xs.size(); i ++) {
                interface_particle_xs[i] = system.get_x()[coating_model.get_interface_particles()[i]];
            }
            std::stringstream interface_particle_dump_name;
            interface_particle_dump_name << "run/interface_particles_" << n / dump_period;
            dump_particles(interface_particle_dump_name.str(), interface_particle_xs, r_part);

        }

        coating_model.set_time(n * dt);
        coating_model.updateCOM(system.get_x());
        coating_model.update_interface_particles(system.get_x());
        system.do_step(dt);
        break_strained_necks(aggregate_model, system.get_x(), k_n_bond, k_t_bond, k_r_bond, k_o_bond, neck_strengths, r_part);
    }

    return 0;
}
