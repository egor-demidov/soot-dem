#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system_neighbor_list.h>

#include "rect_substrate.h"
#include "aggregate.h"
#include "afm_tip.h"
#include "writer.h"
#include "energy.h"
#include "reader.h"
#include "break_neck.h"
#include "remove_overlap.h"
#include "io_common.h"

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;
using rect_substrate_model_t = rect_substrate<Eigen::Vector3d, double>;
using afm_tip_model_t = afm_tip<Eigen::Vector3d, double>;

using binary_force_container_t =
    binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t>;

using unary_force_container_t =
    unary_force_functor_container<Eigen::Vector3d, double, rect_substrate_model_t, afm_tip_model_t>;

using granular_system_t = granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>;

int main(int argc, const char ** argv) {
    if (argc < 2) {
        std::cerr << "Path to the input file must be provided as an argument" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto parameter_store = load_parameters(argv[1]);

    print_header(parameter_store, "afm_necking_fraction");

    if (parameter_store.simulation_type != "afm_necking_fraction") {
        std::cerr << "Parameter file simulation type must be `afm_necking_fraction`" << std::endl;
        exit(EXIT_FAILURE);
    }

    // General simulation parameters
    const double dt = get_real_parameter(parameter_store, "dt");
    const double t_tot = get_real_parameter(parameter_store, "t_tot");
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = get_integer_parameter(parameter_store, "n_dumps");
    const size_t dump_period = n_steps / n_dumps;
    const size_t n_thermo_dumps = get_integer_parameter(parameter_store, "n_force_pts");
    const size_t thermo_dump_period = n_steps / n_thermo_dumps;
    const size_t neighbor_list_update_period = get_integer_parameter(parameter_store, "neighbor_update_period");
    const double substrate_size = get_real_parameter(parameter_store, "substrate_size");
    const long n_overlap_iter = get_integer_parameter(parameter_store, "n_overlap_iter");
    const double v_afm = get_real_parameter(parameter_store, "v_afm");
    const double tip_x0 = get_real_parameter(parameter_store, "tip_x0");
    const double tip_y0 = get_real_parameter(parameter_store, "tip_y0");
    const double tip_z0 = get_real_parameter(parameter_store, "tip_z0");
    const double t_reversal = get_real_parameter(parameter_store, "t_reversal");

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

    // Parameter for the tip friction model
    const double k_n_tip = get_real_parameter(parameter_store, "k_n_tip");
    const double gamma_n_tip = get_real_parameter(parameter_store, "gamma_n_tip");
    const double k_t_tip = get_real_parameter(parameter_store, "k_t_tip");
    const double gamma_t_tip = get_real_parameter(parameter_store, "gamma_t_tip");
    const double mu_t_tip = get_real_parameter(parameter_store, "mu_t_tip");
    const double phi_t_tip = get_real_parameter(parameter_store, "phi_t_tip");
    const double k_r_tip = get_real_parameter(parameter_store, "k_r_tip");
    const double gamma_r_tip = get_real_parameter(parameter_store, "gamma_r_tip");
    const double mu_r_tip = get_real_parameter(parameter_store, "mu_r_tip");
    const double phi_r_tip = get_real_parameter(parameter_store, "phi_r_tip");
    const double k_o_tip = get_real_parameter(parameter_store, "k_o_tip");
    const double gamma_o_tip = get_real_parameter(parameter_store, "gamma_o_tip");
    const double mu_o_tip = get_real_parameter(parameter_store, "mu_o_tip");
    const double phi_o_tip = get_real_parameter(parameter_store, "phi_o_tip");

    // Parameters for the tip Van der Waals model
    const double A_tip = get_real_parameter(parameter_store, "A_tip");
    const double h0_tip = get_real_parameter(parameter_store, "h0_tip");

    // Parameter for AFM transfer function
    const double omega_0_trs = get_real_parameter(parameter_store, "omega_trans");

    // Necking fraction
    const double frac_necks = get_real_parameter(parameter_store, "frac_necks");

    // Substrate vertices
    const std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> substrate_vertices {
            {-substrate_size / 2.0, -substrate_size / 2.0, 0.0},
            {substrate_size / 2.0, -substrate_size / 2.0 , 0.0},
            {substrate_size / 2.0 , substrate_size / 2.0, 0.0},
            {-substrate_size / 2.0, substrate_size / 2.0, 0.0}
    };

    const Eigen::Vector3d afm_tip_offset {tip_x0, tip_y0, tip_z0};

    // Afm tip base vertices
    const std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> afm_base_vertices {
        Eigen::Vector3d {-6.0 * r_part, -3.0 * r_part, 25.0 * r_part} + afm_tip_offset,
        Eigen::Vector3d {6.0 * r_part, -3.0 * r_part, 25.0 * r_part} + afm_tip_offset,
        Eigen::Vector3d {0.0, 6.0 * r_part, 25.0 * r_part} + afm_tip_offset
    };

    // Afm tip peak vertex
    const Eigen::Vector3d afm_peak_vertex = Eigen::Vector3d {
        0.0, 0, 20.0 * r_part
    } + afm_tip_offset;

    // Afm tip initial velocity
//    const double afm_tip_t_max = /*2.0 * */t_tot * 0.018 /*2.1e-07*/;
    auto afm_tip_v = [v_afm, t_reversal] (double t) -> Eigen::Vector3d {
        return Eigen::Vector3d{0, 0, -v_afm} + Eigen::Vector3d{0, 0, 2.0 * v_afm}
            * (0.5 + 0.5 * tanh(100000000.0 * (t - t_reversal)));
    };

    // Declare a buffer for AFM force-dispalcement data
    std::vector<double> displacement_buffer, force_buffer;
    displacement_buffer.reserve(n_thermo_dumps);
    force_buffer.reserve(n_thermo_dumps);

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0 = load_aggregate(parameter_store);

    remove_overlap(x0, r_part, d_crit, n_overlap_iter);

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

    // Create an instance of rectangular substrate model
    rect_substrate_model_t substrate_model {substrate_vertices, x0.size(), k_n_substrate, gamma_n_substrate, k_t_substrate,
                                            gamma_t_substrate, mu_t_substrate, phi_t_substrate, k_r_substrate, gamma_r_substrate, mu_r_substrate,
                                            phi_r_substrate, k_o_substrate, gamma_o_substrate,
                                            mu_o_substrate, phi_o_substrate, A_substrate, h0_substrate, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0};

    // Create an instance of afm tip model
    afm_tip_model_t afm_tip_model {afm_base_vertices, afm_peak_vertex, afm_tip_v(0.0), x0.size(), k_n_tip, gamma_n_tip, k_t_tip,
                                   gamma_t_tip, mu_t_tip, phi_t_tip, k_r_tip, gamma_r_tip, mu_r_tip,
                                   phi_r_tip, k_o_tip, gamma_o_tip,
                                   mu_o_tip, phi_o_tip, A_tip, h0_tip, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0};


    binary_force_container_t binary_force_functors {aggregate_model};

    unary_force_container_t unary_force_functors {substrate_model, afm_tip_model};

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

    for (size_t i = n_necks; i > target_n_necks; i --) {
        break_random_neck(aggregate_model.get_bonded_contacts(), x0.size());
    }

    state_printer_t state_printer(system.get_x(), system.get_v(), system.get_theta(), system.get_omega(), mass, inertia, n_dumps);

    double f_apprent  = 0.0;
    double df_apparent = 0.0;

    afm_tip_model.toggle_force_accumulation();

    std::filesystem::create_directory("run");

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % neighbor_list_update_period) {
            system.update_neighbor_list();
        }
        if (n % dump_period == 0) {

            std::cout << state_printer << std::endl;

            dump_particles("run", n / dump_period, system.get_x(), r_part);
            dump_necks("run", n / dump_period, system.get_x(), aggregate_model.get_bonded_contacts(), r_part);
            substrate_model.dump_mesh("run", n / dump_period);
            afm_tip_model.dump_mesh("run", n / dump_period);
        }

        system.do_step(dt);

        double f_applied = afm_tip_model.force_accumulator;
        f_apprent += df_apparent * dt;
        df_apparent += omega_0_trs * omega_0_trs * (f_applied - f_apprent - 2.0 / omega_0_trs * df_apparent) * dt;
        afm_tip_model.force_accumulator = 0.0;

        if (n % thermo_dump_period == 0) {
            force_buffer.emplace_back(f_apprent);
            displacement_buffer.emplace_back(afm_tip_model.get_tip_position()[2]);
        }

        afm_tip_model.update_positions(dt);
        afm_tip_model.update_velocities(afm_tip_v(double(n) * dt));
    }

    std::ofstream ofs("force-displacement.csv");

    ofs << "dz, F\n";
    for (size_t i = 0; i < displacement_buffer.size(); i ++) {
        ofs << displacement_buffer[i] << ", " << force_buffer[i] << "\n";
    }

    return 0;
}
