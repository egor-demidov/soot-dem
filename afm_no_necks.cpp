#include <iostream>
#include <fstream>
#include <vector>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/contact_force/contact_force.h>
#include <libgran/granular_system/granular_system.h>

#include "rect_substrate.h"
#include "afm_tip.h"
#include "writer.h"
#include "energy.h"
#include "reader.h"

using vdw_model_t = hamaker_functor<Eigen::Vector3d, double>;
using contact_model_t = contact_force_functor<Eigen::Vector3d, double>;
using rect_substrate_model_t = rect_substrate<Eigen::Vector3d, double>;
using afm_tip_model_t = afm_tip<Eigen::Vector3d, double>;

using binary_force_container_t =
    binary_force_functor_container<Eigen::Vector3d, double, vdw_model_t, contact_model_t>;

using unary_force_container_t =
    unary_force_functor_container<Eigen::Vector3d, double, rect_substrate_model_t, afm_tip_model_t>;

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

    // Parameter for AFM transfer fucntion
    const double omega_0_trs = 70.0e-3 * 2.0 * M_PI * 1e9;

    // Substrate vertices
    const std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> substrate_vertices {
        {-25.0 * r_part, -25.0 * r_part, 0.0},
        {25.0 * r_part, -25.0 * r_part, 0.0},
        {25.0 * r_part, 25.0 * r_part, 0.0},
        {-25.0 * r_part, 25.0 * r_part, 0.0}
    };

    const Eigen::Vector3d afm_tip_offset {8.0 * r_part, 7.0 * r_part, -5.0 * r_part};

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
    const double afm_tip_t_max = /*2.0 * */t_tot * 0.4 /*2.1e-07*/;
    auto afm_tip_v = [t_tot, afm_tip_t_max] (double t) -> Eigen::Vector3d {
        return Eigen::Vector3d{0, 0, -1.0} + Eigen::Vector3d{0, 0, 2.0}
            * (0.5 + 0.5 * tanh(100000000.0 * (t - afm_tip_t_max)));
    };

    // Declare a buffer for AFM force-dispalcement data
    std::vector<double> displacement_buffer, force_buffer;
    displacement_buffer.reserve(n_thermo_dumps);
    force_buffer.reserve(n_thermo_dumps);

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0 = load_vtk_aggregate("../afm_indentation_no_necks/particles_init.vtk", r_part);

    // Fill the remaining buffers with zeros
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());
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

    // Create an instance of afm tip model
    afm_tip_model_t afm_tip_model {afm_base_vertices, afm_peak_vertex, afm_tip_v(0.0), x0.size(), k, gamma_n, k,
        gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
        mu_o, phi, A, h0, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0};


    binary_force_container_t binary_force_functors {vdw_model, contact_model};

    unary_force_container_t unary_force_functors {substrate_model, afm_tip_model};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
        step_handler_instance;

    granular_system_t system(x0,
        v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
        step_handler_instance, binary_force_functors, unary_force_functors);

    double f_apprent  = 0.0;
    double df_apparent = 0.0;

    afm_tip_model.toggle_force_accumulation();

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump " << n / dump_period << "/" << n_dumps << std::endl;
            dump_particles("run", n / dump_period, system.get_x(), r_part);
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