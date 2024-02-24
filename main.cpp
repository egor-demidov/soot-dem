#include <iostream>
#include <fstream>
#include <vector>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>

#include "rect_substrate.h"
#include "aggregate.h"
#include "afm_tip.h"

using aggregate_model_t = aggregate<Eigen::Vector3d, double>;
using rect_substrate_model_t = rect_substrate<Eigen::Vector3d, double>;
using afm_tip_model_t = afm_tip<Eigen::Vector3d, double>;

using binary_force_container_t =
    binary_force_functor_container<Eigen::Vector3d, double, aggregate_model_t>;

using unary_force_container_t =
    unary_force_functor_container<Eigen::Vector3d, double, rect_substrate_model_t, afm_tip_model_t>;

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>;

void dump_particle_positions(std::string const & dir, size_t count,
    std::vector<Eigen::Vector3d> const & x, double r_part) {

    std::stringstream out_file_name;
    out_file_name << dir << "/particles_" << count << ".csv";
    std::ofstream ofs(out_file_name.str());

    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        exit(EXIT_FAILURE);
    }

    ofs << "x, y, z\n";
    for (auto const & point : x) {
        ofs << point[0] / r_part << ", " << point[1] / r_part << ", " << point[2] / r_part << "\n";
    }
}

std::vector<Eigen::Vector3d> load_mackowski_aggregate(std::string const & path, double r_part) {
    std::ifstream ifs(path);

    if (!ifs.good()) {
        std::cerr << "Unable to read the aggregate file: " << path << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::vector<Eigen::Vector3d> x0;

    while (getline(ifs, line)) {
        if (!line.empty()) {
            std::istringstream oss(line);
            double _, x, y, z;
            oss >> _ >> x >> y >> z;
            x0.emplace_back(x * r_part, y * r_part, z * r_part);
        }
    }

    return x0;
}

int main() {
    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 8.0e-7;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 500;
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

    // Parameters for the bonded contact model
    const double k_bond = 1000.0;
    const double gamma_n_bond = 2.0*sqrt(2.0*mass*k_bond);
    const double gamma_t_bond = 0.2 * gamma_n_bond;
    const double gamma_r_bond = 0.05 * gamma_n_bond;
    const double gamma_o_bond = 0.05 * gamma_n_bond;
    const double d_crit = 1.0e-9; // Critical separation

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

    // Afm tip base vertices
    const std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> afm_base_vertices {
        {2.0 * r_part, -3.0 * r_part, 30.0 * r_part},
        {8.0 * r_part, -3.0 * r_part, 30.0 * r_part},
        {5.0 * r_part, 3.0 * r_part, 30.0 * r_part}
    };

    // Afm tip peak vertex
    const Eigen::Vector3d afm_peak_vertex {
        5.0 * r_part, 0, 25.0 * r_part
    };

    // Afm tip initial velocity
    const double afm_tip_t_max = 1.8e-07;
    auto afm_tip_v = [t_tot, afm_tip_t_max] (double t) -> Eigen::Vector3d {
        return Eigen::Vector3d{0, 0, -1.0} + Eigen::Vector3d{0, 0, 2.0}
            * (0.5 + 0.5 * tanh(100000000.0 * (t - afm_tip_t_max)));
    };

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    x0 = load_mackowski_aggregate("../mackowski_aggregates/aggregate_2.txt", r_part);

    std::transform(x0.begin(), x0.end(), x0.begin(), [r_part] (auto const & x) {
        return x + Eigen::Vector3d::UnitZ() * 14.0 * r_part;
    });

    // Fill the remaining buffers with zeros
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d {0, 0, -1.0});
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    aggregate_model_t aggregate_model(k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
      mu_o, phi, k_bond, gamma_n_bond, k_bond, gamma_t_bond, k_bond, gamma_r_bond, k_bond, gamma_o_bond,
      d_crit, A, h0, x0, x0.size(), r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of rectangular substrate model
    rect_substrate_model_t substrate_model {substrate_vertices, x0.size(), k, gamma_n, k,
        gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
        mu_o, phi, A, h0, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0};

    // Create an instance of afm tip model
    afm_tip_model_t afm_tip_model {afm_base_vertices, afm_peak_vertex, afm_tip_v(0.0), x0.size(), k, gamma_n, k,
        gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
        mu_o, phi, A, h0, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0};


    binary_force_container_t binary_force_functors {aggregate_model};

    unary_force_container_t unary_force_functors {substrate_model, afm_tip_model};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
        step_handler_instance;

    granular_system_t system(x0,
        v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
        step_handler_instance, binary_force_functors, unary_force_functors);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump " << n / dump_period << "/" << n_dumps << std::endl;
            dump_particle_positions("run", n / dump_period, system.get_x(), r_part);
            substrate_model.dump_mesh("run", n / dump_period);
            afm_tip_model.dump_mesh("run", n / dump_period);
        }
        system.do_step(dt);
        afm_tip_model.update_positions(dt);
        afm_tip_model.update_velocities(afm_tip_v(double(n) * dt));
    }

    return 0;
}
