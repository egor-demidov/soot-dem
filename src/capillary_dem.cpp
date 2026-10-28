/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <stdexcept>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system_neighbor_list.h>
#include <libgran/contact_force/contact_force.h>

#include "coating_force.h"
#include "aggregate.h"

#include "coating_force.h"
#include "aggregate.h"
#include "break_neck.h"
#include "aggregate_stats.h"
#include "random_engine.h"
#include "parameter_loader.h"
#include "io_common.h"
#include "writer.h"

#include "dem_lb_parameters.h"

struct fluid_load_functor {
    using Vec = Eigen::Vector3d;
    using Field = std::vector<Vec>;

    Field forces;   // N
    Field torques;  // N m

    double mass;
    double inertia;

    fluid_load_functor(size_t count, double mass_, double inertia_)
        : forces(count, Vec::Zero()),
          torques(count, Vec::Zero()),
          mass(mass_),
          inertia(inertia_)
    {}

    std::pair<Vec, Vec> operator()(
        size_t i,
        const Field& /* positions */,
        const Field& /* velocities */,
        const Field& /* orientations */,
        const Field& /* angular velocities */,
        double /* time */
    ) const {
        return {
            forces[i] / mass,
            torques[i] / inertia
        };
    }
};

using contact_model_t =
    contact_force_functor<Eigen::Vector3d, double>;

using hamaker_model_t =
    hamaker_functor<Eigen::Vector3d, double>;

using binary_force_container_t =
    binary_force_functor_container<
        Eigen::Vector3d, double,
        contact_model_t, hamaker_model_t>;
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double, fluid_load_functor>;

using granular_system_t = granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>;

class fluid_loaded_system : public granular_system_t {
public:
    using granular_system_t::granular_system_neighbor_list;

    void set_fluid_loads(
        fluid_load_functor& loads,
        const std::vector<Eigen::Vector3d>& new_forces,
        const std::vector<Eigen::Vector3d>& new_torques
    ) {
        const size_t count = this->get_x().size();

        if (new_forces.size() != count ||
            new_torques.size() != count) {
            throw std::invalid_argument(
                "Expected one fluid force and torque per particle"
            );
            }

        // Validate everything before modifying the state.
        for (size_t i = 0; i < count; ++i) {
            if (!new_forces[i].allFinite() ||
                !new_torques[i].allFinite()) {
                throw std::invalid_argument(
                    "Non-finite fluid force or torque"
                );
                }
        }

        if (started_) {
            for (size_t i = 0; i < count; ++i) {
                const Eigen::Vector3d deltaAcceleration =
                    (new_forces[i] - loads.forces[i]) / loads.mass;
                const Eigen::Vector3d deltaAngularAcceleration =
                    (new_torques[i] - loads.torques[i]) / loads.inertia;
                // The integrator stores half-step velocities. Replacing a cached
                // force must not change the physical synchronized velocity.
                this->v[i] -= 0.5 * lastStepDt_ * deltaAcceleration;
                this->omega[i] -= 0.5 * lastStepDt_ * deltaAngularAcceleration;
                this->a[i] += deltaAcceleration;
                this->alpha[i] += deltaAngularAcceleration;
            }
        }

        loads.forces = new_forces;
        loads.torques = new_torques;
    }

    void do_step(double dt) {
        if (!std::isfinite(dt) || dt <= 0.0 || (started_ && dt != lastStepDt_))
            throw std::invalid_argument("Fluid-loaded half-step integrator requires a fixed positive timestep");
        granular_system_t::do_step(dt);
        lastStepDt_ = dt;
        started_ = true;
    }

    Eigen::Vector3d synchronized_velocity(
        std::size_t i,
        double dt
    ) const {
        Eigen::Vector3d velocity = this->get_v().at(i);

        if (started_)
            velocity += 0.5 * dt * this->get_a().at(i);

        return velocity;
    }

    Eigen::Vector3d synchronized_angular_velocity(
        std::size_t i,
        double dt
    ) const {
        Eigen::Vector3d omega = this->get_omega().at(i);

        if (started_)
            omega += 0.5 * dt * this->get_alpha().at(i);

        return omega;
    }

    void apply_fluid_impulse(
    std::size_t particle,
    const Eigen::Vector3d& linearImpulse,
    const Eigen::Vector3d& angularImpulse,
    double particleMass,
    double particleInertia)
    {
        if (!linearImpulse.allFinite() ||
            !angularImpulse.allFinite() ||
            !std::isfinite(particleMass) ||
            !std::isfinite(particleInertia) ||
            particleMass <= 0.0 ||
            particleInertia <= 0.0)
        {
            throw std::invalid_argument(
                "Invalid fluid impulse or particle properties");
        }

        this->v.at(particle) += linearImpulse / particleMass;

        this->omega.at(particle) +=
            angularImpulse / particleInertia;
    }

private:
    bool started_ = false;
    double lastStepDt_ = 0.0;
};

int main(int argc, char ** argv) {
    if (argc != 3)
    {
        std::cerr
            << "Usage: " << argv[0]
            << " DEM.xml LB.prm\n";

        return EXIT_FAILURE;
    }

    auto parameter_store = load_parameters(argv[1]);

    print_header(parameter_store, "capillary_dem");

    if (parameter_store.simulation_type != "capillary_dem") {
        std::cerr << "Parameter file simulation type must be `capillary_dem`" << std::endl;
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

    // Parameters for the Van der Waals model
    const double A = get_real_parameter(parameter_store, "A");
    const double h0 = get_real_parameter(parameter_store, "h0");

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    const long n_particles =
    get_integer_parameter(parameter_store, "n_particles");

    if (n_particles != 1 && n_particles != 2) {
        std::cerr << "n_particles must be 1 or 2\n";
        return EXIT_FAILURE;
    }

    const double initial_gap =
        get_real_parameter(parameter_store, "initial_gap");

    const double approach_speed =
        get_real_parameter(parameter_store, "approach_speed");

    if (initial_gap < 0.0) {
        std::cerr << "Use a nonnegative initial_gap for this test\n";
        return EXIT_FAILURE;
    }

    const Eigen::Vector3d center = Eigen::Vector3d::Zero();

    x0.assign(
        static_cast<size_t>(n_particles),
        center
    );

    v0.assign(x0.size(), Eigen::Vector3d::Zero());
    theta0.assign(x0.size(), Eigen::Vector3d::Zero());
    omega0.assign(x0.size(), Eigen::Vector3d::Zero());

    if (n_particles == 2) {
        const double half_separation =
            r_part + 0.5 * initial_gap;

        x0[0].x() -= half_separation;
        x0[1].x() += half_separation;

        v0[0].x() = approach_speed;
        v0[1].x() = -approach_speed;
    }

    // Fill the remaining buffers with zeros
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    contact_model_t contact_model {
        x0.size(),
        k_n, gamma_n,
        k_t, gamma_t, mu_t, phi_t,
        k_r, gamma_r, mu_r, phi_r,
        k_o, gamma_o, mu_o, phi_o,
        r_part,
        mass,
        inertia,
        dt,
        Eigen::Vector3d::Zero(),
        0.0
    };

    hamaker_model_t hamaker_model {
        A,
        h0,
        r_part,
        mass,
        Eigen::Vector3d::Zero(),
        0.0
    };

    binary_force_container_t binary_force_functors {
        contact_model,
        hamaker_model
    };

    fluid_load_functor fluid_loads(
    x0.size(),
    mass,
    inertia
);

    unary_force_container_t unary_force_functors {
        fluid_loads
    };

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance;

    fluid_loaded_system system(x0.size(), r_verlet, x0,
                             v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                             step_handler_instance, binary_force_functors, unary_force_functors);

    const long diameterCells = get_integer_parameter(
    parameter_store,
    "lb_cells_per_particle_diameter"
);

    if (diameterCells <= 0)
    {
        std::cerr
            << "lb_cells_per_particle_diameter must be positive\n";

        return EXIT_FAILURE;
    }

    DemLbParameters lbParameters {
        .particleRadiusSI = r_part,
        .particleMassSI = mass,
        .particleInertiaSI = inertia,
        .demTimeStepSI = dt,

        .liquidDensitySI = get_real_parameter(
            parameter_store,
            "liquid_density"
        ),

        .liquidDynamicViscositySI = get_real_parameter(
            parameter_store,
            "liquid_dynamic_viscosity"
        ),

        .liquidSurfaceTensionSI = get_real_parameter(
            parameter_store,
            "liquid_surface_tension"
        ),

        .liquidContactAngleDeg = get_real_parameter(
            parameter_store,
            "liquid_contact_angle_deg"
        ),

        .cellsPerParticleDiameter =
            static_cast<unsigned>(diameterCells)
    };

    if (system.get_x().size() != 1)
    {
        throw std::runtime_error(
            "This coupling stage currently requires one particle"
        );
    }

    bool printedFirstFluidLoad = false;

    lbParameters.submitFluidLoadSI =
        [&](const Eigen::Vector3d& forceSI,
            const Eigen::Vector3d& torqueSI)
        {
            system.set_fluid_loads(
                fluid_loads,
                {forceSI},
                {torqueSI}
            );

            if (!printedFirstFluidLoad)
            {
                std::cout
                    << "DEM received fluid force [N]: "
                    << forceSI.transpose()
                    << "\nDEM received fluid torque [N m]: "
                    << torqueSI.transpose()
                    << '\n';

                printedFirstFluidLoad = true;
            }
        };

    // Pass only the LB configuration to waLBerla.
    char* lbArguments[] = {
        argv[0],
        argv[2],
        nullptr
    };

    lbParameters.particlePositionSI = [&]() -> Eigen::Vector3d
    {
        return system.get_x().at(0);
    };

    lbParameters.particleVelocitySI =
        [&]() -> Eigen::Vector3d
        {
            return system.synchronized_velocity(0, dt);
        };

    lbParameters.particleAngularVelocitySI =
        [&]() -> Eigen::Vector3d
        {
            return system.synchronized_angular_velocity(0, dt);
        };

    std::size_t completedDemSteps = 0;

    lbParameters.advanceDemStep = [&]()
    {
        // Cheap for the current one-particle case.
        system.update_neighbor_list();

        system.do_step(dt);

        ++completedDemSteps;

        if (completedDemSteps == 1 ||
            completedDemSteps % 100 == 0)
        {
            std::cout
                << "DEM step=" << completedDemSteps
                << " timeSI=" << double(completedDemSteps) * dt
                << " positionSI="
                << system.get_x().at(0).transpose()
                << " velocitySI="
                << system.synchronized_velocity(0, dt).transpose()
                << '\n';
        }
    };

    lbParameters.applyRemapImpulseSI =
        [&](const Eigen::Vector3d& linearImpulseSI,
            const Eigen::Vector3d& angularImpulseSI)
        {
            system.apply_fluid_impulse(
                0,
                linearImpulseSI,
                angularImpulseSI,
                mass,
                inertia);
        };

    const int status = runCapillaryLb(
        2,
        lbArguments,
        lbParameters);

    std::cout
        << "Final DEM state:"
        << " steps=" << completedDemSteps
        << " timeSI=" << double(completedDemSteps) * dt
        << " positionSI=" << system.get_x().at(0).transpose()
        << " velocitySI="
        << system.synchronized_velocity(0, dt).transpose()
        << '\n';

    return status;
}
