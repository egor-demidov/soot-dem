/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/granular_system/granular_system_neighbor_list.h>
#include <libgran/hamaker_force/hamaker_force.h>

#include "aggregate.h"
#include "capillary_prm.h"
#include "capillary_setup.h"
#include "coating_force.h"
#include "lb_units.h"

#include "aggregate.h"
#include "aggregate_stats.h"
#include "break_neck.h"
#include "coating_force.h"
#include "io_common.h"
#include "parameter_loader.h"
#include "random_engine.h"
#include "writer.h"

#include "dem_lb_parameters.h"

struct fluid_load_functor {
    using Vec = Eigen::Vector3d;
    using Field = std::vector<Vec>;

    Field forces;  // N
    Field torques; // N m

    double mass;
    double inertia;

    fluid_load_functor(size_t count, double mass_, double inertia_)
        : forces(count, Vec::Zero()), torques(count, Vec::Zero()), mass(mass_), inertia(inertia_) {}

    std::pair<Vec, Vec> operator()(size_t i, const Field & /* positions */,
                                   const Field & /* velocities */, const Field & /* orientations */,
                                   const Field & /* angular velocities */, double /* time */
    ) const {
        return {forces[i] / mass, torques[i] / inertia};
    }
};

using contact_model_t = contact_force_functor<Eigen::Vector3d, double>;

using hamaker_model_t = hamaker_functor<Eigen::Vector3d, double>;

using binary_force_container_t =
    binary_force_functor_container<Eigen::Vector3d, double, contact_model_t, hamaker_model_t>;
using unary_force_container_t =
    unary_force_functor_container<Eigen::Vector3d, double, fluid_load_functor>;

using granular_system_t =
    granular_system_neighbor_list<Eigen::Vector3d, double, rotational_velocity_verlet_half,
                                  rotational_step_handler, binary_force_container_t,
                                  unary_force_container_t>;

class fluid_loaded_system : public granular_system_t {
  public:
    using granular_system_t::granular_system_neighbor_list;

    void set_fluid_loads(fluid_load_functor &loads, const std::vector<Eigen::Vector3d> &new_forces,
                         const std::vector<Eigen::Vector3d> &new_torques) {
        const size_t count = this->get_x().size();

        if (new_forces.size() != count || new_torques.size() != count) {
            throw std::invalid_argument("Expected one fluid force and torque per particle");
        }

        // Validate everything before modifying the state.
        for (size_t i = 0; i < count; ++i) {
            if (!new_forces[i].allFinite() || !new_torques[i].allFinite()) {
                throw std::invalid_argument("Non-finite fluid force or torque");
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
            throw std::invalid_argument("Fluid-loaded half-step integrator requires "
                                        "a fixed positive timestep");
        granular_system_t::do_step(dt);
        lastStepDt_ = dt;
        started_ = true;
    }

    Eigen::Vector3d synchronized_velocity(std::size_t i, double dt) const {
        Eigen::Vector3d velocity = this->get_v().at(i);

        if (started_)
            velocity += 0.5 * dt * this->get_a().at(i);

        return velocity;
    }

    Eigen::Vector3d synchronized_angular_velocity(std::size_t i, double dt) const {
        Eigen::Vector3d omega = this->get_omega().at(i);

        if (started_)
            omega += 0.5 * dt * this->get_alpha().at(i);

        return omega;
    }

    void apply_fluid_impulse(std::size_t particle, const Eigen::Vector3d &linearImpulse,
                             const Eigen::Vector3d &angularImpulse, double particleMass,
                             double particleInertia) {
        if (!linearImpulse.allFinite() || !angularImpulse.allFinite() ||
            !std::isfinite(particleMass) || !std::isfinite(particleInertia) ||
            particleMass <= 0.0 || particleInertia <= 0.0) {
            throw std::invalid_argument("Invalid fluid impulse or particle properties");
        }

        this->v.at(particle) += linearImpulse / particleMass;

        this->omega.at(particle) += angularImpulse / particleInertia;
    }

  private:
    bool started_ = false;
    double lastStepDt_ = 0.0;
};

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " DEM.xml\n";

        return EXIT_FAILURE;
    }

    auto parameter_store = load_parameters(argv[1]);

    std::cout << "capillary_dem: " << argv[1] << '\n';

    if (parameter_store.simulation_type != "capillary_dem") {
        std::cerr << "Parameter file simulation type must be `capillary_dem`" << std::endl;
        exit(EXIT_FAILURE);
    }

    // General simulation parameters
    const double dt = get_real_parameter(parameter_store, "dt");
    const double t_tot = get_real_parameter(parameter_store, "t_tot");
    if (!(std::isfinite(dt) && dt > 0 && std::isfinite(t_tot) && t_tot > 0))
        throw std::runtime_error("dt and t_tot must be finite and positive");
    const double stepsExact = t_tot / dt;
    if (stepsExact > 1e9 || std::abs(stepsExact - std::round(stepsExact)) > 1e-6)
        throw std::runtime_error("t_tot must be a whole number of DEM steps (at most 1e9)");
    const auto n_steps = static_cast<long>(std::llround(stepsExact));
    const long neighbor_update_period =
        get_integer_parameter(parameter_store, "neighbor_update_period");
    if (neighbor_update_period <= 0)
        throw std::runtime_error("neighbor_update_period must be positive");

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

    const long resolution =
        get_integer_parameter(parameter_store, "lb_cells_per_particle_diameter");
    if (resolution < 4 || resolution > 1024)
        throw std::runtime_error("LB resolution must be 4..1024 cells per diameter");
    auto setup = std::make_shared<CapillarySetup>(
        CapillarySetup::load(argv[1], r_part, static_cast<unsigned>(resolution)));
    x0 = setup->positions;
    v0 = setup->velocities;

    // Fill the remaining buffers with zeros
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    contact_model_t contact_model{x0.size(), k_n,   gamma_n, k_t,     gamma_t,
                                  mu_t,      phi_t, k_r,     gamma_r, mu_r,
                                  phi_r,     k_o,   gamma_o, mu_o,    phi_o,
                                  r_part,    mass,  inertia, dt,      Eigen::Vector3d::Zero(),
                                  0.0};

    hamaker_model_t hamaker_model{A, h0, r_part, mass, Eigen::Vector3d::Zero(), 0.0};

    binary_force_container_t binary_force_functors{contact_model, hamaker_model};

    fluid_load_functor fluid_loads(x0.size(), mass, inertia);

    unary_force_container_t unary_force_functors{fluid_loads};

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;

    fluid_loaded_system system(x0.size(), r_verlet, x0, v0, theta0, omega0, 0.0,
                               Eigen::Vector3d::Zero(), 0.0, step_handler_instance,
                               binary_force_functors, unary_force_functors);

    const long diameterCells =
        get_integer_parameter(parameter_store, "lb_cells_per_particle_diameter");

    if (diameterCells <= 0) {
        std::cerr << "lb_cells_per_particle_diameter must be positive\n";

        return EXIT_FAILURE;
    }

    DemLbParameters lbParameters{
        .particleRadiusSI = r_part,
        .particleMassSI = mass,
        .particleInertiaSI = inertia,
        .demTimeStepSI = dt,

        .liquidDensitySI = get_real_parameter(parameter_store, "liquid_density"),

        .liquidDynamicViscositySI = get_real_parameter(parameter_store, "liquid_dynamic_viscosity"),

        .liquidSurfaceTensionSI = get_real_parameter(parameter_store, "liquid_surface_tension"),

        .liquidContactAngleDeg = get_real_parameter(parameter_store, "liquid_contact_angle_deg"),

        .cellsPerParticleDiameter = static_cast<unsigned>(diameterCells)};

    lbParameters.particleCount = system.get_x().size();
    lbParameters.setup = setup;

    lbParameters.submitFluidLoadSI = [&](const DemLbParameters::Field &forcesSI,
                                         const DemLbParameters::Field &torquesSI) {
        const auto count = system.get_x().size();

        if (forcesSI.size() != count || torquesSI.size() != count) {
            throw std::runtime_error("Expected one fluid force and torque per particle");
        }

        system.set_fluid_loads(fluid_loads, forcesSI, torquesSI);
    };

    // Pass only the LB configuration to waLBerla.
    const auto units = makeLbUnits(lbParameters);
    if (static_cast<unsigned long long>(n_steps) * units.lbStepsPerDemStep > 1000000000ULL)
        throw std::runtime_error("Too many LB steps");
    const auto prm = writeCapillaryPrm(*setup, static_cast<unsigned long long>(n_steps) *
                                                   units.lbStepsPerDemStep);
    std::string prmPath = prm.string();
    char *lbArguments[] = {argv[0], prmPath.data(), nullptr};
    std::cout << "Input: particles=" << x0.size() << " liquid=" << setup->liquidMode
              << " dxSI=" << setup->dx << " boxLowerSI=" << setup->lower.transpose()
              << " boxUpperSI=" << setup->upper.transpose() << '\n';
    std::filesystem::create_directories("dem-out");
    if (setup->demFrequency)
        for (const auto &entry : std::filesystem::directory_iterator("dem-out")) {
            const auto name = entry.path().filename().string();
            if (entry.is_regular_file() && name.starts_with("particles_") &&
                entry.path().extension() == ".vtk")
                std::filesystem::remove(entry.path());
        }
    if (setup->demFrequency)
        dump_particles("dem-out", 0, system.get_x(), system.get_v(), system.get_a(),
                       system.get_omega(), system.get_alpha(), r_part);

    lbParameters.particlePositionSI = [&](std::size_t i) -> Eigen::Vector3d {
        return system.get_x().at(i);
    };

    lbParameters.particleVelocitySI = [&](std::size_t i) -> Eigen::Vector3d {
        return system.synchronized_velocity(i, dt);
    };

    lbParameters.particleAngularVelocitySI = [&](std::size_t i) -> Eigen::Vector3d {
        return system.synchronized_angular_velocity(i, dt);
    };

    std::size_t completedDemSteps = 0;

    lbParameters.advanceDemStep = [&]() {
        if (completedDemSteps % neighbor_update_period == 0)
            system.update_neighbor_list();
        system.do_step(dt);
        ++completedDemSteps;

        if (setup->demFrequency && completedDemSteps % setup->demFrequency == 0) {
            auto velocity = system.get_v();
            auto angularVelocity = system.get_omega();
            for (std::size_t i = 0; i < velocity.size(); ++i) {
                velocity[i] = system.synchronized_velocity(i, dt);
                angularVelocity[i] = system.synchronized_angular_velocity(i, dt);
            }
            dump_particles("dem-out", completedDemSteps / setup->demFrequency, system.get_x(),
                           velocity, system.get_a(), angularVelocity, system.get_alpha(), r_part);
        }
    };

    lbParameters.applyRemapImpulseSI = [&](std::size_t i, const Eigen::Vector3d &linearImpulseSI,
                                           const Eigen::Vector3d &angularImpulseSI) {
        if (i >= system.get_x().size()) {
            throw std::out_of_range("Invalid particle index for remapping impulse");
        }

        system.apply_fluid_impulse(i, linearImpulseSI, angularImpulseSI, mass, inertia);
    };

    const int status = runCapillaryLb(2, lbArguments, lbParameters);

    std::cout << "Final DEM state:"
              << " steps=" << completedDemSteps << " timeSI=" << double(completedDemSteps) * dt
              << " positionSI=" << system.get_x().at(0).transpose()
              << " velocitySI=" << system.synchronized_velocity(0, dt).transpose() << '\n';

    return status;
}
