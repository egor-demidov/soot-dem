// Validation adapter for SPHinXsys v1.0-beta.08. No modifications to upstream physics.
#include "sphinxsys.h"
#include "time_window.h"
#include "viscous_stability.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tbb/global_control.h>
using namespace SPH;
constexpr double pi = 3.14159265358979323846;
struct Options {
    std::string kind = "sessile", viscosity_update = "acoustic";
    double dx = .12, theta = 60, end = 2, sigma = .008, sound = 10, gas = .001, viscosity = .05,
           gap = .3, dt_scale = 1., snapshot_interval = 0.;
    int steps = 1000000, dump = 100, record = 10, threads = 2, audit = 0;
};
Options parse(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; i++) {
        std::string k = argv[i];
        if (k == "--help") {
            std::cout << "--case laplace|sessile|bridge|moving --dx H --theta DEG --end-time T "
                         "--steps N --dump-every N --record-every N --threads N --gap S --sigma S "
                         "--sound-speed C --gas-density R --viscosity MU --dt-scale F --snapshot-interval T --audit-every N --viscosity-update acoustic|advection\n";
            std::exit(0);
        }
        if (i + 1 == argc)
            throw std::runtime_error("Missing value for " + k);
        std::string v = argv[++i];
        if (k == "--case")
            o.kind = v;
        else if (k == "--dx")
            o.dx = std::stod(v);
        else if (k == "--theta")
            o.theta = std::stod(v);
        else if (k == "--end-time")
            o.end = std::stod(v);
        else if (k == "--steps")
            o.steps = std::stoi(v);
        else if (k == "--dump-every")
            o.dump = std::stoi(v);
        else if (k == "--record-every")
            o.record = std::stoi(v);
        else if (k == "--threads")
            o.threads = std::stoi(v);
        else if (k == "--gap")
            o.gap = std::stod(v);
        else if (k == "--sigma")
            o.sigma = std::stod(v);
        else if (k == "--sound-speed")
            o.sound = std::stod(v);
        else if (k == "--gas-density")
            o.gas = std::stod(v);
        else if (k == "--viscosity")
            o.viscosity = std::stod(v);
        else if (k == "--viscosity-update")
            o.viscosity_update = v;
        else if (k == "--audit-every")
            o.audit = std::stoi(v);
        else if (k == "--snapshot-interval")
            o.snapshot_interval = std::stod(v);
        else if (k == "--dt-scale")
            o.dt_scale = std::stod(v);
        else
            throw std::runtime_error("Unknown argument " + k);
    }
    if (o.viscosity_update != "acoustic" && o.viscosity_update != "advection")
        throw std::runtime_error("Unknown viscosity update schedule");
    if (o.kind != "laplace" && o.kind != "sessile" && o.kind != "bridge" && o.kind != "moving")
        throw std::runtime_error("Unknown case");
    if (!(o.dx > 0 && o.dx <= .3 && o.theta > 0 && o.theta < 180 && o.end > 0 && o.sigma > 0 &&
          o.sound > 0 && o.gas > 0 && o.viscosity > 0 && o.gap > 0 && o.dt_scale > 0 &&
          o.dt_scale <= 1 && o.snapshot_interval >= 0 && std::isfinite(o.snapshot_interval)) ||
        o.steps < 1 || o.dump < 0 || o.record < 1 || o.threads < 1 || o.audit < 0)
        throw std::runtime_error("Invalid parameter");
    return o;
}
// Analytic implicit regions. Closest-point operation is used only by library geometry queries.
class Region : public Shape {
    std::function<double(const Vecd &)> sdf_;
    BoundingBox bounds_;

  public:
    Region(std::string name, std::function<double(const Vecd &)> sdf, BoundingBox b)
        : Shape(name), sdf_(sdf), bounds_(b) {}
    bool checkContain(const Vecd &p, bool included = true) override {
        return included ? sdf_(p) <= 0 : sdf_(p) < 0;
    }
    Vecd findClosestPoint(const Vecd &p) override {
        Vecd q = p;
        for (int k = 0; k < 20; k++) {
            double f = sdf_(q);
            Vecd g;
            double e = 1e-6;
            for (int a = 0; a < 3; a++) {
                Vecd d = Vecd::Zero();
                d[a] = e;
                g[a] = (sdf_(q + d) - sdf_(q - d)) / (2 * e);
            }
            if (g.norm() < 1e-12)
                return q;
            q -= f * g / g.squaredNorm();
            if (std::abs(f) < 1e-10)
                break;
        }
        return q;
    }

  protected:
    BoundingBox findBounds() override { return bounds_; }
};
double box(const Vecd &p, const Vecd &half) {
    Vecd q = p.cwiseAbs() - half;
    return q.cwiseMax(0).norm() + std::min(0., q.maxCoeff());
}
struct Stats {
    double m = 0, volume = 0, kinetic = 0, speed = 0, rhoerr = 0;
    Vecd p = Vecd::Zero(), l = Vecd::Zero();
};
Stats stats(BaseParticles &b, double rho0) {
    Stats s;
    for (size_t i = 0; i < b.total_real_particles_; i++) {
        double m = b.mass_[i], rho = b.rho_[i];
        if (!std::isfinite(rho) || rho <= 0 || !b.vel_[i].allFinite() || !b.pos_[i].allFinite())
            throw std::runtime_error("Nonfinite SPH state");
        s.m += m;
        s.volume += m / rho;
        s.p += m * b.vel_[i];
        s.l += b.pos_[i].cross(m * b.vel_[i]);
        s.kinetic += .5 * m * b.vel_[i].squaredNorm();
        s.speed = std::max(s.speed, b.vel_[i].norm());
        s.rhoerr = std::max(s.rhoerr, std::abs(rho / rho0 - 1));
    }
    return s;
}
int run(Options o) {
    tbb::global_control threads(tbb::global_control::max_allowed_parallelism, o.threads);
    bool bridge = o.kind == "bridge" || o.kind == "moving", moving = o.kind == "moving",
         sessile = o.kind == "sessile";
    double H = 1.6, R = .6, BW = 4 * o.dx;
    Vecd centers[2] = {Vecd(0, 0, -R - o.gap / 2), Vecd(0, 0, R + o.gap / 2)};
    BoundingBox domain(Vecd(-H - BW, -H - BW, -H - BW), Vecd(H + BW, H + BW, H + BW));
    auto solid = [&](const Vecd &p) {
        double d = std::max(box(p, Vecd(H + BW, H + BW, H + BW)), -box(p, Vecd(H, H, H)));
        if (sessile)
            d = std::min(d, std::max(box(p, Vecd(H, H, H)), p[2]));
        if (bridge)
            for (auto &c : centers)
                d = std::min(d, (p - c).norm() - R);
        return d;
    };
    auto liquid0 = [&](const Vecd &p) {
        if (bridge) {
            Vecd q = p;
            q[0] -= .04;
            return std::max(std::hypot(q[0], q[1]) - .4, std::abs(p[2]) - .95);
        }
        return p.norm() - R;
    };
    auto liquid = [&](const Vecd &p) { return std::max(liquid0(p), -solid(p)); };
    auto gas = [&](const Vecd &p) {
        return std::max({box(p, Vecd(H, H, H)), -solid(p), -liquid0(p)});
    };
    SPHSystem system(domain, o.dx);
    GlobalStaticVariables::physical_time_ = 0.;
    IOEnvironment io(system);
    FluidBody water(system, makeShared<Region>("Liquid", liquid, domain));
    water.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(1., o.sound,
                                                                              o.viscosity);
    water.generateParticles<ParticleGeneratorLattice>();
    FluidBody air(system, makeShared<Region>("Gas", gas, domain));
    air.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(o.gas, o.sound,
                                                                            o.viscosity * .001);
    air.generateParticles<ParticleGeneratorLattice>();
    SolidBody wall(system, makeShared<Region>("Solid", solid, domain));
    wall.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall.generateParticles<ParticleGeneratorLattice>();
    auto &w = dynamic_cast<FluidParticles &>(water.getBaseParticles());
    auto &a = dynamic_cast<FluidParticles &>(air.getBaseParticles());
    auto &s = dynamic_cast<SolidParticles &>(wall.getBaseParticles());
    ComplexRelation wc(water, {&air});
    ContactRelation ww(water, {&wall});
    ComplexRelation ac(air, {&water});
    ContactRelation aw(air, {&wall});
    ContactRelation sw(wall, {&water, &air});
    SimpleDynamics<NormalDirectionFromBodyShape> normals(wall);
    SimpleDynamics<TimeStepInitialization> initw(water), inita(air);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> densityw(
        ww, wc.getInnerRelation());
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> densitya(aw, ac);
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> transporta(aw, ac);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> advw(water, 1.), adva(air, 1.);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> acousticw(water), acoustica(air);
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressurew(
        ww, wc.getInnerRelation());
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> relaxw(ww,
                                                                             wc.getInnerRelation());
    Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall> pressurea(
        aw, ac, 2.);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall> relaxa(aw, ac);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscousw(ww, wc),
        viscousa(aw, ac);
    InteractionWithUpdate<fluid_dynamics::FreeSurfaceIndicationComplex> detect(
        wc.getInnerRelation(), ww);
    InteractionDynamics<fluid_dynamics::ColorFunctionGradientComplex> gradient(
        wc.getInnerRelation(), ww);
    InteractionDynamics<fluid_dynamics::ColorFunctionGradientInterpolationInner> interpolate(
        wc.getInnerRelation());
    InteractionDynamics<fluid_dynamics::SurfaceNormWithWall> wetting(ww, o.theta * pi / 180.);
    InteractionDynamics<fluid_dynamics::SurfaceTensionAccelerationInner> tension(
        wc.getInnerRelation(), o.sigma);
    InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluidRiemann> reaction(sw);
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscousreaction(sw);
    BodyStatesRecordingToVtp output(io, system.real_bodies_);
    w.addVariableToWrite<Real>("Pressure");
    a.addVariableToWrite<Real>("Pressure");
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    normals.exec();
    std::vector<int> owner(s.total_real_particles_, -1);
    std::vector<Vecd> offset(s.total_real_particles_);
    for (size_t i = 0; i < owner.size(); i++)
        if (bridge)
            for (int k = 0; k < 2; k++)
                if ((s.pos_[i] - centers[k]).norm() <= R + 1e-9) {
                    owner[i] = k;
                    offset[i] = s.pos_[i] - centers[k];
                }
    Vecd V[2] = {Vecd::Zero(), Vecd::Zero()}, omega[2] = {Vecd::Zero(), Vecd::Zero()};
    Eigen::Quaterniond Q[2] = {Eigen::Quaterniond::Identity(), Eigen::Quaterniond::Identity()};
    double drymass = 2. * 4 * pi * R * R * R / 3., inertia = .4 * drymass * R * R;
    std::ofstream hist("history.csv"), loads("loads.csv");
    hist << std::setprecision(17);
    loads << std::setprecision(17);
    hist << "step,time,liquid_mass,liquid_volume,kinetic,max_speed,mach,max_liquid_density_error,"
            "px,py,pz,lx,ly,lz,pressure_inside,pressure_outside,laplace_pressure_reference,raw_"
            "capillary_fx,raw_capillary_fy,raw_capillary_fz,support_impulse_x,support_impulse_y,"
            "support_impulse_z,support_angular_impulse_x,support_angular_impulse_y,support_angular_"
            "impulse_z,liquid_max_speed,gas_max_speed,gas_density_error,liquid_kinetic,gas_kinetic,capillary_impulse_x,capillary_impulse_y,capillary_impulse_z,capillary_angular_impulse_x,capillary_angular_impulse_y,capillary_angular_impulse_z,max_viscous_rate,viscous_dt_limit,last_dt\n";
    loads << "step,time,body,fx,fy,fz,tx,ty,tz,cx,cy,cz,vx,vy,vz,wx,wy,wz\n";
    Vecd capsum = Vecd::Zero(), supportP = Vecd::Zero(), supportL = Vecd::Zero();
    Vecd capillaryP = Vecd::Zero(), capillaryL = Vecd::Zero();
    std::vector<Vecd> capillary_force(w.total_real_particles_, Vecd::Zero());
    std::ofstream audit, outliers;
    if (o.audit) {
        audit.open("momentum_stages.csv");
        outliers.open("gas_outliers.csv");
        audit << std::setprecision(17)
              << "step,time,dt,stage,liquid_px,liquid_py,liquid_pz,gas_px,gas_py,gas_pz,"
                 "fluid_lx,fluid_ly,fluid_lz,prior_fx,prior_fy,prior_fz,wall_fx,wall_fy,wall_fz,capillary_fx,capillary_fy,capillary_fz\n";
        outliers << std::setprecision(17)
                 << "step,time,rank,index,x,y,z,vx,vy,vz,speed,mach,rho,p,ax,ay,az,prior_ax,prior_ay,prior_az,solid_sdf,nearest_liquid_distance,cached_gas_neighbors,cached_liquid_neighbors,cached_wall_neighbors\n";
    }
    double max_viscous_rate = 0., viscous_dt_limit = 0., last_dt = 0.;
    std::vector<double> water_viscous_rows(w.total_real_particles_), air_viscous_rows(a.total_real_particles_);
    auto cache_viscous_rows = [&](ComplexRelation &phase, ContactRelation &boundary,
                                  double mu, double other_mu, double h, std::vector<double> &rows) {
        auto &inner = phase.getInnerRelation().inner_configuration_;
        auto &cross = phase.getContactRelation().contact_configuration_[0];
        auto &wall_neighbors = boundary.contact_configuration_[0];
        for (size_t i = 0; i < rows.size(); i++) {
            double sum = 0.;
            auto add = [&](const Neighborhood &nb, double viscosity, double factor) {
                for (size_t j = 0; j < nb.current_size_; j++)
                    sum += viscous_rate_numerator(viscosity, nb.dW_ijV_j_[j], nb.r_ij_[j], h, factor);
            };
            add(inner[i], mu, 2.);
            add(cross[i], harmonic_viscosity(mu, other_mu), 4.);
            add(wall_neighbors[i], mu, 4.);
            rows[i] = sum;
        }
    };
    int step = 0;
    double time = 0;
    auto record = [&]() {
        Stats ws = stats(w, 1.), as = stats(a, o.gas);
        Vecd P = ws.p + as.p, L = ws.l + as.l;
        double K = ws.kinetic + as.kinetic;
        if (moving)
            for (int k = 0; k < 2; k++) {
                P += drymass * V[k];
                L += centers[k].cross(drymass * V[k]) + inertia * omega[k];
                K += .5 * drymass * V[k].squaredNorm() + .5 * inertia * omega[k].squaredNorm();
            }
        double pin = 0, pout = 0;
        int ni = 0, no = 0;
        for (size_t i = 0; i < w.total_real_particles_; i++)
            if (!bridge && w.pos_[i].norm() < R * .55 && (!sessile || w.pos_[i][2] > .2)) {
                pin += w.p_[i];
                ni++;
            }
        for (size_t i = 0; i < a.total_real_particles_; i++)
            if (a.pos_[i].norm() > .95 && a.pos_[i].norm() < 1.25) {
                pout += a.p_[i];
                no++;
            }
        hist << step << ',' << time << ',' << ws.m << ',' << ws.volume << ',' << K << ','
             << std::max(ws.speed, as.speed) << ',' << std::max(ws.speed, as.speed) / o.sound << ','
             << ws.rhoerr;
        for (auto vec : {P, L})
            for (int j = 0; j < 3; j++)
                hist << ',' << vec[j];
        hist << ',' << (ni ? pin / ni : NAN) << ',' << (no ? pout / no : NAN) << ','
             << 2 * o.sigma / R;
        for (auto vec : {capsum, supportP, supportL})
            for (int j = 0; j < 3; j++)
                hist << ',' << vec[j];
        hist << ',' << ws.speed << ',' << as.speed << ',' << as.rhoerr << ','
             << ws.kinetic << ',' << as.kinetic;
        for (auto vec : {capillaryP, capillaryL})
            for (int j = 0; j < 3; j++) hist << ',' << vec[j];
        hist << ',' << max_viscous_rate << ',' << viscous_dt_limit << ',' << last_dt << '\n';
        hist.flush();
        if (o.audit) {
            std::vector<size_t> indices(a.total_real_particles_);
            for (size_t i = 0; i < indices.size(); i++) indices[i] = i;
            const size_t count = std::min(size_t(8), indices.size());
            std::partial_sort(indices.begin(), indices.begin() + count, indices.end(),
                              [&](size_t i, size_t j) { return a.vel_[i].squaredNorm() > a.vel_[j].squaredNorm(); });
            for (size_t rank = 0; rank < count; rank++) {
                const size_t i = indices[rank];
                double distance = std::numeric_limits<double>::infinity();
                for (size_t j = 0; j < w.total_real_particles_; j++)
                    distance = std::min(distance, (a.pos_[i] - w.pos_[j]).norm());
                outliers << step << ',' << time << ',' << rank << ',' << i;
                for (auto vec : {a.pos_[i], a.vel_[i]})
                    for (int j = 0; j < 3; j++) outliers << ',' << vec[j];
                outliers << ',' << a.vel_[i].norm() << ',' << a.vel_[i].norm()/o.sound
                         << ',' << a.rho_[i] << ',' << a.p_[i];
                for (auto vec : {a.acc_[i], a.acc_prior_[i]})
                    for (int j = 0; j < 3; j++) outliers << ',' << vec[j];
                outliers << ',' << solid(a.pos_[i]) << ',' << distance
                         << ',' << ac.getInnerRelation().inner_configuration_[i].current_size_
                         << ',' << ac.getContactRelation().contact_configuration_[0][i].current_size_
                         << ',' << aw.contact_configuration_[0][i].current_size_ << '\n';
            }
            outliers.flush();
        }
    };
    std::ofstream snapshots("liquid_snapshots.csv");
    snapshots << std::setprecision(17) << "step,time,file\n";
    int last_snapshot = -1;
    auto snapshot = [&](bool write_vtp = false) {
        if (last_snapshot == step) return;
        last_snapshot = step;
        snapshots << step << "," << time << ",liquid_" << step << ".csv\n";
        snapshots.flush();
        std::ofstream f("liquid_" + std::to_string(step) + ".csv");
        f << std::setprecision(17) << "x,y,z,mass,rho,p,surface\n";
        for (size_t i = 0; i < w.total_real_particles_; i++)
            f << w.pos_[i][0] << ',' << w.pos_[i][1] << ',' << w.pos_[i][2] << ',' << w.mass_[i]
              << ',' << w.rho_[i] << ',' << w.p_[i] << ',' << w.surface_indicator_[i] << '\n';
        if (write_vtp)
            output.writeToFile(step);
    };
    record();
    if (o.dump || o.snapshot_interval > 0)
        snapshot(o.dump > 0);
    double next_snapshot = o.snapshot_interval;
    while (time < o.end && step < o.steps) {
        initw.exec();
        inita.exec();
        double Dt = std::min({advw.exec(), adva.exec(), o.end - time});
        densityw.exec();
        densitya.exec();
        transporta.exec();
        if (o.viscosity_update == "advection") {
            // Retained only to reproduce the previous adapter for an A/B test.
            viscousw.exec();
            viscousa.exec();
        }
        cache_viscous_rows(wc, ww, o.viscosity, o.viscosity*.001,
                          water.sph_adaptation_->ReferenceSmoothingLength(), water_viscous_rows);
        cache_viscous_rows(ac, aw, o.viscosity*.001, o.viscosity,
                          air.sph_adaptation_->ReferenceSmoothingLength(), air_viscous_rows);
        detect.exec();
        gradient.exec();
        interpolate.exec();
        wetting.exec();
        std::vector<Vecd> before(w.acc_prior_.begin(),
                                 w.acc_prior_.begin() + w.total_real_particles_);
        tension.exec();
        capsum.setZero();
        for (size_t i = 0; i < w.total_real_particles_; i++) {
            capillary_force[i] = w.mass_[i] * (w.acc_prior_[i] - before[i]);
            capsum += capillary_force[i];
        }
        // Explicit capillary restriction in addition to upstream acoustic/advection limits.
        double capdt = .15 * std::sqrt((1. + o.gas) * std::pow(o.dx, 3) / (2 * pi * o.sigma));
        const double window_end = advection_endpoint(time, Dt, o.end);
        while (time < window_end && step < o.steps) {
            if (o.viscosity_update == "acoustic") {
                // Reset before adding forces: viscosity must not accumulate across substeps.
                initw.exec();
                inita.exec();
                for (size_t i = 0; i < w.total_real_particles_; i++)
                    w.acc_prior_[i] += capillary_force[i] / w.mass_[i];
                // Both phases read the same pre-pressure velocities.
                viscousw.exec();
                viscousa.exec();
            }
            max_viscous_rate = 0.;
            for (size_t i = 0; i < water_viscous_rows.size(); i++)
                max_viscous_rate = std::max(max_viscous_rate, viscous_row_rate(water_viscous_rows[i], w.rho_[i]));
            for (size_t i = 0; i < air_viscous_rows.size(); i++)
                max_viscous_rate = std::max(max_viscous_rate, viscous_row_rate(air_viscous_rows[i], a.rho_[i]));
            viscous_dt_limit = viscous_step_limit(max_viscous_rate);
            const double applied_viscous_limit = o.viscosity_update == "acoustic" ?
                o.dt_scale * viscous_dt_limit : std::numeric_limits<double>::infinity();
            const double limit = std::min({o.dt_scale * acousticw.exec(),
                                           o.dt_scale * acoustica.exec(), o.dt_scale * capdt, applied_viscous_limit});
            const double dt = checked_acoustic_step(time, window_end, limit);
            last_dt = dt;
            // Sample the library's own load before the fluid update. No capillary-force patch is
            // added.
            reaction.exec();
            viscousreaction.exec();
            Vecd F[2] = {Vecd::Zero(), Vecd::Zero()}, T[2] = {Vecd::Zero(), Vecd::Zero()};
            for (size_t i = 0; i < owner.size(); i++) {
                Vecd f = reaction.getForceFromFluid()[i] + viscousreaction.getForceFromFluid()[i];
                if (owner[i] >= 0) {
                    int k = owner[i];
                    F[k] += f;
                    T[k] += (s.pos_[i] - centers[k]).cross(f);
                }
                if (!moving || owner[i] < 0) {
                    supportP += dt * f;
                    supportL += dt * s.pos_[i].cross(f);
                }
            }
            // Pure observation: no correction is applied to a particle or the flow.
            const bool sample_audit = o.audit && step % o.audit == 0;
            Vecd wall_force = Vecd::Zero();
            if (sample_audit)
                for (size_t i = 0; i < owner.size(); i++)
                    wall_force += reaction.getForceFromFluid()[i] + viscousreaction.getForceFromFluid()[i];
            auto audit_stage = [&](const char *stage) {
                if (!sample_audit) return;
                const Stats ws = stats(w, 1.), as = stats(a, o.gas);
                const Vecd angular = ws.l + as.l;
                Vecd prior = Vecd::Zero();
                for (size_t i = 0; i < w.total_real_particles_; i++) prior += w.mass_[i]*w.acc_prior_[i];
                for (size_t i = 0; i < a.total_real_particles_; i++) prior += a.mass_[i]*a.acc_prior_[i];
                audit << step << ',' << time << ',' << dt << ',' << stage;
                for (auto vec : {ws.p, as.p, angular, prior, wall_force, capsum})
                    for (int j = 0; j < 3; j++) audit << ',' << vec[j];
                audit << '\n';
            };
            capillaryP += dt * capsum;
            for (size_t i = 0; i < w.total_real_particles_; i++)
                capillaryL += dt * w.pos_[i].cross(capillary_force[i]);
            audit_stage("before_pressure");
            pressurew.exec(dt);
            audit_stage("after_liquid_pressure");
            pressurea.exec(dt);
            audit_stage("after_gas_pressure");
            relaxw.exec(dt);
            audit_stage("after_liquid_density");
            relaxa.exec(dt);
            audit_stage("after_gas_density");
            if (sample_audit) audit.flush();
            if (moving)
                for (int k = 0; k < 2; k++) {
                    V[k] += dt * F[k] / drymass;
                    omega[k] += dt * T[k] / inertia;
                    centers[k] += dt * V[k];
                    double angle = dt * omega[k].norm();
                    if (angle > 0)
                        Q[k] = Eigen::Quaterniond(Eigen::AngleAxisd(angle, omega[k].normalized())) *
                               Q[k];
                    Q[k].normalize();
                }
            if (moving)
                for (size_t i = 0; i < owner.size(); i++)
                    if (owner[i] >= 0) {
                        int k = owner[i];
                        Vecd r = Q[k] * offset[i];
                        s.pos_[i] = centers[k] + r;
                        s.n_[i] = r.normalized();
                        s.vel_[i] = V[k] + omega[k].cross(r);
                        s.acc_[i] = F[k] / drymass + (T[k] / inertia).cross(r) +
                                    omega[k].cross(omega[k].cross(r));
                    }
            if (moving && (centers[1] - centers[0]).norm() <= 2 * R)
                throw std::runtime_error("Moving audit reached sphere contact; DEM contact law is "
                                         "not part of this validation adapter");
            time = std::min(window_end, time + dt);
            step++;
            GlobalStaticVariables::physical_time_ = time;
            if (step % o.record == 0) {
                record();
                if (bridge)
                    for (int k = 0; k < 2; k++) {
                        loads << step << ',' << time << ',' << k;
                        for (auto v : {F[k], T[k], centers[k], V[k], omega[k]})
                            for (int j = 0; j < 3; j++)
                                loads << ',' << v[j];
                        loads << '\n';
                    }
                loads.flush();
            }
            const bool vtp_due = o.dump && step % o.dump == 0;
            const bool csv_due = o.snapshot_interval > 0 && time >= next_snapshot;
            if (vtp_due || csv_due)
                snapshot(vtp_due);
            if (csv_due)
                next_snapshot = (std::floor(time / o.snapshot_interval) + 1) * o.snapshot_interval;
        }
        water.updateCellLinkedList();
        air.updateCellLinkedList();
        if (moving)
            wall.updateCellLinkedList();
        wc.updateConfiguration();
        ww.updateConfiguration();
        ac.updateConfiguration();
        aw.updateConfiguration();
        sw.updateConfiguration();
    }
    record();
    snapshot(o.dump > 0);
    std::ofstream meta("run.json");
    meta << std::setprecision(17) << "{\"case\":\"" << o.kind << "\",\"steps\":" << step
         << ",\"time\":" << time
         << ",\"completed_time\":" << (time >= o.end - 1e-12 ? "true" : "false")
         << ",\"dx\":" << o.dx << ",\"theta\":" << o.theta << ",\"gap\":" << o.gap
         << ",\"sigma\":" << o.sigma << ",\"radius\":0.6,\"gas_density\":" << o.gas
         << ",\"dt_scale\":" << o.dt_scale << ",\"sound_speed\":" << o.sound
         << ",\"liquid_particles\":" << w.total_real_particles_
         << ",\"viscosity_update\":\"" << o.viscosity_update
         << "\",\"capillary_reaction_verified\":false}\n";
    std::cout << "Completed " << step << " acoustic steps, t=" << time
              << ". Physical validation requires converged reference comparisons.\n";
    return 0;
}
int main(int argc, char **argv) {
    try {
        return run(parse(argc, argv));
    } catch (const std::exception &e) {
        std::cerr << "capillary_validation: " << e.what() << '\n';
        return 1;
    }
}
