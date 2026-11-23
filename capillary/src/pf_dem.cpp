// General moving-sphere DEM driver using the same PF-WCSPH operators as the
// SPHinXsys validation adapter. No gas, cell remapping, or momentum projection.
#include "pf_dem_io.h"
#include "pf_dem_add_liquid.h"
#include "pf_density_diagnostics.h"
#ifdef PF_WITH_SOOT_DEM
#include "pf_soot_dem.h"
#endif
#ifdef PF_DEM_USE_TBB
#include "pf_parallel.h"
#include <tbb/global_control.h>
#endif
#include <chrono>
#include <iostream>
#include <memory>
using namespace pf_dem;
double norm(Vector a){return std::sqrt(dot(a,a));}
double kinetic_energy(const System &system) {
    double energy=0.;
    for(size_t i=0;i<system.nf;++i){
        const auto &s=system.state[i];
        energy+=.5*s.mass*dot(s.v,s.v);
    }
    for(const auto &b:system.bodies)
        energy+=.5*b.mass*dot(b.velocity,b.velocity)+.5*dot(b.angular_momentum,b.angular_momentum)/b.inertia();
    return energy;
}
// Bulk stress from the conservative pair virial, which includes both the
// EOS pressure force and the cohesive force. Use the sampled core volume.
void write_droplet_pressure(std::ostream &f,const System &s,const Rates &rates,double volume) {
    constexpr double pi=3.14159265358979323846;
    const double radius=std::cbrt(3*volume/(4*pi));
    const double core_radius=radius-2*s.parameters.support();
    Vector center{};
    double mass=0.;
    for(size_t i=0;i<s.nf;++i) {
        const auto &a=s.state[i];
        center=add(center,scale(a.x,a.mass));
        mass+=a.mass;
    }
    center=scale(center,1./mass);
    size_t count=0;
    double virial=0.,core_volume=0.,eos=0.;
    if(core_radius>s.parameters.dx) {
        for(size_t i=0;i<s.nf;++i) {
            const auto &a=s.state[i];
            const auto offset=subtract(a.x,center);
            if(dot(offset,offset)>=core_radius*core_radius)continue;
            const double v=a.mass/a.rho;
            core_volume+=v;
            virial+=rates.conservative_virial[i];
            eos+=v*pressure(a.rho,s.parameters);
            ++count;
        }
    }
    const double missing=std::numeric_limits<double>::quiet_NaN();
    f<<','<<radius<<','<<core_radius<<','<<count
     <<','<<(count?virial/(3*core_volume):missing)<<','<<(count?eos/core_volume:missing);
}
void check_wall(const System &s,const InitialGeometry &g) {
    const auto &wall=*g.wall;
    const auto &center=g.droplets.front().center;
    const double dx=s.parameters.dx;
    const double interior=std::ceil(wall.half_width/dx)*dx-s.parameters.support();
    for(size_t i=0;i<s.nf;++i) {
        const auto &x=s.state[i].x;
        const double gap=x[2]-wall.z;
        if(gap<-.5*dx)
            throw std::runtime_error("Sessile liquid penetrated the wall by more than dx/2");
        if(gap<s.parameters.support() &&
           std::max(std::abs(x[0]-center[0]),std::abs(x[1]-center[1]))>interior)
            throw std::runtime_error("Sessile liquid approached the plate edge. Increase wall half_width and rerun");
    }
}
int main(int argc,char **argv) try {
    std::string input,resume;double target=-1,equilibration=-1;long long limit=-1;bool density_check=true,legacy_overlap=false;
    std::optional<double> add_volume,post_hold,saturation_ratio,kelvin_length;
    std::optional<double> ke_target_fraction,ke_hold_time;
    std::optional<std::uint64_t> addition_seed;
    std::optional<std::string> addition_mode;
    for(int i=1;i<argc;++i){std::string k=argv[i];
        if(k=="--help"){std::cout<<"pf_dem --input prepared.dat [--resume checkpoint.dat] [--end-time T] [--steps N] [--disable-density-check] [--legacy-overlap-mask] [--fluid-equilibration-time SECONDS] [--add-liquid-volume M3] [--liquid-addition-mode MODE] [--post-addition-equilibration-time SECONDS]\nOptional early stop: --ke-target-fraction FRACTION [--ke-hold-time SECONDS].\nCurvature-based mode requires --saturation-ratio NU --kelvin-length METERS and accepts --liquid-addition-seed SEED.\nPrepare input from XML with tools/run_pf_dem.py. Parameters remain provisional until calibrated.\n";return 0;}
        if(k=="--disable-density-check"){density_check=false;continue;}
        if(k=="--legacy-overlap-mask"){legacy_overlap=true;continue;}
        if(i+1==argc)throw std::runtime_error("Missing value for "+k);
        std::string v=argv[++i];
        if(k=="--input")input=v;else if(k=="--resume")resume=v;
        else if(k=="--end-time")target=std::stod(v);else if(k=="--steps")limit=std::stoll(v);
        else if(k=="--fluid-equilibration-time"){equilibration=std::stod(v);if(!(equilibration>=0&&std::isfinite(equilibration)))throw std::runtime_error("Invalid fluid equilibration duration");}
        else if(k=="--liquid-addition-mode")addition_mode=v;
        else if(k=="--ke-target-fraction"||k=="--ke-hold-time"){
            size_t parsed=0;
            const double value=std::stod(v,&parsed);
            if(parsed!=v.size()||!std::isfinite(value)||value<0||
               (k=="--ke-target-fraction" && !(value>0 && value<1)))
                throw std::runtime_error("Invalid value for "+k);
            if(k=="--ke-target-fraction")ke_target_fraction=value;
            else ke_hold_time=value;
        }
        else if(k=="--liquid-addition-seed"){
            if(v.empty() || v.find_first_not_of("0123456789")!=std::string::npos)
                throw std::runtime_error("Invalid liquid addition seed");
            addition_seed=std::stoull(v);
        }
        else if(k=="--add-liquid-volume"||k=="--post-addition-equilibration-time"||
                k=="--saturation-ratio"||k=="--kelvin-length"){
            size_t parsed=0;
            const double value=std::stod(v,&parsed);
            if(parsed!=v.size() || !std::isfinite(value) || value<0 ||
               ((k=="--add-liquid-volume"||k=="--saturation-ratio") && value==0))
                throw std::runtime_error("Invalid value for "+k);
            if(k=="--add-liquid-volume")add_volume=value;
            else if(k=="--post-addition-equilibration-time")post_hold=value;
            else if(k=="--saturation-ratio")saturation_ratio=value;
            else kelvin_length=value;
        }
        else throw std::runtime_error("Unknown option "+k);
    }
    if(add_volume && resume.empty())throw std::runtime_error("--add-liquid-volume requires --resume");
    if(post_hold && !add_volume)throw std::runtime_error("--post-addition-equilibration-time requires --add-liquid-volume");
    if(addition_mode && !add_volume)throw std::runtime_error("--liquid-addition-mode requires --add-liquid-volume");
    if(ke_hold_time && !ke_target_fraction)throw std::runtime_error("--ke-hold-time requires --ke-target-fraction");
    const std::string placement=addition_mode.value_or("liquid-vapor");
    const auto mode=parse_addition_mode(placement);
    std::optional<KelvinAddition> kelvin;
    if(mode==AdditionMode::curvature_based){
        if(!saturation_ratio || !kelvin_length)
            throw std::runtime_error("Curvature-based addition requires --saturation-ratio and --kelvin-length");
        kelvin=KelvinAddition{*saturation_ratio,*kelvin_length,addition_seed.value_or(0)};
    } else if(saturation_ratio || kelvin_length || addition_seed)
        throw std::runtime_error("Kelvin inputs and the addition seed require --liquid-addition-mode curvature-based");
    auto a=read_input(input);
    if(legacy_overlap&&resume.empty())throw std::runtime_error("--legacy-overlap-mask requires --resume");
    MechanicalFactory factory;
#ifdef PF_WITH_SOOT_DEM
    factory=SootDemLaw::restore;
#else
    if(a.soot.enabled)throw std::runtime_error("Soot neck/contact inputs require the pf_soot_dem target; rebuild with PF_WITH_SOOT_DEM=ON");
#endif
    auto system=resume.empty()?initialize(a.p,a.geometry,a.contact):System::load(resume,factory,legacy_overlap);
#ifdef PF_WITH_SOOT_DEM
    if(a.soot.enabled && resume.empty())system.mechanics=std::make_shared<SootDemLaw>(a.soot,system.bodies);
    if(a.soot.enabled!=bool(system.mechanics))throw std::runtime_error("Checkpoint mechanical backend differs from XML");
#endif
    if(system.mechanics && (a.contact.stiffness>0||a.contact.damping>0||a.contact.friction>0))
        throw std::runtime_error("Disable built-in contacts when soot-dem supplies contacts");
    if(resume.empty())system.configuration_signature=a.signature;
    else if(system.configuration_signature!=a.signature)throw std::runtime_error("Restart input geometry/material/contact configuration differs from checkpoint");
    const bool calibration=system.bodies.empty();
    const bool free_droplet=calibration && !a.geometry.wall;
    if(calibration && equilibration>0)
        throw std::runtime_error("Fluid equilibration requires DEM particles");
    if(resume.empty()){if(equilibration>0)system.begin_equilibration(equilibration);}
    else if(equilibration>=0&&equilibration!=system.equilibration_end)
        throw std::runtime_error("XML equilibration duration differs from checkpoint; equilibration does not restart on resume");
    if(target!=-1){
        if(!(std::isfinite(target)&&target>0))throw std::runtime_error("Invalid target time");
        a.end_time=ke_target_fraction?std::min(a.end_time,target):target;
    }
    if(limit!=-1)a.max_steps=limit;
    if(!(a.end_time>system.time&&std::isfinite(a.end_time))||a.max_steps<1)
        throw std::runtime_error("Target time must exceed checkpoint time; steps must be positive");
    if(add_volume){
        const auto report=add_surface_liquid(system,*add_volume,a.geometry.max_samples,mode,kelvin);
        if(post_hold.value_or(0)>0)system.hold_bodies(*post_hold);
        std::ofstream audit("liquid_addition.json");
        audit<<std::setprecision(17)
             <<"{\n  \"step\": "<<system.step
             <<",\n  \"time_s\": "<<system.time
             <<",\n  \"requested_delta_volume_m3\": "<<*add_volume
             <<",\n  \"added_volume_m3\": "<<report.added_volume
             <<",\n  \"volume_before_m3\": "<<report.volume_before
             <<",\n  \"added_mass_kg\": "<<report.added_mass
             <<",\n  \"added_particles\": "<<report.particles
             <<",\n  \"layers\": "<<report.layers
             <<",\n  \"eligible_sites\": "<<report.eligible_sites
             <<",\n  \"addition_mode\": "<<std::quoted(placement)
             <<",\n  \"added_particle_mass_scale\": "<<report.mass_scale
             <<",\n  \"release_time_s\": "<<(system.equilibration_active?system.equilibration_end:system.time);
        if(kelvin)
            audit<<",\n  \"saturation_ratio\": "<<kelvin->saturation_ratio
                 <<",\n  \"kelvin_length_m\": "<<kelvin->kelvin_length
                 <<",\n  \"addition_seed\": "<<kelvin->seed
                 <<",\n  \"resolved_solid_vapor_sites\": "<<report.resolved_solid_vapor_sites
                 <<",\n  \"resolved_liquid_vapor_sites\": "<<report.resolved_liquid_vapor_sites
                 <<",\n  \"eligible_solid_vapor_sites\": "<<report.eligible_solid_vapor_sites
                 <<",\n  \"eligible_liquid_vapor_sites\": "<<report.eligible_sites-report.eligible_solid_vapor_sites
                 <<",\n  \"added_solid_vapor_particles\": "<<report.added_solid_vapor_particles
                 <<",\n  \"added_liquid_vapor_particles\": "<<report.particles-report.added_solid_vapor_particles;
        audit<<"\n}\n";
        audit.close();
        if(!audit)throw std::runtime_error("Failed liquid addition report");
        std::cout<<"Added "<<report.particles<<" liquid particles at time="<<system.time
                 <<" with initial delta volume="<<report.added_volume<<" m^3\n";
    }
#ifdef PF_DEM_USE_TBB
    tbb::global_control parallelism(tbb::global_control::max_allowed_parallelism,a.threads);
    system.evaluator=[threads=a.threads](const auto&s,const auto&g,const auto&p){return pf_parallel::evaluate(s,g,p,threads);};
    const int threads=a.threads;
#else
    const int threads=1;
    if(a.threads>1)std::cerr<<"pf_dem: built without TBB; running serially\n";
#endif
    std::cout<<"pf_dem: "<<system.bodies.size()<<" DEM particles, "<<system.nf<<" liquid samples, "<<system.sites.size()<<" boundary samples; threads="<<threads<<". Fluid parameters provisional.\n";
    std::ofstream hist("dem_history.csv"),body("dem_particles.csv");hist<<std::setprecision(17);body<<std::setprecision(17);
#ifdef PF_WITH_SOOT_DEM
    const auto soot_law=std::dynamic_pointer_cast<SootDemLaw>(system.mechanics);
    std::ofstream neck_forces;
    if(soot_law){
        neck_forces.open("dem_neck_forces.csv");
        neck_forces<<std::setprecision(17)
                   <<"step,time,particle_i,particle_j,normal_force_N,sliding_force_N,rolling_force_N,torsional_force_N\n";
        if(!neck_forces)throw std::runtime_error("Failed neck forces CSV output");
    }
#endif
    // Historical max_density_error column means deviation from reference rho0,
    // not discretization error or excess above a calibrated equilibrium density.
    hist<<"step,time,liquid_mass,liquid_volume,linear_error,angular_error,kinetic_energy,builtin_contact_energy,max_density_error,max_speed,neighbor_builds,neighbor_reuses,neck_count,fluid_equilibrating,coupled_time";
    if(free_droplet)hist<<",volume_radius,core_radius,core_particle_count,virial_pressure,core_eos_pressure";
    hist<<'\n';
    body<<"step,time,particle,x,y,z,vx,vy,vz,lx,ly,lz,qw,qx,qy,qz,fluid_fx,fluid_fy,fluid_fz,fluid_tx,fluid_ty,fluid_tz,contact_fx,contact_fy,contact_fz,contact_tx,contact_ty,contact_tz,radius,mass,fixed,dem_fx,dem_fy,dem_fz,dem_tx,dem_ty,dem_tz\n";
    double last_dt=0.;
    std::ofstream density("density_extrema.csv");write_density_header(density);
    auto record=[&](const Stage &r){double mass=0,volume=0,err=0,speed=0;
        for(size_t i=0;i<system.nf;++i){const auto &s=system.state[i];mass+=s.mass;volume+=s.mass/s.rho;err=std::max(err,std::abs(s.rho/system.parameters.rho0-1));speed=std::max(speed,norm(s.v));}
        for(size_t j=0;j<system.bodies.size();++j){const auto &b=system.bodies[j];
            body<<system.step<<','<<system.time<<','<<j;
            for(auto v:{b.center,b.velocity,b.angular_momentum})for(double x:v)body<<','<<x;
            for(double x:b.orientation)body<<','<<x;
            for(auto v:{r.hydrodynamic[j].force,r.hydrodynamic[j].torque,r.contact.loads[j].force,r.contact.loads[j].torque})for(double x:v)body<<','<<x;
            body<<','<<b.radius<<','<<b.mass<<','<<b.fixed;
            for(auto v:{r.mechanical[j].force,r.mechanical[j].torque})for(double x:v)body<<','<<x;
            body<<'\n';}
        hist<<system.step<<','<<system.time<<','<<mass<<','<<volume<<','<<norm(system.momentum_error())<<','<<norm(system.momentum_error(true))<<','<<kinetic_energy(system)<<','<<r.contact.energy<<','<<err<<','<<speed<<','<<system.neighbor_builds()<<','<<system.neighbor_reuses()<<','<<(system.mechanics?system.mechanics->neck_count():0)<<','<<int(system.equilibration_active)<<','<<std::max(0.,system.time-system.equilibration_end);
        if(free_droplet)write_droplet_pressure(hist,system,r.fluid,volume);
        hist<<'\n';
        hist.flush();body.flush();if(!hist||!body)throw std::runtime_error("Failed DEM CSV output");
#ifdef PF_WITH_SOOT_DEM
        if(soot_law){
            soot_law->write_neck_forces(neck_forces,system.bodies,system.step,system.time);
            neck_forces.flush();
            if(!neck_forces)throw std::runtime_error("Failed neck forces CSV output");
        }
#endif
        auto [lo,hi]=density_extreme_indices(system);
        auto low=density_probe(system,r,lo),high=density_probe(system,r,hi);
        write_density_probe(density,system,last_dt,"min",low);write_density_probe(density,system,last_dt,"max",high);
        density.flush();if(!density)throw std::runtime_error("Failed density diagnostics output");
        if((density_check&&err>.3) || speed>.2*system.parameters.sound){
            write_vtp(system,"failure_liquid.vtp");
            if(!calibration)write_sphere_vtp(system,"failure_particles.vtp");
            std::ostringstream message;message<<std::setprecision(17)
                <<"DEM-SPH reference-density/Mach operating screen exceeded (density deviation is not an error estimate): rhoMin/rho0="<<low.rho_ratio
                <<" rhoMax/rho0="<<high.rho_ratio<<" maxSpeed="<<speed
                <<" Mach="<<speed/system.parameters.sound
                <<"; see density_extrema.csv and failure_*.vtp; checkpoint is the latest saved checkpoint";
            throw std::runtime_error(message.str());
        }
        std::cout<<"DEM-SPH step="<<system.step<<" time="<<system.time<<" phase="<<(system.equilibration_active?"fluid-equilibration":"coupled")<<" linearError="<<norm(system.momentum_error())<<" angularError="<<norm(system.momentum_error(true))<<'\n';
    };
    auto dump=[&](){std::ostringstream n;n<<std::setw(10)<<std::setfill('0')<<system.step;
        write_vtp(system,"liquid_"+n.str()+".vtp");
        if(!calibration)write_sphere_vtp(system,"particles_"+n.str()+".vtp");
        if(system.mechanics)write_necks_vtp(system,"necks_"+n.str()+".vtp");};
    std::cout<<"Reference-density stopping check: "<<(density_check?"enabled":"disabled")<<"; Mach and invalid-state checks remain enabled\n";
    if(a.geometry.wall) {
        check_wall(system,a.geometry);
        if(a.dump_every)write_wall_vtp(*a.geometry.wall,a.geometry.droplets.front(),a.p.dx);
    }
    double peak_ke=0.,current_ke=0.;
    std::optional<double> below_since;
    auto observe_ke=[&](){
        if(!ke_target_fraction || system.equilibration_active)return false;
        current_ke=kinetic_energy(system);
        if(!std::isfinite(current_ke))throw std::runtime_error("Nonfinite kinetic energy");
        peak_ke=std::max(peak_ke,current_ke);
        if(peak_ke==0 || current_ke/peak_ke>*ke_target_fraction){
            below_since.reset();
            return false;
        }
        if(!below_since)below_since=system.time;
        return system.time-*below_since>=ke_hold_time.value_or(0.);
    };
    observe_ke();
    auto rates=system.rates();record(rates);if(a.dump_every)dump();
    system.save("dem_checkpoint.dat");const auto start=std::chrono::steady_clock::now();
    long long initial_step=system.step;
    bool ke_reached=false;
    while(system.time<a.end_time && system.step-initial_step<a.max_steps && !ke_reached){
        double dt=std::min(a.end_time-system.time,a.dt_scale*system.timestep(rates));
        bool was_equilibrating=system.equilibration_active;
        if(was_equilibrating)dt=std::min(dt,system.equilibration_end-system.time);
        if(!(dt>0&&std::isfinite(dt))||system.time+dt==system.time)throw std::runtime_error("Invalid/nonadvancing DEM timestep");
        system.advance(dt,&rates);last_dt=dt;
        if(a.geometry.wall) {
            try {check_wall(system,a.geometry);}
            catch(const std::runtime_error &) {
                write_vtp(system,"failure_liquid.vtp");
                throw;
            }
        }
        rates=system.rates();
        bool released=was_equilibrating&&!system.equilibration_active;
        if(released)std::cout<<"Fluid equilibration complete at time="<<system.time<<"; DEM coupling enabled\n";
        ke_reached=observe_ke();
        bool final=system.time>=a.end_time||system.step-initial_step==a.max_steps||ke_reached;
        if(a.checkpoint_every && (system.step%a.checkpoint_every==0||final||released))system.save("dem_checkpoint.dat");
        const bool frame=a.dump_every && (system.step%a.dump_every==0||final||released);
        if(system.step%a.record_every==0||final||released||(calibration&&frame))record(rates);
        if(frame)dump();
    }
    system.save("dem_checkpoint.dat");
    const bool time_reached=system.time>=a.end_time;
    const char *reason=time_reached?"time_limit":(ke_reached?"kinetic_energy":"step_limit");
    if(ke_reached)std::cout<<"Kinetic energy target reached at time="<<system.time
                           <<"; KE/peak="<<current_ke/peak_ke<<'\n';
    std::ofstream meta("dem_run.json");
    meta<<std::setprecision(17)<<"{\"physics_validated\":false,\"completed_time\":"<<(time_reached?"true":"false")
        <<",\"termination_reason\":"<<std::quoted(reason)<<",\"time\":"<<system.time
        <<",\"step\":"<<system.step<<",\"threads\":"<<threads
        <<",\"wall_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if(ke_target_fraction){
        meta<<",\"ke_target_fraction\":"<<*ke_target_fraction<<",\"ke_hold_time\":"<<ke_hold_time.value_or(0.)
            <<",\"peak_kinetic_energy\":"<<peak_ke<<",\"final_ke_fraction\":";
        if(peak_ke>0)meta<<current_ke/peak_ke;else meta<<"null";
    }
    meta<<"}\n";
    if(!meta)throw std::runtime_error("Failed DEM metadata output");
    return time_reached||ke_reached?0:2;
} catch(const std::exception &e){std::cerr<<"pf_dem: "<<e.what()<<'\n';return 1;}
