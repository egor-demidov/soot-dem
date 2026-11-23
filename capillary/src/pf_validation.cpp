// Gas-free PF-SPH validation. Native SPHinXsys storage, neighbor search and VTP I/O;
// separate published WCSPH operators, without upstream multiphase/CSF dynamics.
#include "sphinxsys.h"
#include "pf_liquid.h"
#include "pf_bridge.h"
#include "pf_parallel.h"
#include "pf_neighbors.h"
#include "pf_restart.h"
#include "pf_snapshot.h"
#include <tbb/global_control.h>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <string>
using namespace SPH;
namespace pf=pf_liquid;
struct Options {
    pf::Parameters p;
    std::string kind="droplet";
    double radius=.6,end=5.,dt_scale=1.,snapshot_interval=0.,checkpoint_interval=0.;
    double sphere_radius=.2,sphere_density=2.,sphere_vx=.0,sphere_omega_z=.0;
    double bridge_gap=.24,contact_stiffness=0.;
    bool fixed_spheres=true;
    std::string resume, resume_snapshot;
    std::string neighbor_mode="verlet";
    double neighbor_skin_factor=.25;
    int steps=10000000,dump=500,record=50,threads=2;
};
Options parse(int argc,char **argv) {
    Options o;
    for(int i=1;i<argc;++i) {
        std::string k=argv[i];
        if(k=="--help") {
            std::cout<<"--case droplet|sessile|immersed|bridge --bridge-gap G --fixed-spheres 0|1 --contact-stiffness K --sphere-radius R --sphere-density RHO --sphere-vx V --sphere-omega-z W --dx D --radius R --sff S --sfs S --sound-speed C --viscosity MU --end-time T --steps N --dt-scale F --dump-every N --record-every N --threads N --snapshot-interval T --checkpoint-interval T --resume FILE --resume-snapshot CSV --neighbor-mode native|cells|verlet --neighbor-skin-factor F\n";
            std::exit(0);
        }
        if(i+1==argc) throw std::runtime_error("Missing value for "+k);
        std::string v=argv[++i];
        if(k=="--case") o.kind=v;
        else if(k=="--dx") o.p.dx=std::stod(v);
        else if(k=="--radius") o.radius=std::stod(v);
        else if(k=="--bridge-gap") o.bridge_gap=std::stod(v);
        else if(k=="--contact-stiffness") o.contact_stiffness=std::stod(v);
        else if(k=="--fixed-spheres") {
            if(v!="0" && v!="1") throw std::runtime_error("--fixed-spheres must be 0 or 1");
            o.fixed_spheres=v=="1";
        }
        else if(k=="--sphere-radius") o.sphere_radius=std::stod(v);
        else if(k=="--sphere-density") o.sphere_density=std::stod(v);
        else if(k=="--sphere-vx") o.sphere_vx=std::stod(v);
        else if(k=="--sphere-omega-z") o.sphere_omega_z=std::stod(v);
        else if(k=="--sff") o.p.sff=std::stod(v);
        else if(k=="--sfs") o.p.sfs=std::stod(v);
        else if(k=="--sound-speed") o.p.sound=std::stod(v);
        else if(k=="--viscosity") o.p.mu=std::stod(v);
        else if(k=="--end-time") o.end=std::stod(v);
        else if(k=="--steps") o.steps=std::stoi(v);
        else if(k=="--dt-scale") o.dt_scale=std::stod(v);
        else if(k=="--dump-every") o.dump=std::stoi(v);
        else if(k=="--record-every") o.record=std::stoi(v);
        else if(k=="--threads") o.threads=std::stoi(v);
        else if(k=="--snapshot-interval") o.snapshot_interval=std::stod(v);
        else if(k=="--checkpoint-interval") o.checkpoint_interval=std::stod(v);
        else if(k=="--neighbor-mode") o.neighbor_mode=v;
        else if(k=="--neighbor-skin-factor") o.neighbor_skin_factor=std::stod(v);
        else if(k=="--resume") o.resume=v;
        else if(k=="--resume-snapshot") o.resume_snapshot=v;
        else throw std::runtime_error("Unknown option "+k);
    }
    o.p.validate();
    if((o.kind!="droplet"&&o.kind!="sessile"&&o.kind!="immersed"&&o.kind!="bridge")||!(o.radius>o.p.dx&&std::isfinite(o.radius))||
       !(o.end>0&&std::isfinite(o.end))||!(o.dt_scale>0&&o.dt_scale<=1)||
       o.steps<1||o.dump<0||o.record<1||o.threads<1||
       (o.neighbor_mode!="native"&&o.neighbor_mode!="cells"&&o.neighbor_mode!="verlet") ||
       !(o.neighbor_skin_factor>=0.&&o.neighbor_skin_factor<=1.&&std::isfinite(o.neighbor_skin_factor)) ||
       !(o.snapshot_interval>=0&&std::isfinite(o.snapshot_interval))||
       !(o.checkpoint_interval>=0&&std::isfinite(o.checkpoint_interval))||
       (!o.resume.empty()&&!o.resume_snapshot.empty()) ||
       !(o.sphere_radius>0&&std::isfinite(o.sphere_radius)) ||
       !(o.sphere_density>0&&std::isfinite(o.sphere_density)) ||
       !std::isfinite(o.sphere_vx) || !std::isfinite(o.sphere_omega_z) ||
       !(o.bridge_gap>0&&std::isfinite(o.bridge_gap)) ||
       !(o.contact_stiffness>=0&&std::isfinite(o.contact_stiffness)) ||
       (o.kind=="immersed" && (o.sphere_radius+2*o.p.dx>=o.radius || !o.resume_snapshot.empty())) ||
       (o.kind=="bridge" && (o.sphere_radius<o.p.dx || !o.resume_snapshot.empty()))) throw std::runtime_error("Invalid PF run options");
    return o;
}
class Region:public Shape {
    static pf::Vector vec_bridge(const Vecd &x) {return {x[0],x[1],x[2]};}
    double radius_,extent_,bottom_,depth_,inner_;
    bool solid_,sessile_,spherical_solid_,bridge_=false;
    pf_bridge::Geometry bridge_geometry_;
public:
    Region(std::string name,double radius,double extent,double bottom,double depth,
           bool solid,bool sessile,double inner=0.,bool spherical_solid=false)
        :Shape(name),radius_(radius),extent_(extent),bottom_(bottom),depth_(depth),
         inner_(inner),solid_(solid),sessile_(sessile),spherical_solid_(spherical_solid){}
    Region(std::string name,const pf_bridge::Geometry &geometry,bool solid)
      :Shape(name),radius_(geometry.liquid_radius),extent_(geometry.solid_x_extent()),
       bottom_(-geometry.solid_x_extent()),depth_(geometry.sphere_radius),inner_(0.),
       solid_(solid),sessile_(false),spherical_solid_(false),bridge_(true),bridge_geometry_(geometry) {
        bridge_geometry_.validate();
    }
    bool checkContain(const Vecd &x,bool=true) override {
        if(bridge_) return solid_?bridge_geometry_.in_solid(vec_bridge(x)):
                                 bridge_geometry_.in_liquid(vec_bridge(x));
        if(solid_) {
            if(spherical_solid_) return x.squaredNorm()<radius_*radius_;
            return std::abs(x[0])<extent_&&std::abs(x[1])<extent_&&x[2]<0&&x[2]>-depth_;
        }
        return x.squaredNorm()<radius_*radius_ && x.squaredNorm()>inner_*inner_ &&
               (!sessile_||x[2]>0);
    }
    Vecd findClosestPoint(const Vecd &x) override {
        if(bridge_) {
            if(solid_) {
                const auto origin=bridge_geometry_.center(x[0]<0?0:1);
                Vecd v(x[0]-origin[0],x[1],x[2]);
                if(v.norm()==0) return Vecd(origin[0]+bridge_geometry_.sphere_radius,0,0);
                return Vecd(origin[0]+bridge_geometry_.sphere_radius*v[0]/v.norm(),
                            bridge_geometry_.sphere_radius*v[1]/v.norm(),
                            bridge_geometry_.sphere_radius*v[2]/v.norm());
            }
            const double transverse=std::hypot(x[1],x[2]);
            if(transverse==0) return Vecd(x[0],bridge_geometry_.liquid_radius,0);
            const double factor=bridge_geometry_.liquid_radius/transverse;
            return Vecd(x[0],x[1]*factor,x[2]*factor);
        }
        if(solid_&&!spherical_solid_) return Vecd(x[0],x[1],0.);
        const double r=x.norm();
        if(r==0) return Vecd(radius_,0,0);
        const double boundary=(inner_>0 && std::abs(r-inner_)<std::abs(r-radius_))?inner_:radius_;
        return boundary*x/r;
    }
protected:
    BoundingBox findBounds() override {
        if(bridge_) {
            const double length=solid_?bridge_geometry_.solid_x_extent():bridge_geometry_.half_cylinder_length();
            const double r=solid_?bridge_geometry_.sphere_radius:bridge_geometry_.liquid_radius;
            return BoundingBox(Vecd(-length,-r,-r),Vecd(length,r,r));
        }
        if(spherical_solid_) return BoundingBox(Vecd(-radius_,-radius_,-radius_),Vecd(radius_,radius_,radius_));
        return BoundingBox(Vecd(-extent_,-extent_,bottom_),Vecd(extent_,extent_,extent_));
    }
};
pf::Vector vec(const Vecd &x) {return {x[0],x[1],x[2]};}
Vecd vec(const pf::Vector &x) {return Vecd(x[0],x[1],x[2]);}
int run(const Options &o) {
    tbb::global_control threads(tbb::global_control::max_allowed_parallelism,o.threads);
    const bool sessile=o.kind=="sessile",immersed=o.kind=="immersed",bridge=o.kind=="bridge";
    const bool rigid=immersed||bridge;
    pf_bridge::Geometry geometry;
    geometry.sphere_radius=o.sphere_radius;geometry.gap=o.bridge_gap;geometry.liquid_radius=o.radius;
    if(bridge) geometry.validate();
    double H=o.p.support();
    double domain_radius=bridge?std::max(geometry.solid_x_extent(),geometry.half_cylinder_length()):o.radius;
    double extent=std::ceil((domain_radius+2*H+1)/o.p.dx)*o.p.dx;
    BoundingBox bounds(Vecd(-extent,-extent,-extent),Vecd(extent,extent,extent));
    SPHSystem system(bounds,o.p.dx);
    GlobalStaticVariables::physical_time_=0.;
    IOEnvironment io(system);
    auto fluid_shape=bridge?makeShared<Region>("Liquid",geometry,false):
        makeShared<Region>("Liquid",o.radius,extent,-extent,H,false,sessile,
                           immersed?o.sphere_radius+.5*o.p.dx:0.);
    FluidBody liquid(system,fluid_shape);
    liquid.defineAdaptationRatios(2.); // Wendland support 2h=4dx, BEFORE cell-list construction.
    liquid.defineParticlesAndMaterial<FluidParticles,WeaklyCompressibleFluid>(o.p.rho0,o.p.sound,o.p.mu);
    liquid.generateParticles<ParticleGeneratorLattice>();
    auto &w=dynamic_cast<FluidParticles&>(liquid.getBaseParticles());
    std::unique_ptr<SolidBody> wall;
    std::unique_ptr<ContactRelation> contact;
    BaseParticles *b=nullptr;
    if(sessile||rigid) {
        // Fixed plate, single immersed sphere, or two distinct moving/fixed spheres.
        double plate=std::ceil((o.radius+2*H)/o.p.dx)*o.p.dx;
        auto solid_shape=bridge?makeShared<Region>("Solid",geometry,true):
            makeShared<Region>("Solid",immersed?o.sphere_radius:o.radius,
                               plate,-H,H,true,true,0.,immersed);
        wall=std::make_unique<SolidBody>(system,solid_shape);
        wall->defineAdaptationRatios(2.);
        wall->defineParticlesAndMaterial<SolidParticles,Solid>();
        wall->generateParticles<ParticleGeneratorLattice>();
        b=&wall->getBaseParticles();
        contact=std::make_unique<ContactRelation>(liquid,RealBodyVector{wall.get()});
    }
    InnerRelation inner(liquid);
    system.initializeSystemCellLinkedLists(); system.initializeSystemConfigurations();
    const size_t nf=w.total_real_particles_,ns=b?b->total_real_particles_:0;
    if(nf==0 || (rigid && ns==0)) throw std::runtime_error("No liquid/DEM solid samples");
    std::vector<pf::Sample> state; state.reserve(nf+ns);
    for(size_t i=0;i<nf;++i) state.push_back({vec(w.pos_[i]),{},o.p.rho0,o.p.rho0*std::pow(o.p.dx,3),false});
    for(size_t i=0;i<ns;++i) state.push_back({vec(b->pos_[i]),{},o.p.rho0,o.p.rho0*std::pow(o.p.dx,3),true});
    std::vector<pf_rigid::Sphere> spheres;
    std::vector<pf_rigid::Site> sites;
    if(rigid) {
        const size_t body_count=bridge?2:1;
        for(size_t j=0;j<body_count;++j) {
            pf_rigid::Sphere sphere;
            sphere.radius=o.sphere_radius;
            sphere.mass=(4./3.)*Pi*std::pow(o.sphere_radius,3)*o.sphere_density;
            sphere.center=bridge?geometry.center(j):pf::Vector{};
            sphere.fixed=bridge && o.fixed_spheres;
            sphere.velocity={bridge?(j==0?-o.sphere_vx:o.sphere_vx):o.sphere_vx,0.,0.};
            sphere.angular_momentum={0.,0.,sphere.inertia()*o.sphere_omega_z};
            if(sphere.fixed) {sphere.velocity={};sphere.angular_momentum={};}
            spheres.push_back(sphere);
        }
        std::array<size_t,2> counts{0,0};
        for(size_t i=nf;i<state.size();++i) {
            size_t j=bridge?(state[i].x[0]<0?0:1):0;
            sites.push_back({i,j,pf::subtract(state[i].x,spheres[j].center)});
            counts[j]++;
        }
        if(bridge && (counts[0]==0||counts[1]==0))
            throw std::runtime_error("Bridge solid sampling omitted one sphere; reduce dx");
        pf_rigid::validate_sites(state,spheres,sites);
        pf_rigid::sync_sites(state,spheres,sites);
    }
    pf::Vector wallP{},wallL{},referenceP{},referenceL{};
    if(rigid) {
        referenceP=pf_rigid::total_momentum(state,spheres);
        referenceL=pf_rigid::total_momentum(state,spheres,true);
    }
    double time=0,lastdt=0; long long step=0;
    bool legacy_rebased=false;
    if(!o.resume.empty()) {
        auto c=pf_restart::read_file(o.resume);
        pf_restart::check_compatibility(c,o.kind,o.p,o.radius,state,nf,&spheres,
                                        rigid?&sites:nullptr,bridge?o.bridge_gap:0.,
                                        bridge?o.contact_stiffness:0.);
        state=std::move(c.state); time=c.time; step=c.step; lastdt=c.lastdt;
        wallP=c.wallP; wallL=c.wallL;
        referenceP=c.referenceP;referenceL=c.referenceL;
        if(rigid) spheres=std::move(c.spheres);
        std::cout<<"Resuming complete PF checkpoint: "<<o.resume<<" at t="<<time<<" step="<<step<<"\n";
    } else if(!o.resume_snapshot.empty()) {
        auto recovered=pf_snapshot::load_file(o.resume_snapshot,state,nf);
        state=std::move(recovered.second);
        time=recovered.first.time; step=recovered.first.step;
        wallP=pf::scale(pf_snapshot::fluid_momentum(state,nf,false),-1.);
        wallL=pf::scale(pf_snapshot::fluid_momentum(state,nf,true),-1.);
        legacy_rebased=true;
        std::cerr<<"WARNING: recovered fluid state at t="<<time<<" step="<<step
                 <<"; historical wall impulse missing from legacy CSV, momentum ledger rebased here. "
                 <<"Do not compare absolute cumulative impulse with the uninterrupted run.\n";
    }
    if(time>=o.end) throw std::runtime_error("PF resume time is not before requested end-time");
    GlobalStaticVariables::physical_time_=time;
    auto sync=[&](const auto &s) {
        for(size_t i=0;i<s.size();++i) {
            for(int a=0;a<3;++a) if(std::abs(s[i].x[a])>=extent-H)
                throw std::runtime_error("Liquid reached neighbor-search domain margin; enlarge domain");
            if(i<nf) {
                w.pos_[i]=vec(s[i].x); w.vel_[i]=vec(s[i].v); w.rho_[i]=s[i].rho;
                w.Vol_[i]=s[i].mass/s[i].rho; w.p_[i]=pf::pressure(s[i].rho,o.p);
            } else if(rigid) {
                b->pos_[i-nf]=vec(s[i].x);
                b->vel_[i-nf]=vec(s[i].v);
            }
        }
        liquid.setNewlyUpdated(); // Direct state updates must notify native VTP I/O.
        if(rigid) wall->setNewlyUpdated();
    };
    struct Timings {
        double graph_s=0,force_s=0;
        size_t graph_calls=0,force_calls=0,pairs_seen=0,graph_builds=0,graph_reuses=0;
    };
    Timings perf;
    const auto wallclock_begin=std::chrono::steady_clock::now();
    auto elapsed=[](auto begin) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    };
    // The native graph is retained as a reference, but it rebuilds both cell
    // lists/configurations at every midpoint and endpoint evaluation.
    // "cells" builds an independent cell graph each call; "verlet" reuses
    // candidate pairs until any sample (including moving DEM sites) moves
    // by skin/2 from the last rebuild. All interactions remain cutoff at H.
    pf_neighbors::Verlet cached(nf,H,H*o.neighbor_skin_factor);
    pf::Pairs native_pairs;
    bool audited=false;
    auto graph=[&](const auto &s)->const pf::Pairs& {
        const auto started=std::chrono::steady_clock::now();
        for(const auto &sample:s) for(int axis=0;axis<3;++axis)
            if(!std::isfinite(sample.x[axis]) || std::abs(sample.x[axis])>=extent-H)
                throw std::runtime_error("PF sample reached neighbor-search domain margin");
        if(o.neighbor_mode=="native") {
            sync(s);
            liquid.updateCellLinkedList(); inner.updateConfiguration();
            if(contact) {
                if(rigid) wall->updateCellLinkedList();
                contact->updateConfiguration();
            }
            native_pairs.clear(); native_pairs.reserve(nf*140);
            for(size_t i=0;i<nf;++i) {
                const auto &nb=inner.inner_configuration_[i];
                for(size_t k=0;k<nb.current_size_;++k)
                    if(nb.j_[k]>i) native_pairs.emplace_back(i,nb.j_[k]);
                if(contact) {
                    const auto &nbc=contact->contact_configuration_[0][i];
                    for(size_t k=0;k<nbc.current_size_;++k)
                        native_pairs.emplace_back(i,nf+nbc.j_[k]);
                }
            }
            std::sort(native_pairs.begin(),native_pairs.end());
            ++perf.graph_builds;
        } else if(o.neighbor_mode=="cells") {
            native_pairs=pf_neighbors::build_cells(s,nf,H);
            ++perf.graph_builds;
        } else {
            const size_t prior=cached.builds();
            cached.graph(s);
            if(cached.builds()!=prior) ++perf.graph_builds;
            else ++perf.graph_reuses;
        }
        const pf::Pairs &pairs=o.neighbor_mode=="verlet"?cached.candidates():native_pairs;
        // The legacy small-system audit is deliberately retained for both
        // native and cached implementations, using the true physical cutoff.
        if(!audited && s.size()<4000) {
            const auto expected=pf_neighbors::brute_force(s,nf,H);
            pf::Pairs actual;
            for(auto pair:pairs) {
                const auto d=pf::subtract(s[pair.first].x,s[pair.second].x);
                if(pf::dot(d,d)<H*H) actual.push_back(pair);
            }
            if(actual!=expected || std::adjacent_find(pairs.begin(),pairs.end())!=pairs.end())
                throw std::runtime_error("PF neighbor graph missed or duplicated a pair");
            std::cout<<"Neighbor audit passed ("<<o.neighbor_mode<<"): "<<actual.size()<<" active pairs\n";
        }
        audited=true;
        ++perf.graph_calls; perf.graph_s+=elapsed(started); perf.pairs_seen+=pairs.size();
        return pairs;
    };
    w.addVariableToWrite<Real>("Pressure");
    BodyStatesRecordingToVtp output(io,system.real_bodies_);
    std::ofstream perf_log("pf_perf.csv"); perf_log<<std::setprecision(17)
        <<"step,time,wall_seconds,graph_seconds,evaluate_seconds,graph_calls,evaluate_calls,total_pairs,avg_graph_ms,avg_evaluate_ms,neighbor_builds,neighbor_reuses\n";
    std::ofstream hist("pf_history.csv"); hist<<std::setprecision(17);
    hist<<"step,time,mass,volume,kinetic,pair_energy,max_speed,max_density_error,linear_error,angular_error,px,py,pz,lx,ly,lz,wall_fx,wall_fy,wall_fz,wall_tx,wall_ty,wall_tz,volume_radius,virial_radius,virial_pressure,dt,core_eos_pressure,core_mean_density,core_particle_count\n";
    auto diagnostics=[&](const pf::Rates &rates) {
        double mass=0,volume=0,ke=0,speed=0,rhoerr=0; pf::Vector P{},L{},center{};
        for(size_t i=0;i<nf;++i) {
            const auto &a=state[i]; mass+=a.mass; volume+=a.mass/a.rho;
            ke+=.5*a.mass*pf::dot(a.v,a.v); speed=std::max(speed,std::sqrt(pf::dot(a.v,a.v)));
            rhoerr=std::max(rhoerr,std::abs(a.rho/o.p.rho0-1));
            P=pf::add(P,pf::scale(a.v,a.mass)); L=pf::add(L,pf::cross(a.x,pf::scale(a.v,a.mass)));
            center=pf::add(center,pf::scale(a.x,a.mass));
        }
        center=pf::scale(center,1/mass);
        double R=std::cbrt(3*volume/(4*Pi)),rv=R-2*H,virial=0;
        double core_eos=0., core_rho=0.;size_t core_count=0;
        if(!sessile && !rigid && rv>o.p.dx) {
            for(size_t i=0;i<nf;++i) if(pf::dot(pf::subtract(state[i].x,center),pf::subtract(state[i].x,center))<rv*rv) {
                virial+=rates.conservative_virial[i];
                core_eos+=pf::pressure(state[i].rho,o.p);
                core_rho+=state[i].rho;
                ++core_count;
            }
            virial/=4*Pi*rv*rv*rv; // 3*V, eq.3.2 directed-half sum
            if(core_count>0) {core_eos/=core_count;core_rho/=core_count;}
        } else {rv=std::numeric_limits<double>::quiet_NaN(); virial=rv;}
        if(!core_count) {
            core_eos=std::numeric_limits<double>::quiet_NaN();
            core_rho=std::numeric_limits<double>::quiet_NaN();
        }
        pf::RigidLoad wall_load;
        for(size_t i=nf;i<state.size();++i) wall_load.add_sample(state[i].x,rates.force[i],{});
        auto conservedP=rigid?pf::add(pf_rigid::total_momentum(state,spheres),wallP):pf::add(P,wallP);
        auto conservedL=rigid?pf::add(pf_rigid::total_momentum(state,spheres,true),wallL):pf::add(L,wallL);
        double ep=std::sqrt(pf::dot(pf::subtract(conservedP,referenceP),pf::subtract(conservedP,referenceP)));
        double el=std::sqrt(pf::dot(pf::subtract(conservedL,referenceL),pf::subtract(conservedL,referenceL)));
        hist<<step<<','<<time<<','<<mass<<','<<volume<<','<<ke<<','<<rates.pair_energy<<','<<speed<<','<<rhoerr<<','<<ep<<','<<el;
        for(auto v:{P,L,wall_load.force,wall_load.torque}) for(double a:v) hist<<','<<a;
        hist<<','<<R<<','<<rv<<','<<virial<<','<<lastdt<<','<<core_eos<<','<<core_rho<<','<<core_count<<'\n'; hist.flush();
        perf_log<<step<<','<<time<<','<<elapsed(wallclock_begin)<<','<<perf.graph_s<<','<<perf.force_s
                <<','<<perf.graph_calls<<','<<perf.force_calls<<','<<perf.pairs_seen
                <<','<<(perf.graph_calls?1000*perf.graph_s/perf.graph_calls:0.)
                <<','<<(perf.force_calls?1000*perf.force_s/perf.force_calls:0.)
                <<','<<perf.graph_builds<<','<<perf.graph_reuses<<'\n';
        perf_log.flush();
        std::cout<<"PF liquid: step="<<step<<" time="<<time<<" maxSpeed="<<speed<<" densityError="<<rhoerr<<" linearError="<<ep<<" angularError="<<el<<'\n';
        if(speed/o.p.sound>.2||rhoerr>.2) throw std::runtime_error("PF state outside validation limits; inspect history before changing inputs");
    };
    auto write_samples=[&](const std::string &filename) {
        std::ofstream positions(filename); positions<<std::setprecision(17)<<"x,y,z,rho,vx,vy,vz\n";
        for(size_t i=0;i<nf;++i) {
            const auto &a=state[i]; positions<<a.x[0]<<','<<a.x[1]<<','<<a.x[2]<<','<<a.rho;
            for(double v:a.v) positions<<','<<v;
            positions<<'\n';
        }
        if(!positions) throw std::runtime_error("Cannot write liquid snapshot "+filename);
    };
    std::ofstream snapshot_index;
    double next_snapshot=o.snapshot_interval;
    auto snapshot=[&]() {
        std::string filename="pf_snapshot_"+std::to_string(step)+".csv";
        write_samples(filename);
        snapshot_index<<step<<','<<std::setprecision(17)<<time<<','<<filename<<'\n';
        snapshot_index.flush();
        if(!snapshot_index) throw std::runtime_error("Cannot write snapshot index");
    };
    if(o.snapshot_interval>0) {
        snapshot_index.open("pf_snapshots.csv"); snapshot_index<<"step,time,file\n"; snapshot();
    }
    auto evaluator=[&](const auto &samples,const auto &pairs,const auto &parameters) {
        auto started=std::chrono::steady_clock::now();
        auto out=pf_parallel::evaluate(samples,pairs,parameters,static_cast<unsigned>(o.threads));
        ++perf.force_calls; perf.force_s+=elapsed(started);
        return out;
    };
    std::ofstream rigid_log;
    std::ofstream bridge_log;
    if(rigid) {
        rigid_log.open("pf_rigid.csv");
        rigid_log<<std::setprecision(17)
                 <<"step,time,cx,cy,cz,vx,vy,vz,wx,wy,wz,lx,ly,lz,fx,fy,fz,tx,ty,tz,contact_energy,body_id,gap\n";
    }
    if(bridge) {
        bridge_log.open("pf_bridge.csv");
        bridge_log<<std::setprecision(17)
            <<"step,time,separation,gap,neck_particles,wetting_left,wetting_right,liquid_connected,min_liquid_clearance,left_fx,right_fx,"
            <<"left_pf_fx,right_pf_fx,left_pressure_viscous_fx,right_pressure_viscous_fx,contact_energy\n";
    }
    auto write_rigid=[&](const pf::Rates &r) {
        if(!rigid) return;
        auto loads=pf_rigid::total_loads(state,r,spheres,sites,o.contact_stiffness);
        const auto contact_load=pf_rigid::sphere_contacts(spheres,o.contact_stiffness);
        const double gap=spheres.size()==2?
            std::sqrt(pf::dot(pf::subtract(spheres[0].center,spheres[1].center),
                              pf::subtract(spheres[0].center,spheres[1].center)))
            -spheres[0].radius-spheres[1].radius : 0.;
        for(size_t j=0;j<spheres.size();++j) {
            const auto &a=spheres[j];const auto &f=loads[j];
            rigid_log<<step<<','<<time;
            for(auto v:{a.center,a.velocity,a.omega(),a.angular_momentum,f.force,f.torque})
                for(auto x:v) rigid_log<<','<<x;
            rigid_log<<','<<contact_load.energy<<','<<j<<','<<gap<<'\n';
        }
        rigid_log.flush();
        if(bridge) {
            auto pairs=graph(state);
            auto d=pf_bridge::diagnostics(state,pairs,o.p,geometry,spheres,sites);
            bridge_log<<step<<','<<time<<','<<gap+2*o.sphere_radius<<','<<gap<<','
                      <<d.neck_particles<<','<<d.wetting_left<<','<<d.wetting_right<<','
                      <<(d.liquid_connected?1:0)<<','<<d.liquid_min_clearance<<','
                      <<loads[0].force[0]<<','<<loads[1].force[0]<<','
                      <<d.pf_solid[0].force[0]<<','<<d.pf_solid[1].force[0]<<','
                      <<loads[0].force[0]-d.pf_solid[0].force[0]-contact_load.loads[0].force[0]<<','
                      <<loads[1].force[0]-d.pf_solid[1].force[0]-contact_load.loads[1].force[0]<<','
                      <<contact_load.energy<<'\n';
            bridge_log.flush();
            if(d.liquid_min_clearance<-o.p.dx*.25)
                throw std::runtime_error("Bridge liquid penetrated sphere by more than dx/4");
        }
    };
    auto checkpoint=[&]() {
        pf_restart::Checkpoint c;
        c.kind=o.kind; c.parameters=o.p; c.radius=o.radius;
        c.time=time; c.lastdt=lastdt; c.step=step;
        c.bridge_gap=bridge?o.bridge_gap:0.; c.contact_stiffness=bridge?o.contact_stiffness:0.;
        c.nf=nf; c.ns=ns; c.wallP=wallP; c.wallL=wallL; c.state=state;
        c.referenceP=referenceP;c.referenceL=referenceL;c.spheres=spheres;
        pf_restart::write_atomic("pf_checkpoint.dat",c);
    };
    double next_checkpoint=o.checkpoint_interval>0?
        (std::floor(time/o.checkpoint_interval)+1)*o.checkpoint_interval:0.;
    auto rates=evaluator(state,graph(state),o.p);
    std::cout<<"Gas-free PF-SPH: fluid="<<nf<<" solid="<<ns<<" H/dx=4 sff="<<o.p.sff<<" sfs="<<o.p.sfs<<" (strengths require calibration)\n";
    diagnostics(rates);write_rigid(rates); if(o.dump) {sync(state);output.writeToFile(step);}
    if(o.checkpoint_interval>0) checkpoint();
    while(time<o.end && step<o.steps) {
        double dt=std::min(o.end-time,o.dt_scale*(rigid ?
            pf_rigid::timestep(state,spheres,sites,rates,o.p,o.contact_stiffness) : pf::timestep(state,rates,o.p)));
        if(!(dt>0&&time+dt>time)) throw std::runtime_error("PF timestep underflow");
        auto midpoint=rigid ? pf_rigid::advance(state,spheres,sites,o.p,dt,graph,&rates,evaluator,o.contact_stiffness)
                            : pf::advance(state,o.p,dt,graph,&rates,evaluator);
        if(!rigid) for(size_t i=nf;i<state.size();++i) {
            wallP=pf::add(wallP,pf::scale(midpoint.force[i],dt));
            wallL=pf::add(wallL,pf::scale(pf::cross(state[i].x,midpoint.force[i]),dt));
        }
        // A fixed DEM body exchanges impulse with its external support. Include
        // this reaction in the conserved fluid + movable-body momentum ledger.
        if(rigid) for(const auto &site:sites) if(spheres[site.body].fixed) {
            const auto &f=midpoint.force[site.particle];
            wallP=pf::add(wallP,pf::scale(f,dt));
            wallL=pf::add(wallL,pf::scale(pf::cross(state[site.particle].x,f),dt));
        }
        time=dt==o.end-time?o.end:time+dt; ++step; lastdt=dt;
        rates=evaluator(state,graph(state),o.p); // current acceleration for next restriction
        GlobalStaticVariables::physical_time_=time;
        if(step%o.record==0||time>=o.end||step==o.steps) {diagnostics(rates);write_rigid(rates);}
        if(o.snapshot_interval>0&&(time>=next_snapshot||time>=o.end||step==o.steps)) {
            snapshot();
            next_snapshot=(std::floor(time/o.snapshot_interval)+1)*o.snapshot_interval;
        }
        if(o.dump&&(step%o.dump==0||time>=o.end||step==o.steps)) {
            sync(state);output.writeToFile(step);
        }
        if(o.checkpoint_interval>0&&(time>=next_checkpoint||time>=o.end||step==o.steps)) {
            checkpoint();
            next_checkpoint=(std::floor(time/o.checkpoint_interval)+1)*o.checkpoint_interval;
        }
    }
    write_samples("pf_final_liquid.csv");
    std::ofstream summary("pf_run.json"); summary<<std::setprecision(17)
        <<"{\"model\":\"Whebell2025-PF-WCSPH\",\"case\":\""<<o.kind<<"\",\"completed_time\":"<<(time>=o.end?"true":"false")
        <<",\"time\":"<<time<<",\"steps\":"<<step<<",\"dx\":"<<o.p.dx<<",\"radius\":"<<o.radius
        <<",\"resumed_from\":\""<<(o.resume.empty()?o.resume_snapshot:o.resume)
        <<"\",\"legacy_momentum_rebased\":"<<(legacy_rebased?"true":"false")
        <<",\"sphere_radius\":"<<(rigid?o.sphere_radius:0.)
        <<",\"sphere_density\":"<<(rigid?o.sphere_density:0.)
        <<",\"sphere_count\":"<<spheres.size()
        <<",\"bridge_gap\":"<<(bridge?o.bridge_gap:0.)
        <<",\"contact_stiffness\":"<<(bridge?o.contact_stiffness:0.)
        <<",\"fixed_spheres\":"<<(bridge&&o.fixed_spheres?"true":"false")
        <<",\"sff\":"<<o.p.sff<<",\"sfs\":"<<o.p.sfs<<",\"sound_speed\":"<<o.p.sound<<",\"viscosity\":"<<o.p.mu<<",\"snapshot_interval\":"<<o.snapshot_interval
        <<",\"neighbor_mode\":\""<<o.neighbor_mode<<"\""
        <<",\"neighbor_skin_factor\":"<<o.neighbor_skin_factor
        <<",\"neighbor_builds\":"<<perf.graph_builds
        <<",\"neighbor_reuses\":"<<perf.graph_reuses
        <<",\"rho0\":"<<o.p.rho0<<",\"fluid_samples\":"<<nf<<",\"solid_samples\":"<<ns<<"}\n";
    if(time<o.end) {std::cerr<<"Step limit reached before end time\n"; return 2;}
    return 0;
}
int main(int argc,char **argv) {
    try {return run(parse(argc,argv));}
    catch(const std::exception &e) {std::cerr<<"pf_capillary_validation: "<<e.what()<<'\n';return 1;}
}
