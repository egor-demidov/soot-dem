#pragma once
// N-body DEM / gas-free PF-SPH coupling. Quadrature samples carry no DEM mass.
// Forces are evaluated once per fluid pair; their reciprocal solid endpoint
// supplies both the body force and its moment. No global momentum projection.
#include "pf_neighbors.h"
#include "pf_rigid.h"
#include "pf_mechanics.h"
#include <functional>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

namespace pf_dem {
using namespace pf_rigid;
struct ContactModel {
    double stiffness=100., damping=0., friction=0., slip_speed=.01;
    void validate() const {
        for(double v:{stiffness,damping,friction})
            if(!(v>=0&&std::isfinite(v))) throw std::invalid_argument("Invalid DEM contact coefficient");
        if(!(slip_speed>0&&std::isfinite(slip_speed))) throw std::invalid_argument("Invalid friction regularization speed");
    }
};
struct ContactResult {
    std::vector<Load> loads;
    double energy=0., power=0.;
    double dt=std::numeric_limits<double>::infinity();
};
// Reversible Hertz normal elasticity + optional unilateral normal dashpot and
// regularized kinetic Coulomb friction. No history spring or bonded neck model.
inline ContactResult contacts(const std::vector<Sphere> &b,const ContactModel &m,double dx) {
    m.validate(); ContactResult out;out.loads.resize(b.size());
    for(size_t i=0;i<b.size();++i) for(size_t j=i+1;j<b.size();++j) {
        auto d=subtract(b[j].center,b[i].center);double r=std::sqrt(dot(d,d));
        double overlap=b[i].radius+b[j].radius-r;
        if(overlap<=0 || (b[i].fixed&&b[j].fixed)) continue;
        if(r==0) throw std::runtime_error("Coincident DEM centers");
        const auto n=scale(d,1/r);
        // One common point, even for overlapping spheres. Distinct nominal
        // surface points would introduce a spurious couple under friction.
        auto ai=scale(n,b[i].radius-.5*overlap);
        auto aj=subtract(add(b[i].center,ai),b[j].center);
        const auto vi=b[i].fixed?Vector{}:add(b[i].velocity,cross(b[i].omega(),ai));
        const auto vj=b[j].fixed?Vector{}:add(b[j].velocity,cross(b[j].omega(),aj));
        const auto dv=subtract(vi,vj);const double vn=dot(dv,n);
        const auto vt=subtract(dv,scale(n,vn));double speed=std::sqrt(dot(vt,vt));
        const double elastic=m.stiffness*std::pow(overlap,1.5);
        const double fn=std::max(0.,elastic+m.damping*vn);
        const auto ft=scale(vt,-m.friction*fn/std::max(m.slip_speed,speed));
        const auto fi=add(scale(n,-fn),ft),fj=scale(fi,-1);
        out.loads[i].force=add(out.loads[i].force,fi);
        out.loads[j].force=add(out.loads[j].force,fj);
        out.loads[i].torque=add(out.loads[i].torque,cross(ai,fi));
        out.loads[j].torque=add(out.loads[j].torque,cross(aj,fj));
        out.energy+=.4*m.stiffness*std::pow(overlap,2.5);
        out.power+=-(fn-elastic)*vn+dot(ft,vt);
        double inverse_mass=(b[i].fixed?0:1/b[i].mass)+(b[j].fixed?0:1/b[j].mass);
        if(inverse_mass>0) {
            double meff=1/inverse_mass;
            double tangent=1.5*m.stiffness*std::sqrt(std::max(dx,overlap));
            if(tangent>0) out.dt=std::min(out.dt,.08*std::sqrt(meff/tangent));
            if(m.damping>0) out.dt=std::min(out.dt,.1*meff/m.damping);
            double inverse_tangent_mass=inverse_mass+
                (b[i].fixed?0:dot(ai,ai)/b[i].inertia())+
                (b[j].fixed?0:dot(aj,aj)/b[j].inertia());
            if(m.friction*fn>0)
                out.dt=std::min(out.dt,.1*m.slip_speed/(m.friction*fn*inverse_tangent_mass));
        }
    }
    // Limit approaching bodies before the first finite-overlap evaluation too.
    if(m.stiffness>0) for(const auto &a:b) if(!a.fixed)
        out.dt=std::min(out.dt,.04*std::sqrt(a.mass/(1.5*m.stiffness*std::sqrt(dx))));
    return out;
}
struct Stage {
    Rates fluid;
    std::vector<Load> hydrodynamic,total,external,mechanical;
    ContactResult contact;
};
using FluidEvaluator=std::function<Rates(const std::vector<Sample>&,const Pairs&,const Parameters&)>;
// Callback is a PURE force/torque evaluation, not a second integrator. It may
// supply neck elasticity or other external DEM laws at each synchronized stage.
// Its result must contain contacts/external loads only, never fluid loads.
using ExtraLoads=std::function<std::vector<Load>(const std::vector<Sphere>&,double)>;
class System {
    pf_neighbors::Verlet graph_;
public:
    System(const System&)=delete;
    System& operator=(const System&)=delete;
    System(System&&)=default;
    System& operator=(System&&)=default;
    Parameters parameters;
    ContactModel contact_model;
    std::vector<Sample> state;
    std::vector<Sphere> bodies;
    std::vector<Site> sites;
    size_t nf;
    double time=0.;long long step=0;
    Vector supportP{},supportL{},referenceP{},referenceL{},appliedP{},appliedL{};
    std::string configuration_signature="-";
    FluidEvaluator evaluator=[](const auto &s,const auto &g,const auto &p){return evaluate(s,g,p);};
    ExtraLoads extra_loads;
    std::shared_ptr<MechanicalLaw> mechanics;
    // Set an additional explicit stability limit when an external law is used.
    bool legacy_overlap_mask=false; // explicit continuation of pre-v3 checkpoints only
    double external_max_dt=std::numeric_limits<double>::infinity();
    // A temporary clamp preserves initial motion for release and accounts for
    // its impulse in the support ledger, independently of permanent fixity.
    struct ReleaseState {bool fixed=false;Vector velocity{},angular_momentum{};};
    double equilibration_end=0.;
    bool equilibration_active=false;
    std::vector<ReleaseState> release_states;
    void begin_equilibration(double duration){
        if(!(duration>0&&std::isfinite(duration))||time!=0||equilibration_end!=0)
            throw std::runtime_error("Equilibration must be configured once at fresh time zero");
        hold_bodies(duration);
    }
    void hold_bodies(double duration){
        const double end=time+duration;
        if(bodies.empty() || equilibration_active || !release_states.empty() ||
           !(duration>0 && std::isfinite(end) && end>time))
            throw std::runtime_error("Cannot start a new fluid hold in this checkpoint state");
        release_states.reserve(bodies.size());
        equilibration_end=end;equilibration_active=true;
        for(auto &b:bodies){release_states.push_back({b.fixed,b.velocity,b.angular_momentum});
            auto momentum=scale(b.velocity,b.mass);
            supportP=add(supportP,momentum);
            supportL=add(supportL,add(cross(b.center,momentum),b.angular_momentum));
            b.fixed=true;b.velocity={};b.angular_momentum={};}
        sync_sites(state,bodies,sites);
    }
    void finish_equilibration(){
        if(!equilibration_active||time<equilibration_end)return;
        for(size_t j=0;j<bodies.size();++j){auto &b=bodies[j];const auto &saved=release_states.at(j);
            b.fixed=saved.fixed;b.velocity=saved.velocity;b.angular_momentum=saved.angular_momentum;
            auto momentum=scale(b.velocity,b.mass);
            supportP=subtract(supportP,momentum);
            supportL=subtract(supportL,add(cross(b.center,momentum),b.angular_momentum));}
        release_states.clear();equilibration_active=false;sync_sites(state,bodies,sites);
    }

    System(Parameters p,std::vector<Sample> s,size_t fluid_count,
           std::vector<Sphere> b,std::vector<Site> mapping,ContactModel cm={},double skin_factor=.25)
      :graph_(fluid_count,p.support(),skin_factor*p.support()),parameters(p),contact_model(cm),
       state(std::move(s)),bodies(std::move(b)),sites(std::move(mapping)),nf(fluid_count) {
        p.validate();cm.validate();
        if(nf==0||!(skin_factor>=0&&std::isfinite(skin_factor)))
            throw std::invalid_argument("Invalid liquid sample count or neighbor skin");
        pf_neighbors::validate_layout(state,nf);
        if(bodies.empty()) {
            if(!sites.empty())throw std::runtime_error("Unexpected rigid sites in calibration");
            for(size_t i=nf;i<state.size();++i)
                if(state[i].v!=Vector{})throw std::runtime_error("Flat wall samples must be stationary");
        } else {
            validate_sites(state,bodies,sites);
        }
        std::vector<size_t> count(bodies.size());
        for(const auto &site:sites){++count[site.body];}
        for(auto n:count) if(n==0)throw std::runtime_error("Unresolved DEM sphere: no boundary samples");
        for(const auto &a:bodies) if(a.fixed && (dot(a.velocity,a.velocity)>0||dot(a.angular_momentum,a.angular_momentum)>0))
            throw std::runtime_error("Fixed DEM bodies must have zero velocity and spin");
        sync_sites(state,bodies,sites);
        referenceP=total_momentum(state,bodies);referenceL=total_momentum(state,bodies,true);
    }
    Stage stage(const std::vector<Sample> &s,const std::vector<Sphere> &b,double t) {
        // Interior sites integrate solid VOLUME, not exposed surface area.
        // Partition unity chi_k(x)/sum_j chi_j(x) counts the union once instead
        // of deleting the intersection from every sphere. Keep samples/mass
        // immutable; derive transient quadrature weights at each force stage.
        std::vector<std::vector<size_t>> overlaps(b.size());
        bool any_overlap=false;
        for(size_t i=0;i<b.size();++i)for(size_t j=i+1;j<b.size();++j){
            auto d=subtract(b[i].center,b[j].center);double reach=b[i].radius+b[j].radius;
            if(dot(d,d)<reach*reach){overlaps[i].push_back(j);overlaps[j].push_back(i);any_overlap=true;}}
        std::vector<Sample> weighted;
        if(any_overlap){weighted=s;
            for(const auto &site:sites){size_t cover=1;
                for(size_t j:overlaps[site.body]){auto d=subtract(s[site.particle].x,b[j].center);
                    if(dot(d,d)<b[j].radius*b[j].radius)++cover;}
                weighted[site.particle].quadrature_weight=1./cover;}}
        const auto &candidates=graph_.graph(s);
        Pairs legacy_pairs;
        if(legacy_overlap_mask&&any_overlap){legacy_pairs.reserve(candidates.size());
            for(auto ij:candidates)if(weighted[ij.first].quadrature_weight==1.&&weighted[ij.second].quadrature_weight==1.)legacy_pairs.push_back(ij);}
        Stage out;out.fluid=evaluator(legacy_overlap_mask?s:(any_overlap?weighted:s),
            legacy_overlap_mask&&any_overlap?legacy_pairs:candidates,parameters);
        out.hydrodynamic=loads_from_fluid(s,out.fluid,b,sites);
        out.contact=contacts(b,contact_model,parameters.dx);out.total=out.hydrodynamic;
        auto external=extra_loads?extra_loads(b,t):std::vector<Load>(b.size());
        out.external=external;
        out.mechanical=mechanics&&!equilibration_active?mechanics->evaluate(b,t):std::vector<Load>(b.size());
        if(out.mechanical.size()!=b.size())throw std::runtime_error("Mechanical load count mismatch");
        if(external.size()!=b.size())throw std::runtime_error("DEM callback load count mismatch");
        for(size_t j=0;j<b.size();++j) {
            if(!finite(external[j].force)||!finite(external[j].torque))throw std::runtime_error("Nonfinite external DEM load");
            out.total[j].force=add(out.total[j].force,add(out.mechanical[j].force,add(out.contact.loads[j].force,external[j].force)));
            out.total[j].torque=add(out.total[j].torque,add(out.mechanical[j].torque,add(out.contact.loads[j].torque,external[j].torque)));
        }
        return out;
    }
    Stage rates(){return stage(state,bodies,time);}
    double timestep(const Stage &r) const {
        double dt=std::min({pf_liquid::timestep(state,r.fluid,parameters),r.contact.dt,external_max_dt});
        if(equilibration_active)return std::min(pf_liquid::timestep(state,r.fluid,parameters),external_max_dt);
        if(mechanics)dt=std::min(dt,mechanics->timestep(bodies));
        for(size_t i=0;i<bodies.size();++i) if(!bodies[i].fixed) {
            const auto &b=bodies[i];double w=std::sqrt(dot(b.omega(),b.omega()));
            double speed=std::sqrt(dot(b.velocity,b.velocity))+w*b.radius;
            dt=std::min(dt,.15*parameters.dx/(parameters.sound+speed));
            double acc=std::sqrt(dot(r.total[i].force,r.total[i].force))/b.mass+
                       std::sqrt(dot(r.total[i].torque,r.total[i].torque))*b.radius/b.inertia();
            if(acc>0)dt=std::min(dt,.15*std::sqrt(parameters.dx/acc));
            if(w>0)dt=std::min(dt,.1/w);
        }
        return dt;
    }
    // Synchronized midpoint evaluation. Trapezoidal center/particle drift gives
    // exact central-force orbital torque balance at the force-stage positions.
    void advance(double dt,const Stage *initial=nullptr) {
        if(!(dt>0&&std::isfinite(dt)))throw std::invalid_argument("Invalid DEM timestep");
        if(equilibration_active&&dt>equilibration_end-time)throw std::runtime_error("Step crosses equilibration release time");
        // Avoid a rounding-sized extra hold step at the release boundary.
        const bool releases_equilibration=equilibration_active&&dt==equilibration_end-time;
        Stage own;if(!initial){own=rates();initial=&own;}
        auto mid=state;auto mb=bodies;
        for(size_t i=0;i<nf;++i) {
            mid[i].x=add(state[i].x,scale(state[i].v,.5*dt));
            mid[i].v=add(state[i].v,scale(initial->fluid.force[i],.5*dt/state[i].mass));
            mid[i].rho+=.5*dt*initial->fluid.density[i];validate(mid[i]);
        }
        for(size_t j=0;j<bodies.size();++j) if(!bodies[j].fixed) {
            const auto &b=bodies[j];auto &m=mb[j];
            m.center=add(b.center,scale(b.velocity,.5*dt));
            m.velocity=add(b.velocity,scale(initial->total[j].force,.5*dt/b.mass));
            m.angular_momentum=add(b.angular_momentum,scale(initial->total[j].torque,.5*dt));
            m.orientation=normalize(add_quaternion(b.orientation,scale_quaternion(derivative(b.orientation,b.omega()),.5*dt)));
        }
        if(mechanics&&!equilibration_active)mechanics->begin_step(bodies,time,dt);
        struct TrialGuard {
            std::shared_ptr<MechanicalLaw> law;
            ~TrialGuard(){if(law)law->rollback();}
        } guard{equilibration_active?nullptr:mechanics};
        sync_sites(mid,mb,sites);auto r=stage(mid,mb,time+.5*dt);
        auto next=state;auto nb=bodies;
        Vector dp{},dl{},ep{},el{};
        // A flat wall has no DEM body. Its fixed quadrature sites exchange
        // momentum with liquid and deliver the reaction to the support.
        if(bodies.empty()) {
            for(size_t i=nf;i<state.size();++i) {
                const auto impulse=scale(r.fluid.force[i],dt);
                dp=add(dp,impulse);
                dl=add(dl,cross(state[i].x,impulse));
            }
        }
        for(size_t i=0;i<nf;++i) {
            next[i].v=add(state[i].v,scale(r.fluid.force[i],dt/state[i].mass));
            next[i].x=add(state[i].x,scale(add(state[i].v,next[i].v),.5*dt));
            next[i].rho+=dt*r.fluid.density[i];validate(next[i]);
        }
        for(size_t j=0;j<bodies.size();++j) {
            const auto &b=bodies[j];auto &n=nb[j];
            ep=add(ep,scale(r.external[j].force,dt));
            el=add(el,scale(add(cross(mb[j].center,r.external[j].force),r.external[j].torque),dt));
            if(b.fixed) {
                // Impulse delivered to the support is included in diagnostics.
                dp=add(dp,scale(r.total[j].force,dt));
                dl=add(dl,scale(add(cross(mb[j].center,r.total[j].force),r.total[j].torque),dt));
            } else {
                n.velocity=add(b.velocity,scale(r.total[j].force,dt/b.mass));
                n.center=add(b.center,scale(add(b.velocity,n.velocity),.5*dt));
                n.angular_momentum=add(b.angular_momentum,scale(r.total[j].torque,dt));
                n.orientation=normalize(add_quaternion(b.orientation,scale_quaternion(derivative(mb[j].orientation,mb[j].omega()),dt)));
                n.validate();
            }
        }
        sync_sites(next,nb,sites);
        if(mechanics&&!equilibration_active)mechanics->commit(mb,nb,time+.5*dt,dt);
        state.swap(next);bodies.swap(nb);
        supportP=add(supportP,dp);supportL=add(supportL,dl);
        appliedP=add(appliedP,ep);appliedL=add(appliedL,el);time=releases_equilibration?equilibration_end:time+dt;++step;finish_equilibration();
    }
    Vector momentum_error(bool angular=false) const {
        return subtract(subtract(add(total_momentum(state,bodies,angular),angular?supportL:supportP),angular?referenceL:referenceP),angular?appliedL:appliedP);
    }
    // Insert ahead of the solid suffix and retain all body-local site mappings.
    void append_liquid(const std::vector<Sample> &added){
        Vector dp{},dl{};
        for(const auto &a:added){
            if(a.solid)throw std::invalid_argument("Expected liquid addition");
            validate(a);
            const auto impulse=scale(a.v,a.mass);
            dp=add(dp,impulse);
            dl=add(dl,cross(a.x,impulse));
        }
        state.insert(state.begin()+static_cast<std::ptrdiff_t>(nf),added.begin(),added.end());
        for(auto &site:sites)site.particle+=added.size();
        nf+=added.size();
        graph_.set_fluid_count(nf);
        appliedP=add(appliedP,dp);
        appliedL=add(appliedL,dl);
    }
    size_t neighbor_builds() const{return graph_.builds();}
    size_t neighbor_reuses() const{return graph_.reuses();}
    // Full state includes immutable sample ownership and local reference sites.
    // No re-generation from rounded positions on resume; cache is rebuilt.
    void save(const std::string &path) const {
        if(extra_loads)throw std::runtime_error("External DEM law must implement its own checkpoint before saving");
        std::string tmp=path+".tmp";std::ofstream f(tmp);
        if(!f)throw std::runtime_error("Cannot write DEM checkpoint");
        f<<std::setprecision(17)<<"PFDEM "<<(bodies.empty()?6:5)<<'\n'<<std::quoted(configuration_signature)<<'\n';
        const auto &p=parameters;const auto &m=contact_model;
        f<<p.dx<<' '<<p.rho0<<' '<<p.sound<<' '<<p.mu<<' '<<p.sff<<' '<<p.sfs<<'\n';
        f<<m.stiffness<<' '<<m.damping<<' '<<m.friction<<' '<<m.slip_speed<<'\n';
        f<<nf<<' '<<state.size()<<' '<<bodies.size()<<' '<<sites.size()<<' '<<time<<' '<<step<<' '<<int(legacy_overlap_mask)<<'\n';
        for(auto v:{supportP,supportL,referenceP,referenceL,appliedP,appliedL}){for(double x:v)f<<x<<' ';f<<'\n';}
        for(const auto &b:bodies){for(auto v:{b.center,b.velocity,b.angular_momentum})for(double x:v)f<<x<<' ';
            for(double x:b.orientation)f<<x<<' ';
            f<<b.radius<<' '<<b.mass<<' '<<b.fixed<<'\n';}
        for(const auto &a:state){for(auto v:{a.x,a.v})for(double x:v)f<<x<<' ';f<<a.rho<<' '<<a.mass<<' '<<a.solid<<'\n';}
        for(const auto &s:sites)f<<s.particle<<' '<<s.body<<' '<<s.local[0]<<' '<<s.local[1]<<' '<<s.local[2]<<'\n';
        f<<"EQUILIBRATION "<<equilibration_end<<' '<<int(equilibration_active)<<' '<<release_states.size()<<'\n';
        for(const auto &r:release_states){f<<int(r.fixed);for(auto v:{r.velocity,r.angular_momentum})for(double x:v)f<<' '<<x;f<<'\n';}
        f<<std::quoted(mechanics?mechanics->type():"none")<<'\n';
        if(mechanics)mechanics->save(f);
        f.close();if(!f||std::rename(tmp.c_str(),path.c_str())!=0)throw std::runtime_error("Cannot commit DEM checkpoint");
    }
    static System load(const std::string &path,MechanicalFactory factory={},bool allow_legacy_overlap=false) {
        std::ifstream f(path);std::string magic;int version;Parameters p;ContactModel cm;
        size_t nf,n,nb,ns;double t;long long step;
        if(!(f>>magic>>version)||magic!="PFDEM"||(version<1||version>6))throw std::runtime_error("Invalid DEM checkpoint");
        std::string signature;f>>std::quoted(signature);
        f>>p.dx>>p.rho0>>p.sound>>p.mu>>p.sff>>p.sfs;
        f>>cm.stiffness>>cm.damping>>cm.friction>>cm.slip_speed;
        if(!(f>>nf>>n>>nb>>ns>>t>>step)||n>10000000||nb>100000||ns>n||nf>n||
           (version==6?(nb!=0||ns!=0):(nb==0||n-nf!=ns))||
           !(t>=0&&std::isfinite(t))||step<0)throw std::runtime_error("Invalid DEM checkpoint header");
        bool stored_legacy=false;
        if(version>=4){int rule=-1;f>>rule;if(rule!=0&&rule!=1)throw std::runtime_error("Invalid overlap checkpoint rule");stored_legacy=rule;}
        else if(version<3)stored_legacy=allow_legacy_overlap;
        Vector ledger[6];for(auto &v:ledger)for(double &x:v)f>>x;
        std::vector<Sphere> b(nb);std::vector<Sample> s(n);std::vector<Site> sites(ns);
        for(auto &a:b){for(auto *v:{&a.center,&a.velocity,&a.angular_momentum})for(double &x:*v)f>>x;
            for(double &x:a.orientation)f>>x;
            int fixed=-1;f>>a.radius>>a.mass>>fixed;
            if(fixed!=0&&fixed!=1)throw std::runtime_error("Invalid checkpoint fixity");
            a.fixed=fixed;}
        for(auto &a:s){for(auto *v:{&a.x,&a.v})for(double &x:*v)f>>x;
            int solid=-1;f>>a.rho>>a.mass>>solid;if(solid!=0&&solid!=1)throw std::runtime_error("Invalid checkpoint phase");a.solid=solid;validate(a);}
        for(auto &a:sites)f>>a.particle>>a.body>>a.local[0]>>a.local[1]>>a.local[2];
        if(!f)throw std::runtime_error("Truncated DEM checkpoint");
        double eq_end=0.;bool eq_active=false;std::vector<ReleaseState> releases;
        if(version>=5){std::string tag;int active=-1;size_t count=0;
            f>>tag>>eq_end>>active>>count;
            if(!f||tag!="EQUILIBRATION"||!(eq_end>=0&&std::isfinite(eq_end))||(active!=0&&active!=1)||
               (active?(count!=nb||!(eq_end>t)):count!=0))throw std::runtime_error("Invalid equilibration checkpoint metadata");
            eq_active=active;releases.resize(count);
            for(auto &r:releases){int fixed=-1;f>>fixed;r.fixed=fixed;
                for(auto *v:{&r.velocity,&r.angular_momentum})for(double &x:*v)f>>x;
                if(!f||(fixed!=0&&fixed!=1)||!finite(r.velocity)||!finite(r.angular_momentum)||
                   (r.fixed&&(dot(r.velocity,r.velocity)>0||dot(r.angular_momentum,r.angular_momentum)>0)))
                    throw std::runtime_error("Invalid equilibration release state");}
            if(eq_active)for(const auto &a:b)if(!a.fixed)throw std::runtime_error("Equilibrating checkpoint contains moving body");
        }
        std::shared_ptr<MechanicalLaw> law;
        if(version>=2){
            std::string kind;if(!(f>>std::quoted(kind)))throw std::runtime_error("Missing mechanical checkpoint tag");
            if(kind!="none"){
                if(!factory)throw std::runtime_error("Checkpoint requires a soot-dem-enabled executable");
                law=factory(kind,f,b);
            }
        }
        std::string extra;
        if(f>>extra)throw std::runtime_error("Trailing DEM checkpoint data");
        if(version<3&&!allow_legacy_overlap)for(size_t i=0;i<b.size();++i)for(size_t j=i+1;j<b.size();++j){
            auto d=subtract(b[i].center,b[j].center);double reach=b[i].radius+b[j].radius;
            if(dot(d,d)<reach*reach)throw std::runtime_error("Legacy checkpoint has overlapping spheres with obsolete volume masking; use --legacy-overlap-mask to continue the old rule, or start fresh");}

        for(const auto &v:ledger)if(!finite(v))throw std::runtime_error("Invalid checkpoint ledger");
        // Constructor regenerates quadrature from exact local sites. Check that
        // the persisted world state agrees before accepting the checkpoint.
        auto saved=s;System result(p,std::move(s),nf,std::move(b),std::move(sites),cm);
        for(size_t i=nf;i<saved.size();++i)for(int k=0;k<3;++k)
            if(std::abs(saved[i].x[k]-result.state[i].x[k])>1e-12 || std::abs(saved[i].v[k]-result.state[i].v[k])>1e-12)
                throw std::runtime_error("Checkpoint site/body mismatch");
        result.supportP=ledger[0];result.supportL=ledger[1];result.referenceP=ledger[2];result.referenceL=ledger[3];
        result.appliedP=ledger[4];result.appliedL=ledger[5];result.configuration_signature=signature;
        result.equilibration_end=eq_end;result.equilibration_active=eq_active;result.release_states=std::move(releases);
        result.legacy_overlap_mask=stored_legacy;
        result.time=t;result.step=step;result.mechanics=std::move(law);return result;
    }
};
} // namespace pf_dem
