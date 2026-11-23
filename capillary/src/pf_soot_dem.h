#pragma once
// Adapter to actual soot-dem aggregate / libgran force functors. This header is
// compiled only for pf_soot_dem with a build-local history-access overlay.
#include "pf_mechanics.h"
#include "pf_soot_settings.h"
#include <Eigen/Eigen>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <memory>
#include <limits>
#include <iomanip>
#include <pf_soot_build_id.h>
#include <aggregate.h>
#include <break_neck.h>
namespace pf_dem {
class SootDemLaw final:public MechanicalLaw {
public:
    using Vec=Eigen::Vector3d;
    using Model=aggregate<Vec,double>;
    using Spring=std::tuple<Vec,Vec,Vec>;
private:
    SootSettings settings_;
    std::unique_ptr<Model> committed_,trial_;
    double radius_,mass_,inertia_;
    size_t count_;
    static Vec eigen(pf_rigid::Vector v){return {v[0],v[1],v[2]};}
    static pf_rigid::Vector array(const Vec &v){return {v[0],v[1],v[2]};}
    static bool equal(double a,double b){return std::abs(a-b)<=1e-12*std::max(std::abs(a),std::abs(b));}
    std::unique_ptr<Model> construct(const std::vector<pf_rigid::Sphere>&b,bool initialize_graph) const {
        const auto &v=settings_.v;std::vector<Vec>x;for(const auto &a:b)x.push_back(eigen(a.center));
        return std::make_unique<Model>(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],
               v[14],v[15],v[16],v[17],v[18],v[19],v[20],v[21],initialize_graph?v[22]:0.,v[23],v[24],
               x,count_,radius_,mass_,inertia_,0.,Vec::Zero(),0.);
    }
    // Native functors return accelerations, not forces. Each canonical pair is
    // evaluated once, then assembled with exact reciprocal force and moment.
    std::vector<pf_rigid::Load> apply(Model &m,const std::vector<pf_rigid::Sphere>&b,double t) const {
        if(b.size()!=count_)throw std::runtime_error("soot-dem body count changed");
        std::vector<Vec>x,v,w,theta(count_,Vec::Zero());
        for(const auto &a:b){x.push_back(eigen(a.center));v.push_back(a.fixed?Vec::Zero().eval():eigen(a.velocity));w.push_back(a.fixed?Vec::Zero().eval():eigen(a.omega()));}
        std::vector<pf_rigid::Load> loads(count_);
        for(size_t i=0;i<count_;++i)for(size_t j=i+1;j<count_;++j){
            if((x[i]-x[j]).squaredNorm()==0)throw std::runtime_error("Coincident soot-dem particle centers");
            auto response=m(i,j,x,v,theta,w,t);
            Vec fi=response.first*mass_,ti=response.second*inertia_;
            if(!fi.allFinite()||!ti.allFinite())throw std::runtime_error("Nonfinite soot-dem force/torque");
            auto force=array(fi),torque=array(ti);
            loads[i].force=pf_liquid::add(loads[i].force,force);
            loads[i].torque=pf_liquid::add(loads[i].torque,torque);
            loads[j].force=pf_liquid::subtract(loads[j].force,force);
            auto tj=pf_liquid::scale(pf_liquid::add(torque,pf_liquid::cross(pf_liquid::subtract(b[i].center,b[j].center),force)),-1.);
            loads[j].torque=pf_liquid::add(loads[j].torque,tj);
        }
        return loads;
    }
    static void write_springs(std::ostream &f,const std::vector<Spring>&springs){
        f<<springs.size()<<'\n';
        for(const auto &s:springs)for(const auto &v:{std::get<0>(s),std::get<1>(s),std::get<2>(s)})
            f<<v[0]<<' '<<v[1]<<' '<<v[2]<<'\n';
    }
    void read_springs(std::istream &f,std::vector<Spring>&springs){
        size_t n;if(!(f>>n)||n!=count_*count_)throw std::runtime_error("Invalid soot-dem spring count");
        for(auto &s:springs)for(auto *v:{&std::get<0>(s),&std::get<1>(s),&std::get<2>(s)}){
            f>>(*v)[0]>>(*v)[1]>>(*v)[2];if(!f||!v->allFinite())throw std::runtime_error("Invalid soot-dem spring state");}
    }
public:
    SootDemLaw(SootSettings settings,const std::vector<pf_rigid::Sphere>&b,bool initialize=true)
      :settings_(std::move(settings)),radius_(b.at(0).radius),mass_(b[0].mass),inertia_(b[0].inertia()),count_(b.size()) {
        settings_.validate(count_);
        for(const auto &a:b)if(!equal(a.radius,radius_)||!equal(a.mass,mass_))
            throw std::runtime_error("soot-dem aggregate model requires equal primary radius and mass");
        committed_=construct(b,initialize&&settings_.neck_mode==1);
        if(initialize&&settings_.neck_mode!=1){auto &mask=committed_->pf_sinter_model().bonded_contacts;
            std::fill(mask.begin(),mask.end(),false);
            for(auto [i,j]:settings_.edges)mask[i*count_+j]=mask[j*count_+i]=true;}
    }
    std::vector<pf_rigid::Load> evaluate(const std::vector<pf_rigid::Sphere>&b,double t) const override {
        Model temporary(trial_?*trial_:*committed_);temporary.pf_sinter_model().pf_set_dt(0.);
        return apply(temporary,b,t);
    }
    void begin_step(const std::vector<pf_rigid::Sphere>&b,double t,double dt) override {
        auto next=std::make_unique<Model>(*committed_);next->pf_sinter_model().pf_set_dt(.5*dt);
        (void)apply(*next,b,t);trial_=std::move(next);
    }
    void commit(const std::vector<pf_rigid::Sphere>&mid,const std::vector<pf_rigid::Sphere>&end,double t,double dt) override {
        auto next=std::make_unique<Model>(*committed_);next->pf_sinter_model().pf_set_dt(dt);
        (void)apply(*next,mid,t);
        if(settings_.v[26]>0||settings_.v[27]>0){
            std::vector<Vec>x;for(const auto &a:end)x.push_back(eigen(a.center));
            size_t n=0;for(bool active:next->pf_sinter_model().bonded_contacts)if(active)++n;
            std::vector<double> normal(n/2,settings_.v[26]>0?settings_.v[26]:std::numeric_limits<double>::max());
            std::vector<double> shear(n/2,settings_.v[27]>0?settings_.v[27]:std::numeric_limits<double>::max());
            break_strained_necks(*next,x,settings_.v[14],settings_.v[16],settings_.v[18],settings_.v[20],shear,normal,radius_);
            // A broken neck becomes a fresh nonbonded contact; no hidden spring carry-over.
            auto &sn=next->pf_sinter_model();const auto &old=committed_->pf_sinter_model().bonded_contacts;
            for(size_t i=0;i<count_;++i)for(size_t j=i+1;j<count_;++j)if(old[i*count_+j]&&!sn.bonded_contacts[i*count_+j]){
                sn.pf_springs()[i*count_+j]=sn.pf_springs()[j*count_+i]={Vec::Zero(),Vec::Zero(),Vec::Zero()};
                sn.pf_nonbonded().pf_springs()[i*count_+j]=sn.pf_nonbonded().pf_springs()[j*count_+i]={Vec::Zero(),Vec::Zero(),Vec::Zero()};
            }
        }
        next->pf_sinter_model().pf_set_dt(0.);committed_.swap(next);trial_.reset();
    }
    void rollback() noexcept override{trial_.reset();}
    double timestep(const std::vector<pf_rigid::Sphere>&) const override {
        const auto &v=settings_.v;double k=std::max(v[0],v[14]),g=std::max(v[1],v[15]);
        double kt=std::max(v[2],v[16])+std::max(v[6],v[18])+std::max(v[10],v[20]);
        double gt=std::max(v[3],v[17])+std::max(v[7],v[19])+std::max(v[11],v[21]);
        double degree=std::max<size_t>(1,count_-1),dt=v[25];
        double rate=degree*(k/mass_+kt*radius_*radius_/inertia_);
        double drag=degree*(g/mass_+gt*radius_*radius_/inertia_);
        if(rate>0)dt=std::min(dt,.05/std::sqrt(rate));
        if(drag>0)dt=std::min(dt,.05/drag);
        return dt;
    }
    std::string type() const override{return "soot-dem-v1";}
    size_t neck_count() const override{return necks().size();}
    std::vector<std::pair<size_t,size_t>> necks() const override {
        std::vector<std::pair<size_t,size_t>> out;const auto &mask=committed_->pf_sinter_model().bonded_contacts;
        for(size_t i=0;i<count_;++i)for(size_t j=i+1;j<count_;++j)if(mask[i*count_+j])out.emplace_back(i,j);
        return out;
    }
    const Model &model() const{return *committed_;}
    void write_neck_forces(std::ostream &out,const std::vector<pf_rigid::Sphere> &b,
                           long long step,double time) const {
        if(b.size()!=count_)throw std::runtime_error("soot-dem body count changed");
        const auto &sinter=committed_->pf_sinter_model();
        for(const auto &[i,j]:necks()){
            const auto &[sliding,rolling,torsional]=sinter.pf_springs()[i*count_+j];
            const double extension=eigen(pf_liquid::subtract(b[i].center,b[j].center)).norm()-2*radius_;
            out<<step<<','<<time<<','<<i<<','<<j
               <<','<<settings_.v[14]*extension
               <<','<<settings_.v[16]*sliding.norm()
               <<','<<settings_.v[18]*rolling.norm()
               <<','<<settings_.v[20]*torsional.norm()<<'\n';
        }
    }
    void save(std::ostream &f) const override {
        if(trial_)throw std::runtime_error("Cannot checkpoint an uncommitted soot-dem trial");
        f<<std::quoted(PF_SOOT_BACKEND_FINGERPRINT)<<'\n';
        settings_.write(f);const auto &sn=committed_->pf_sinter_model();
        f<<count_*count_<<'\n';for(bool value:sn.bonded_contacts)f<<value<<' ';f<<'\n';
        write_springs(f,sn.pf_springs());write_springs(f,sn.pf_nonbonded().pf_springs());
    }
    static std::shared_ptr<MechanicalLaw> restore(const std::string &kind,std::istream &f,const std::vector<pf_rigid::Sphere>&b){
        if(kind!="soot-dem-v1")throw std::runtime_error("Unsupported mechanical checkpoint backend");
        std::string fingerprint;if(!(f>>std::quoted(fingerprint))||fingerprint!=PF_SOOT_BACKEND_FINGERPRINT)throw std::runtime_error("soot-dem backend code differs from checkpoint; use the same model revision");
        SootSettings settings;settings.read(f,b.size());auto out=std::make_shared<SootDemLaw>(settings,b,false);
        size_t n;if(!(f>>n)||n!=b.size()*b.size())throw std::runtime_error("Invalid soot-dem graph size");
        auto &sn=out->committed_->pf_sinter_model();
        for(size_t i=0;i<n;++i){int active=-1;if(!(f>>active)||(active!=0&&active!=1))throw std::runtime_error("Invalid soot-dem graph state");sn.bonded_contacts[i]=active;}
        for(size_t i=0;i<b.size();++i)for(size_t j=0;j<b.size();++j)
            if((i==j&&sn.bonded_contacts[i*b.size()+j])||sn.bonded_contacts[i*b.size()+j]!=sn.bonded_contacts[j*b.size()+i])throw std::runtime_error("Invalid soot-dem graph symmetry");
        out->read_springs(f,sn.pf_springs());out->read_springs(f,sn.pf_nonbonded().pf_springs());return out;
    }
};
}
