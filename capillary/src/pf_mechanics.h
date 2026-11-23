#pragma once
#include "pf_rigid.h"
#include <memory>
#include <string>
#include <functional>
#include <istream>
#include <ostream>
namespace pf_dem {
// A history-dependent mechanical backend. Force queries must not change its
// committed state. A failed step discards trial history. Only commit advances it.
struct MechanicalLaw {
    virtual ~MechanicalLaw()=default;
    virtual std::vector<pf_rigid::Load> evaluate(const std::vector<pf_rigid::Sphere>&,double) const=0;
    virtual void begin_step(const std::vector<pf_rigid::Sphere>&,double,double)=0;
    virtual void commit(const std::vector<pf_rigid::Sphere>&,const std::vector<pf_rigid::Sphere>&,double,double)=0;
    virtual void rollback() noexcept=0;
    virtual double timestep(const std::vector<pf_rigid::Sphere>&) const=0;
    virtual std::string type() const=0;
    virtual void save(std::ostream&) const=0;
    virtual size_t neck_count() const{return 0;}
    virtual std::vector<std::pair<size_t,size_t>> necks() const{return {};}
};
using MechanicalFactory=std::function<std::shared_ptr<MechanicalLaw>(const std::string&,std::istream&,const std::vector<pf_rigid::Sphere>&)>;
}
