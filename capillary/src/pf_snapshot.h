#pragma once
// Recovery from pre-checkpoint legacy liquid CSVs. These snapshots carry the
// fluid state, but NOT the historical wall impulse integral; the caller must
// explicitly rebase momentum diagnostics at the recovery point.
#include "pf_liquid.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
namespace pf_snapshot {
struct Metadata {long long step=-1; double time=-1;};
inline Metadata find_metadata(std::istream &f,const std::string &filename) {
    std::string header,line;
    if(!std::getline(f,header)||header!="step,time,file")
        throw std::runtime_error("Invalid PF snapshot index header");
    Metadata out;
    bool found=false;
    while(std::getline(f,line)) {
        if(line.empty()) continue;
        std::istringstream r(line);
        std::string step,time,file,extra;
        if(!std::getline(r,step,',')||!std::getline(r,time,',')||
           !std::getline(r,file,',')||std::getline(r,extra,','))
            throw std::runtime_error("Invalid PF snapshot index row");
        if(file==filename) {
            if(found) throw std::runtime_error("Duplicate PF snapshot index entry");
            size_t used=0;
            out.step=std::stoll(step,&used); if(used!=step.size()) throw std::runtime_error("Invalid PF snapshot step");
            out.time=std::stod(time,&used); if(used!=time.size()) throw std::runtime_error("Invalid PF snapshot time");
            found=true;
        }
    }
    if(!found||out.step<0||!(out.time>=0&&std::isfinite(out.time)))
        throw std::runtime_error("PF snapshot file not found in index or invalid time");
    return out;
}
inline std::vector<pf_liquid::Sample> load(std::istream &f,
                                         std::vector<pf_liquid::Sample> state,
                                         size_t nf) {
    std::string header,line;
    if(!std::getline(f,header)||header!="x,y,z,rho,vx,vy,vz")
        throw std::runtime_error("Invalid PF liquid CSV header");
    for(size_t i=0;i<nf;++i) {
        if(!std::getline(f,line)) throw std::runtime_error("Truncated PF liquid CSV");
        std::istringstream row(line);
        double a[7];
        for(int k=0;k<7;++k) {
            std::string item;
            if(!std::getline(row,item,',')) throw std::runtime_error("Truncated PF liquid CSV row");
            size_t used=0;
            a[k]=std::stod(item,&used);
            if(used!=item.size()) throw std::runtime_error("Invalid PF liquid CSV number");
        }
        std::string extra;
        if(std::getline(row,extra,',')) throw std::runtime_error("Extra PF liquid CSV column");
        state[i].x={a[0],a[1],a[2]};
        state[i].rho=a[3]; state[i].v={a[4],a[5],a[6]};
        pf_liquid::validate(state[i]);
        if(state[i].solid) throw std::runtime_error("Solid sample inside PF liquid CSV rows");
    }
    while(std::getline(f,line)) if(!line.empty())
        throw std::runtime_error("Extra PF liquid CSV rows: geometry/particle count changed");
    return state;
}
inline std::pair<Metadata,std::vector<pf_liquid::Sample>> load_file(
       const std::string &filename,const std::vector<pf_liquid::Sample> &initial,size_t nf) {
    namespace fs=std::filesystem;
    const fs::path path(filename);
    std::ifstream index(path.parent_path()/"pf_snapshots.csv");
    if(!index) throw std::runtime_error("Missing pf_snapshots.csv beside recovered snapshot");
    auto meta=find_metadata(index,path.filename().string());
    std::ifstream data(path);
    if(!data) throw std::runtime_error("Cannot open PF liquid CSV snapshot: "+filename);
    return {meta,load(data,initial,nf)};
}
inline pf_liquid::Vector fluid_momentum(const std::vector<pf_liquid::Sample> &s,
                                       size_t nf,bool angular) {
    pf_liquid::Vector out{};
    for(size_t i=0;i<nf;++i) {
        const auto momentum=pf_liquid::scale(s[i].v,s[i].mass);
        out=pf_liquid::add(out,angular?pf_liquid::cross(s[i].x,momentum):momentum);
    }
    return out;
}
} // namespace pf_snapshot
