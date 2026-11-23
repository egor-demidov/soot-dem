#pragma once
// Dependency-free, deterministic spatial hashing and Verlet candidate pairs.
// The liquid solver evaluates the current distance for every candidate and
// ignores r >= support, so a skin changes *when* the search runs, not physics.
// Particles are ordered [fluid ... fluid, solid ... solid].
#include "pf_liquid.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace pf_neighbors {
namespace pf=pf_liquid;
using Index=std::int64_t;
struct Cell {
    Index x=0,y=0,z=0;
    bool operator==(const Cell &o) const noexcept {return x==o.x&&y==o.y&&z==o.z;}
};
struct CellHash {
    size_t operator()(const Cell &c) const noexcept {
        // Hashing affects only lookup speed; candidate ordering is sorted later.
        auto mix=[](std::uint64_t x) {
            x^=x>>30; x*=UINT64_C(0xbf58476d1ce4e5b9);
            x^=x>>27; x*=UINT64_C(0x94d049bb133111eb);
            return x^(x>>31);
        };
        std::uint64_t h=mix(static_cast<std::uint64_t>(c.x));
        h^=mix(static_cast<std::uint64_t>(c.y)+UINT64_C(0x9e3779b97f4a7c15));
        h^=mix(static_cast<std::uint64_t>(c.z)+UINT64_C(0x3c6ef372fe94f82a));
        return static_cast<size_t>(h);
    }
};
inline Cell cell_of(const pf::Vector &position,double width) {
    Cell c;
    Index *fields[3]={&c.x,&c.y,&c.z};
    for(int k=0;k<3;++k) {
        const double q=std::floor(position[k]/width);
        // Retain a safety margin for +/-1 neighboring-cell arithmetic.
        if(!std::isfinite(q) || q<static_cast<double>(std::numeric_limits<Index>::min()/2) ||
           q>static_cast<double>(std::numeric_limits<Index>::max()/2))
            throw std::runtime_error("Neighbor cell index outside representable range");
        *fields[k]=static_cast<Index>(q);
    }
    return c;
}
inline void validate_layout(const std::vector<pf::Sample> &s,size_t fluid_count) {
    if(fluid_count>s.size()) throw std::invalid_argument("Neighbor fluid count exceeds samples");
    for(size_t i=0;i<s.size();++i) {
        if(!pf::finite(s[i].x)) throw std::invalid_argument("Nonfinite neighbor position");
        if(s[i].solid!=(i>=fluid_count))
            throw std::invalid_argument("Neighbor samples must be ordered liquid then solid");
    }
}
inline pf::Pairs build_cells(const std::vector<pf::Sample> &s,size_t fluid_count,double cutoff) {
    if(!(cutoff>0 && std::isfinite(cutoff)))
        throw std::invalid_argument("Invalid neighbor search cutoff");
    validate_layout(s,fluid_count);
    using Grid=std::unordered_map<Cell,std::vector<size_t>,CellHash>;
    Grid grid;
    grid.reserve(s.size()*2+1);
    std::vector<Cell> keys;
    keys.reserve(s.size());
    for(size_t i=0;i<s.size();++i) {
        auto c=cell_of(s[i].x,cutoff);
        grid[c].push_back(i);
        keys.push_back(c);
    }
    const double cutoff2=cutoff*cutoff;
    pf::Pairs pairs;
    for(size_t i=0;i<fluid_count;++i) {
        const auto c=keys[i];
        for(int ix=-1;ix<=1;++ix) for(int iy=-1;iy<=1;++iy) for(int iz=-1;iz<=1;++iz) {
            auto it=grid.find({c.x+ix,c.y+iy,c.z+iz});
            if(it==grid.end()) continue;
            for(auto j:it->second) {
                if(j<=i) continue;
                const auto d=pf::subtract(s[i].x,s[j].x);
                if(pf::dot(d,d)<cutoff2) pairs.emplace_back(i,j);
            }
        }
    }
    // Stable order makes serial and deterministic-shard rate summation reproducible.
    std::sort(pairs.begin(),pairs.end());
    return pairs;
}
inline pf::Pairs brute_force(const std::vector<pf::Sample> &s,size_t fluid_count,double cutoff) {
    if(!(cutoff>0 && std::isfinite(cutoff))) throw std::invalid_argument("Invalid brute-force cutoff");
    validate_layout(s,fluid_count);
    pf::Pairs pairs;
    const double cutoff2=cutoff*cutoff;
    for(size_t i=0;i<fluid_count;++i) for(size_t j=i+1;j<s.size();++j) {
        auto d=pf::subtract(s[i].x,s[j].x);
        if(pf::dot(d,d)<cutoff2) pairs.emplace_back(i,j);
    }
    return pairs;
}
class Verlet {
    size_t fluid_count_;
    double cutoff_,skin_,candidate_cutoff_;
    std::vector<pf::Vector> reference_positions_;
    pf::Pairs candidates_;
    size_t builds_=0,reuses_=0;
    bool force_rebuild_=true;
public:
    Verlet(size_t fluid_count,double cutoff,double skin)
      :fluid_count_(fluid_count),cutoff_(cutoff),skin_(skin),candidate_cutoff_(cutoff+skin) {
        if(!(cutoff>0&&std::isfinite(cutoff))||!(skin>=0&&std::isfinite(skin))||
           !std::isfinite(candidate_cutoff_)) throw std::invalid_argument("Invalid Verlet cutoff/skin");
    }
    const pf::Pairs& graph(const std::vector<pf::Sample> &state) {
        if(state.size()<fluid_count_) throw std::invalid_argument("Incorrect Verlet particle count");
        bool rebuild=force_rebuild_||reference_positions_.size()!=state.size()||skin_==0.;
        if(!rebuild) {
            const double threshold2=skin_*skin_/4.;
            for(size_t i=0;i<state.size();++i) {
                const auto d=pf::subtract(state[i].x,reference_positions_[i]);
                const double dsq=pf::dot(d,d);
                if(!std::isfinite(dsq)) throw std::runtime_error("Nonfinite Verlet displacement");
                // Conservative: relative displacement for a pair <= 2 max individual displacement.
                if(dsq>=threshold2) {rebuild=true;break;}
            }
        }
        if(rebuild) {
            auto next=build_cells(state,fluid_count_,candidate_cutoff_);
            candidates_.swap(next);
            reference_positions_.resize(state.size());
            for(size_t i=0;i<state.size();++i) reference_positions_[i]=state[i].x;
            force_rebuild_=false;
            ++builds_;
        } else ++reuses_;
        return candidates_;
    }
    void invalidate() noexcept {force_rebuild_=true;}
    void set_fluid_count(size_t count) noexcept {
        fluid_count_=count;
        invalidate();
    }
    const pf::Pairs& candidates() const noexcept {return candidates_;}
    size_t builds() const noexcept {return builds_;}
    size_t reuses() const noexcept {return reuses_;}
    double cutoff() const noexcept {return cutoff_;}
    double skin() const noexcept {return skin_;}
};
} // namespace pf_neighbors
