#pragma once
// Deterministic per-shard reciprocal reduction. Only the SPHinXsys adapter
// depends on oneTBB; the pair/continuity physics remains dependency-free.
#include "pf_liquid.h"
#include <tbb/parallel_for.h>
#include <algorithm>
#include <vector>
namespace pf_parallel {
inline pf_liquid::Rates evaluate(const std::vector<pf_liquid::Sample> &s,
                                 const pf_liquid::Pairs &pairs,
                                 const pf_liquid::Parameters &p,unsigned threads) {
    if(threads<=1 || pairs.size()<4096) return pf_liquid::evaluate(s,pairs,p);
    auto P=pf_liquid::pressures(s,p);
    const size_t shards=std::min<size_t>(threads,(pairs.size()+4095)/4096);
    std::vector<pf_liquid::Rates> partial(shards);
    tbb::parallel_for(size_t(0),shards,[&](size_t k) {
        const size_t first=(pairs.size()*k)/shards;
        const size_t last=(pairs.size()*(k+1))/shards;
        partial[k]=pf_liquid::evaluate_range(s,pairs,p,P,first,last);
    });
    auto result=pf_liquid::empty_rates(s.size());
    for(const auto &part:partial) pf_liquid::add_rates(result,part);
    return result;
}
} // namespace pf_parallel
