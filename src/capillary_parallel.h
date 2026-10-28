#pragma once
#include <algorithm>
#include <cstddef>
#include <exception>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// Independent work only. Each worker owns an output slot; callers combine
// results in index order after the join. No exception crosses an OMP boundary.
template<class Function>
void capillaryParallelFor(std::size_t count, Function&& function) {
    int threads=1;
#ifdef _OPENMP
    if (!omp_in_parallel()) threads=int(std::min<std::size_t>(count, std::size_t(omp_get_max_threads())));
#endif
    if (threads<=1) {
        for(std::size_t i=0;i<count;++i) function(i);
        return;
    }
    std::vector<std::exception_ptr> failures(count);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
#endif
    for(std::ptrdiff_t i=0;i<std::ptrdiff_t(count);++i) {
        try { function(std::size_t(i)); }
        catch (...) { failures[std::size_t(i)]=std::current_exception(); }
    }
    for(const auto& failure:failures)if(failure)std::rethrow_exception(failure);
}
