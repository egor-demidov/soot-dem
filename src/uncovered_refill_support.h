#pragma once
#include <array>
#include <cstddef>
#include <utility>
#include <set>
#include <stdexcept>
#include <vector>

namespace capillary_coupling {
using RefillCell = std::array<int,3>;
enum class RefillCellKind { Blocked, NewlyUncovered, SurvivingPhase };
struct UncoveredRefillSupport {
    std::vector<RefillCell> boundary;
    std::size_t visitedUncovered = 0;
    bool enclosed() const { return boundary.empty(); }
};

// Used only when direct refill support is absent. Traverse newly uncovered
// cells, never surviving solids or pre-existing gas/liquid. Return the nearest
// graph shell of surviving phase cells, each exactly once. An exhausted search
// identifies a newly opened, isolated pore containing no resolved liquid mass.
template<class Classify>
UncoveredRefillSupport findUncoveredRefillSupport(
    const RefillCell& start, const std::vector<RefillCell>& directions, Classify classify) {
    if(classify(start)!=RefillCellKind::NewlyUncovered || directions.empty())
        throw std::invalid_argument("Invalid newly uncovered refill seed/stencil");
    std::set<RefillCell> visited{start};
    std::vector<RefillCell> frontier{start};
    while(!frontier.empty()) {
        std::vector<RefillCell> next;
        std::set<RefillCell> boundarySeen;
        UncoveredRefillSupport result;
        for(const auto& cell:frontier)for(const auto& direction:directions) {
            RefillCell neighbor{cell[0]+direction[0],cell[1]+direction[1],cell[2]+direction[2]};
            const auto kind=classify(neighbor);
            if(kind==RefillCellKind::SurvivingPhase) {
                if(boundarySeen.insert(neighbor).second)result.boundary.push_back(neighbor);
            } else if(kind==RefillCellKind::NewlyUncovered && visited.insert(neighbor).second)
                next.push_back(neighbor);
        }
        result.visitedUncovered=visited.size();
        if(!result.boundary.empty())return result;
        frontier=std::move(next);
    }
    return {{},visited.size()};
}
} // namespace capillary_coupling
