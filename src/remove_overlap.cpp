/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include "remove_overlap.h"

void remove_overlap(std::vector<Eigen::Vector3d> & x, double r_part, double d_crit, size_t n_iter) {

    // Build a neighbor list
    std::vector<std::pair<size_t, size_t>> neighbor_list;

    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            double d = (x[j] - x[i]).norm() - 2.0 * r_part;
            if (abs(d) <= d_crit)
                neighbor_list.emplace_back(i, j);
        }
    }

    // Iteratively remove overlaps
    for (size_t i = 0; i < n_iter; i ++) {
        for (auto const & pair : neighbor_list) {
            auto & part_i = x[pair.first];
            auto & part_j = x[pair.second];

            Eigen::Vector3d n = (part_j - part_i).normalized();
            Eigen::Vector3d dist = ((part_j - part_i).norm() - 2.0 * r_part) * n;

            // Move each particle by an equal amount so that they are in point-touch contact
            part_i += dist / 2.0;
            part_j -= dist / 2.0;
        }
    }
}

