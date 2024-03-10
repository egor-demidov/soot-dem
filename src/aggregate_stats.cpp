//
// Created by egor on 3/4/24.
//

#include "aggregate_stats.h"

Eigen::Vector3d center_of_mass(std::vector<Eigen::Vector3d> const & x) {
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    std::for_each(x.begin(), x.end(), [&acc] (auto const & val) {
        acc += val;
    });
    return acc / double(x.size());
}

double r_gyration(std::vector<Eigen::Vector3d> const & x) {
    Eigen::Vector3d x0 = center_of_mass(x);
    double result = 0.0;
    std::for_each(x.begin(), x.end(), [&result, x0] (auto const & val) {
        Eigen::Vector3d dist = val - x0;
        result += dist.dot(dist);
    });
    return sqrt(result / double(x.size()));
}
