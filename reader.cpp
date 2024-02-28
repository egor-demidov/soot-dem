//
// Created by egor on 2/28/24.
//

#include <iostream>
#include <fstream>

#include "reader.h"

std::vector<Eigen::Vector3d> load_mackowski_aggregate(std::string const & path, double r_part) {
    std::ifstream ifs(path);

    if (!ifs.good()) {
        std::cerr << "Unable to read the aggregate file: " << path << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::vector<Eigen::Vector3d> x0;

    while (getline(ifs, line)) {
        if (!line.empty()) {
            std::istringstream oss(line);
            double _, x, y, z;
            oss >> _ >> x >> y >> z;
            x0.emplace_back(x * r_part, y * r_part, z * r_part);
        }
    }

    return x0;
}

std::vector<Eigen::Vector3d> load_vtk_aggregate(std::string const & path, double r_part) {
    std::ifstream ifs(path);

    if (!ifs.good()) {
        std::cerr << "Unable to read the aggregate file: " << path << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string _;
    size_t num_pts;

    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    ifs >> _ >> num_pts >> _;

    std::vector<Eigen::Vector3d> out;
    out.resize(num_pts);

    double x, y, z;
    for (size_t i = 0; i < num_pts; i ++) {
        ifs >> x >> y >> z;
        out[i] = Eigen::Vector3d {x, y, z} * r_part;
    }

    return out;
}
