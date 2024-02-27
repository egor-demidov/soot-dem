//
// Created by egor on 2/26/24.
//

#include <fstream>
#include <iostream>

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
