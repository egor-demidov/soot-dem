//
// Created by egor on 2/26/24.
//

#ifndef SOOT_AFM_READER_H
#define SOOT_AFM_READER_H

#include <vector>
#include <string>

#include <Eigen/Eigen>

std::vector<Eigen::Vector3d> load_mackowski_aggregate(std::string const & path, double r_part);

#endif //SOOT_AFM_READER_H
