//
// Created by egor on 3/6/24.
//

#ifndef SOOT_AFM_IO_COMMON_H
#define SOOT_AFM_IO_COMMON_H

#include <vector>

#include <Eigen/Eigen>

#include "parameter_loader.h"

std::vector<Eigen::Vector3d> load_aggregate(parameter_store_t const & parameter_store);

#endif //SOOT_AFM_IO_COMMON_H
