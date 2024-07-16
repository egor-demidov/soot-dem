/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_READER_H
#define SOOT_AFM_READER_H

#include <vector>
#include <string>

#include <Eigen/Eigen>

std::vector<Eigen::Vector3d> load_mackowski_aggregate(std::string const & path, double r_part);
std::vector<Eigen::Vector3d> load_vtk_aggregate(std::string const & path, double r_part);
std::vector<Eigen::Vector3d> load_flage_aggregate(std::string const & path, double r_part);

#endif //SOOT_AFM_READER_H
