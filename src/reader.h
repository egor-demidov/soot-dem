/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_READER_H
#define SOOT_AFM_READER_H

#include <vector>
#include <string>
#include <filesystem>

#include <Eigen/Eigen>

std::vector<Eigen::Vector3d> load_mackowski_aggregate(std::filesystem::path const & path, double r_part);
std::vector<Eigen::Vector3d> load_vtk_aggregate(std::filesystem::path const & path, double r_part);
std::vector<Eigen::Vector3d> load_flage_aggregate(std::filesystem::path const & path, double r_part);
std::vector<bool> load_necks(std::filesystem::path const & path, size_t n_part);
std::vector<double> load_neck_strengths(std::filesystem::path const & path);

#endif //SOOT_AFM_READER_H
