//
// Created by egor on 2/24/24.
//

#ifndef SOOT_AFM_WRITER_H
#define SOOT_AFM_WRITER_H

#include <vector>
#include <string>

#include <Eigen/Eigen>

void dump_necks(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                std::vector<bool> const & bonded_contacts, double r_part);

void dump_particles(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x, double r_part);

#endif //SOOT_AFM_WRITER_H
