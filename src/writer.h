/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_WRITER_H
#define SOOT_AFM_WRITER_H

#include <vector>
#include <string>

#include <Eigen/Eigen>

bool dump_necks(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                std::vector<bool> const & bonded_contacts, double r_part);

// Overload to write necks with energies
bool dump_necks(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                std::vector<bool> const & bonded_contacts, double r_part, double k_n_bond,
                double k_t_bond, double k_r_bond, double k_o_bond,
                std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> const & contact_springs,
                size_t n_part);

// Overload to write particle positions without velocities and orientations
bool dump_particles(std::string const & name, std::vector<Eigen::Vector3d> const & x, double r_part);

bool dump_particles(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                    std::vector<Eigen::Vector3d> const & v,
                    std::vector<Eigen::Vector3d> const & a,
                    std::vector<Eigen::Vector3d> const & omega,
                    std::vector<Eigen::Vector3d> const & alpha,
                    double r_part);

bool dump_particles(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                    std::vector<Eigen::Vector3d> const & theta,
                    std::vector<Eigen::Vector3d> const & v,
                    std::vector<Eigen::Vector3d> const & a,
                    std::vector<Eigen::Vector3d> const & omega,
                    std::vector<Eigen::Vector3d> const & alpha,
                    double r_part);

#endif //SOOT_AFM_WRITER_H
