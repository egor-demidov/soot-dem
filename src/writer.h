//
// Created by egor on 2/24/24.
//

#ifndef SOOT_AFM_WRITER_H
#define SOOT_AFM_WRITER_H

#include <vector>
#include <string>

#include <Eigen/Eigen>

/*
struct system_representation_t  {

    // GENERAL ATTRIBUTES
    size_t n_part;
    double r_part;
    double mass;
    double inertia;
    double dt;

    // FRICTIONAL CONTACT ATTRIBUTES
    double k;
    double gamma_n;
    double k_t;
    double gamma_t;
    double mu_s;
    double phi_d;
    double k_r;
    double gamma_r;
    double mu_r;
    double phi_r;
    double k_o;
    double gamma_o;
    double mu_o;
    double phi_o;

    // BONDED CONTACT ATTRIBUTES
    double k_bond;                 // Normal stiffness coefficient
    double gamma_n_bond;           // Normal damping coefficient
    double k_t_bond;               // Stiffness coefficient for sticking/sliding
    double gamma_t_bond;           // Damping coefficient for sticking/sliding
    double k_r_bond;               // Stiffness coefficient for rolling
    double gamma_r_bond;           // Damping coefficient for rolling
    double k_o_bond;               // Stiffness for torsion
    double gamma_o_bond;           // Damping coefficient for torsion

    // VAN DER WAALS ATTRIBUTES
    double A;                       // Hamaker constant
    double h0;                      // Hamaker saturation distance

};

struct vector_representation {
    double x, y, z;
};
 */

void dump_necks(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                std::vector<bool> const & bonded_contacts, double r_part);

// Overload to write particle positions without velocities and orientations
void dump_particles(std::string const & name, std::vector<Eigen::Vector3d> const & x, double r_part);

void dump_particles(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                    std::vector<Eigen::Vector3d> const & v,
                    std::vector<Eigen::Vector3d> const & a,
                    std::vector<Eigen::Vector3d> const & omega,
                    std::vector<Eigen::Vector3d> const & alpha,
                    double r_part);

#endif //SOOT_AFM_WRITER_H
