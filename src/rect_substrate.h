/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef RECT_SUBSTRATE_H
#define RECT_SUBSTRATE_H

#include <iostream>
#include <sstream>
#include <fstream>

#include <libgran/contact_force/surface_contact_force.h>
#include <libgran/hamaker_force/surface_hamaker_force.h>
#include <libgran/surface_force/triangular_facet.h>
#include <libgran/granular_system/granular_system.h>

#include "coating_force.h"

template <typename field_value_t, typename real_t>
struct rect_substrate {

    typedef surface_contact_force_functor<field_value_t, real_t> surface_contact_force_t;
    typedef surface_hamaker_functor<field_value_t, real_t> surface_hamaker_t;
    typedef triangular_facet<field_value_t, real_t, surface_contact_force_t, surface_hamaker_t> facet_t;

    rect_substrate(std::tuple<field_value_t, field_value_t, field_value_t, field_value_t> const & vertices,
                    size_t n_part,            // Number of particles in the system
                      real_t k,                 // Normal stiffness coefficient
                      real_t gamma_n,           // Normal damping coefficient
                      real_t k_t,               // Stiffness coefficient for sticking/sliding
                      real_t gamma_t,           // Damping coefficient for sticking/sliding
                      real_t mu_s,              // Static friction coefficient for sticking/sliding
                      real_t phi_d,             // Coulomb coefficient for sticking/sliding
                      real_t k_r,               // Stiffness coefficient for rolling
                      real_t gamma_r,           // Damping coefficient for rolling
                      real_t mu_r,              // Static friction coefficient for rolling
                      real_t phi_r,             // Coulomb coefficient for rolling
                      real_t k_o,               // Stiffness for torsion
                      real_t gamma_o,           // Damping coefficient for torsion
                      real_t mu_o,              // Static friction coefficient for torsion
                      real_t phi_o,             // Coulomb coefficient for torsion
                      real_t A,                 // Hamaker constant
                      real_t h0,                // Hamaker satruation distance
                      real_t r_part,            // Radius of a particle
                      real_t mass,              // Mass
                      real_t inertia,           // Moment of inertia
                      real_t dt,                // Time step for spring update (same as integration time step for 1st order schemes)
                      field_value_t field_zero, // Zero-valued field_value_t
                      real_t real_zero) : force_functors_facet_1 {
                          {n_part, k, gamma_n, k_t, gamma_t, mu_s, phi_d, k_r, gamma_r, mu_r, phi_r, k_o, gamma_o, mu_o, phi_o, r_part, mass, inertia, dt, field_zero, real_zero},
                          {A, h0, r_part, mass, field_zero, real_zero}
                      }, force_functors_facet_2 {
                          {n_part, k, gamma_n, k_t, gamma_t, mu_s, phi_d, k_r, gamma_r, mu_r, phi_r, k_o, gamma_o, mu_o, phi_o, r_part, mass, inertia, dt, field_zero, real_zero},
                          {A, h0, r_part, mass, field_zero, real_zero}
                      }, field_zero(std::move(field_zero)), r_part {r_part}, facet_1 {
                          field_zero, {std::get<0>(vertices), std::get<1>(vertices), std::get<2>(vertices)}, force_functors_facet_1.first, force_functors_facet_1.second
                      }, facet_2 {
                          field_zero, {std::get<2>(vertices), std::get<3>(vertices), std::get<0>(vertices)}, force_functors_facet_2.first, force_functors_facet_2.second
                      } {}

    std::pair<field_value_t, field_value_t> operator () (size_t i,
        std::vector<field_value_t> const & x,
        std::vector<field_value_t> const & v,
        std::vector<field_value_t> const & theta,
        std::vector<field_value_t> const & omega,
        real_t t) {

        auto result = facet_1(i, x, v, theta, omega, t);
        result += facet_2(i, x, v, theta, omega, t);

        return result;
    }

    // Write mesh to an STL file
    void dump_mesh(std::string const & dir, size_t count) {
        std::stringstream ss;
        ss << dir << "/substrate_" << count << ".stl";

        std::ofstream ofs(ss.str());

        if (!ofs.good()) {
            std::cerr << "Unable to create a facet file" << std::endl;
            exit(EXIT_FAILURE);
        }

        ofs << "solid facet\n";

        auto normal = facet_1.get_unit_normal();
        auto vertices = facet_1.get_vertices();

        ofs << "facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
        ofs << "\touter loop\n";
        ofs << "\t\tvertex " << std::get<0>(vertices)[0] / r_part << " " << std::get<0>(vertices)[1] / r_part << " " << std::get<0>(vertices)[2] / r_part << "\n";
        ofs << "\t\tvertex " << std::get<1>(vertices)[0] / r_part << " " << std::get<1>(vertices)[1] / r_part << " " << std::get<1>(vertices)[2] / r_part << "\n";
        ofs << "\t\tvertex " << std::get<2>(vertices)[0] / r_part << " " << std::get<2>(vertices)[1] / r_part << " " << std::get<2>(vertices)[2] / r_part << "\n";
        ofs << "\tendloop\n";
        ofs << "endfacet\n";

        normal = facet_2.get_unit_normal();
        vertices = facet_2.get_vertices();

        ofs << "facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
        ofs << "\touter loop\n";
        ofs << "\t\tvertex " << std::get<0>(vertices)[0] / r_part << " " << std::get<0>(vertices)[1] / r_part << " " << std::get<0>(vertices)[2] / r_part << "\n";
        ofs << "\t\tvertex " << std::get<1>(vertices)[0] / r_part << " " << std::get<1>(vertices)[1] / r_part << " " << std::get<1>(vertices)[2] / r_part << "\n";
        ofs << "\t\tvertex " << std::get<2>(vertices)[0] / r_part << " " << std::get<2>(vertices)[1] / r_part << " " << std::get<2>(vertices)[2] / r_part << "\n";
        ofs << "\tendloop\n";
        ofs << "endfacet\n";

        ofs << "endsolid facet\n";
    }

private:
    std::pair<surface_contact_force_t, surface_hamaker_t> force_functors_facet_1, force_functors_facet_2;
    const field_value_t field_zero;
    const real_t r_part;
    facet_t facet_1, facet_2;
};

template <typename field_value_t, typename real_t>
struct rect_substrate_with_coating {

    typedef surface_coating_functor<field_value_t, real_t> surface_coating_force_t;
    typedef surface_contact_force_functor<field_value_t, real_t> surface_contact_force_t;
    typedef surface_hamaker_functor<field_value_t, real_t> surface_hamaker_t;
    typedef triangular_facet<field_value_t, real_t, surface_contact_force_t, surface_hamaker_t, surface_coating_force_t> facet_t;

    rect_substrate_with_coating(std::tuple<field_value_t, field_value_t, field_value_t, field_value_t> const & vertices,
                   size_t n_part,            // Number of particles in the system
                   real_t k,                 // Normal stiffness coefficient
                   real_t gamma_n,           // Normal damping coefficient
                   real_t k_t,               // Stiffness coefficient for sticking/sliding
                   real_t gamma_t,           // Damping coefficient for sticking/sliding
                   real_t mu_s,              // Static friction coefficient for sticking/sliding
                   real_t phi_d,             // Coulomb coefficient for sticking/sliding
                   real_t k_r,               // Stiffness coefficient for rolling
                   real_t gamma_r,           // Damping coefficient for rolling
                   real_t mu_r,              // Static friction coefficient for rolling
                   real_t phi_r,             // Coulomb coefficient for rolling
                   real_t k_o,               // Stiffness for torsion
                   real_t gamma_o,           // Damping coefficient for torsion
                   real_t mu_o,              // Static friction coefficient for torsion
                   real_t phi_o,             // Coulomb coefficient for torsion
                   real_t A,                 // Hamaker constant
                   real_t h0,                // Hamaker satruation distance
                   real_t r_part,            // Radius of a particle
                   real_t mass,              // Mass
                   real_t inertia,           // Moment of inertia
                   real_t dt,                // Time step for spring update (same as integration time step for 1st order schemes)
                    real_t d_cutoff,         // Coating force cutoff distance
                    real_t magnitude,        // Coating force maximum magnitude
                    real_t drop_rate,        // Coating force drop rate
                   field_value_t field_zero, // Zero-valued field_value_t
                   real_t real_zero) : force_functors_facet_1 {
            {n_part, k, gamma_n, k_t, gamma_t, mu_s, phi_d, k_r, gamma_r, mu_r, phi_r, k_o, gamma_o, mu_o, phi_o, r_part, mass, inertia, dt, field_zero, real_zero},
            {A, h0, r_part, mass, field_zero, real_zero},
            {d_cutoff, magnitude, drop_rate, mass, field_zero}
    }, force_functors_facet_2 {
            {n_part, k, gamma_n, k_t, gamma_t, mu_s, phi_d, k_r, gamma_r, mu_r, phi_r, k_o, gamma_o, mu_o, phi_o, r_part, mass, inertia, dt, field_zero, real_zero},
            {A, h0, r_part, mass, field_zero, real_zero},
            {d_cutoff, magnitude, drop_rate, mass, field_zero}
    }, field_zero(std::move(field_zero)), r_part {r_part}, facet_1 {
            field_zero, {std::get<0>(vertices), std::get<1>(vertices), std::get<2>(vertices)}, std::get<0>(force_functors_facet_1), std::get<1>(force_functors_facet_1), std::get<2>(force_functors_facet_1)
    }, facet_2 {
            field_zero, {std::get<2>(vertices), std::get<3>(vertices), std::get<0>(vertices)}, std::get<0>(force_functors_facet_2), std::get<1>(force_functors_facet_2), std::get<2>(force_functors_facet_2)
    } {}

    std::pair<field_value_t, field_value_t> operator () (size_t i,
                                                         std::vector<field_value_t> const & x,
                                                         std::vector<field_value_t> const & v,
                                                         std::vector<field_value_t> const & theta,
                                                         std::vector<field_value_t> const & omega,
                                                         real_t t) {

        auto result = facet_1(i, x, v, theta, omega, t);
        result += facet_2(i, x, v, theta, omega, t);

        return result;
    }

    // Write mesh to an STL file
    void dump_mesh(std::string const & dir, size_t count) {
        std::stringstream ss;
        ss << dir << "/substrate_" << count << ".stl";

        std::ofstream ofs(ss.str());

        if (!ofs.good()) {
            std::cerr << "Unable to create a facet file" << std::endl;
            exit(EXIT_FAILURE);
        }

        ofs << "solid facet\n";

        auto normal = facet_1.get_unit_normal();
        auto vertices = facet_1.get_vertices();

        ofs << "facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
        ofs << "\touter loop\n";
        ofs << "\t\tvertex " << std::get<0>(vertices)[0] / r_part << " " << std::get<0>(vertices)[1] / r_part << " " << std::get<0>(vertices)[2] / r_part << "\n";
        ofs << "\t\tvertex " << std::get<1>(vertices)[0] / r_part << " " << std::get<1>(vertices)[1] / r_part << " " << std::get<1>(vertices)[2] / r_part << "\n";
        ofs << "\t\tvertex " << std::get<2>(vertices)[0] / r_part << " " << std::get<2>(vertices)[1] / r_part << " " << std::get<2>(vertices)[2] / r_part << "\n";
        ofs << "\tendloop\n";
        ofs << "endfacet\n";

        normal = facet_2.get_unit_normal();
        vertices = facet_2.get_vertices();

        ofs << "facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
        ofs << "\touter loop\n";
        ofs << "\t\tvertex " << std::get<0>(vertices)[0] / r_part << " " << std::get<0>(vertices)[1] / r_part << " " << std::get<0>(vertices)[2] / r_part << "\n";
        ofs << "\t\tvertex " << std::get<1>(vertices)[0] / r_part << " " << std::get<1>(vertices)[1] / r_part << " " << std::get<1>(vertices)[2] / r_part << "\n";
        ofs << "\t\tvertex " << std::get<2>(vertices)[0] / r_part << " " << std::get<2>(vertices)[1] / r_part << " " << std::get<2>(vertices)[2] / r_part << "\n";
        ofs << "\tendloop\n";
        ofs << "endfacet\n";

        ofs << "endsolid facet\n";
    }

private:
    std::tuple<surface_contact_force_t, surface_hamaker_t, surface_coating_force_t> force_functors_facet_1, force_functors_facet_2;
    const field_value_t field_zero;
    const real_t r_part;
    facet_t facet_1, facet_2;
};

#endif //RECT_SUBSTRATE_H
