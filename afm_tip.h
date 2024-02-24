//
// Created by egor on 2/23/24.
//

#ifndef SOOT_AFM_AFM_TIP_H
#define SOOT_AFM_AFM_TIP_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <functional>

#include <libgran/contact_force/surface_contact_force.h>
#include <libgran/hamaker_force/surface_hamaker_force.h>
#include <libgran/surface_force/triangular_facet.h>
#include <libgran/granular_system/granular_system.h>

template <typename field_value_t, typename real_t>
struct afm_tip {

    typedef surface_contact_force_functor<field_value_t, real_t> surface_contact_force_t;
    typedef surface_hamaker_functor<field_value_t, real_t> surface_hamaker_t;
    typedef triangular_facet<field_value_t, real_t, surface_contact_force_t, surface_hamaker_t> facet_t;

    afm_tip(std::tuple<field_value_t, field_value_t, field_value_t> base_vertices,
            field_value_t peak_vertex,
            field_value_t tip_velocity,

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
           real_t real_zero) :
    force_functors_facets {
        std::pair<surface_contact_force_t, surface_hamaker_t>{
            {n_part, k, gamma_n, k_t, gamma_t, mu_s, phi_d, k_r, gamma_r, mu_r, phi_r, k_o, gamma_o, mu_o, phi_o, r_part, mass, inertia, dt, field_zero, real_zero},
            {A, h0, r_part, mass, field_zero, real_zero}
        }, std::pair<surface_contact_force_t, surface_hamaker_t>{
            {n_part, k, gamma_n, k_t, gamma_t, mu_s, phi_d, k_r, gamma_r, mu_r, phi_r, k_o, gamma_o, mu_o, phi_o, r_part, mass, inertia, dt, field_zero, real_zero},
            {A, h0, r_part, mass, field_zero, real_zero}
        }, std::pair<surface_contact_force_t, surface_hamaker_t>{
            {n_part, k, gamma_n, k_t, gamma_t, mu_s, phi_d, k_r, gamma_r, mu_r, phi_r, k_o, gamma_o, mu_o, phi_o, r_part, mass, inertia, dt, field_zero, real_zero},
            {A, h0, r_part, mass, field_zero, real_zero}
        }
    }, r_part {r_part}, mass {mass}, facets {
       facet_t {
            tip_velocity, {std::get<0>(base_vertices), std::get<1>(base_vertices), peak_vertex}, force_functors_facets[0].first, force_functors_facets[0].second
        }, facet_t {
            tip_velocity, {std::get<1>(base_vertices), std::get<2>(base_vertices), peak_vertex}, force_functors_facets[1].first, force_functors_facets[1].second
        }, facet_t {
            tip_velocity, {std::get<2>(base_vertices), std::get<0>(base_vertices), peak_vertex}, force_functors_facets[2].first, force_functors_facets[2].second
        }
    }, force_accumulator {real_zero}, accumulate_forces {false} {

        field_value_t u = std::get<1>(base_vertices) - std::get<0>(base_vertices);
        field_value_t v = std::get<2>(base_vertices) - std::get<0>(base_vertices);

        force_normal_dir = u.cross(v).normalized();
    }

    bool toggle_force_accumulation() {
        accumulate_forces = !accumulate_forces;
        return accumulate_forces;
    }

    std::pair<field_value_t, field_value_t> operator () (size_t i,
                                                         std::vector<field_value_t> const & x,
                                                         std::vector<field_value_t> const & v,
                                                         std::vector<field_value_t> const & theta,
                                                         std::vector<field_value_t> const & omega,
                                                         real_t t) {

        auto result = facets[0](i, x, v, theta, omega, t);
        result += facets[1](i, x, v, theta, omega, t);
        result += facets[2](i, x, v, theta, omega, t);

        if (accumulate_forces) {
            force_accumulator += result.first.dot(force_normal_dir) * mass;
        }

        return result;
    }

    // Write mesh to an STL file
    void dump_mesh(std::string const & dir, size_t count) {
        std::stringstream ss;
        ss << dir << "/tip_" << count << ".stl";

        std::ofstream ofs(ss.str());

        if (!ofs.good()) {
            std::cerr << "Unable to create a facet file" << std::endl;
            exit(EXIT_FAILURE);
        }

        ofs << "solid facet\n";

        for (size_t i = 0; i < std::size(facets); i ++) {
            auto normal = facets[i].get_unit_normal();
            auto vertices = facets[i].get_vertices();

            ofs << "facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
            ofs << "\touter loop\n";
            ofs << "\t\tvertex " << std::get<0>(vertices)[0] / r_part << " " << std::get<0>(vertices)[1] / r_part << " " << std::get<0>(vertices)[2] / r_part << "\n";
            ofs << "\t\tvertex " << std::get<1>(vertices)[0] / r_part << " " << std::get<1>(vertices)[1] / r_part << " " << std::get<1>(vertices)[2] / r_part << "\n";
            ofs << "\t\tvertex " << std::get<2>(vertices)[0] / r_part << " " << std::get<2>(vertices)[1] / r_part << " " << std::get<2>(vertices)[2] / r_part << "\n";
            ofs << "\tendloop\n";
            ofs << "endfacet\n";
        }

        ofs << "endsolid facet\n";
    }

    void update_positions(real_t dt) {
        for (size_t i = 0; i < std::size(facets); i ++) {
            facets[i].update_positions(dt);
        }
    }

    void update_velocities(field_value_t const & v) {
        for (size_t i = 0; i < std::size(facets); i ++) {
            facets[i].v_facet = v;
        }
    }

    real_t force_accumulator;

private:
    // 7 facets total
    std::array<std::pair<surface_contact_force_t, surface_hamaker_t>, 3> force_functors_facets;
    const real_t r_part, mass;
    std::array<facet_t, 3> facets;
    field_value_t force_normal_dir;
    bool accumulate_forces;
};

#endif //SOOT_AFM_AFM_TIP_H
