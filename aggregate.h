//
// Created by egor on 2/23/24.
//

#ifndef SOOT_AFM_AGGREGATE_H
#define SOOT_AFM_AGGREGATE_H

#include <fstream>

#include <libgran/sinter_bridge/alt_sinter_bridge.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>

template <typename field_value_t, typename real_t>
struct aggregate {

    typedef contact_force_functor<field_value_t, real_t> contact_model_t;
    typedef alt_sinter_functor<field_value_t, real_t> sinter_model_t;
    typedef hamaker_functor<field_value_t, real_t> hamaker_model_t;

    aggregate(
              /* PARAMETERS FOR THE FRICTION MODEL */
              real_t k_contact,                 // Normal stiffness coefficient
              real_t gamma_n_contact,           // Normal damping coefficient
              real_t k_t_contact,               // Stiffness coefficient for sticking/sliding
              real_t gamma_t_contact,           // Damping coefficient for sticking/sliding
              real_t mu_s_contact,              // Static friction coefficient for sticking/sliding
              real_t phi_d_contact,             // Coulomb coefficient for sticking/sliding
              real_t k_r_contact,               // Stiffness coefficient for rolling
              real_t gamma_r_contact,           // Damping coefficient for rolling
              real_t mu_r_contact,              // Static friction coefficient for rolling
              real_t phi_r_contact,             // Coulomb coefficient for rolling
              real_t k_o_contact,               // Stiffness for torsion
              real_t gamma_o_contact,           // Damping coefficient for torsion
              real_t mu_o_contact,              // Static friction coefficient for torsion
              real_t phi_o_contact,             // Coulomb coefficient for torsion

              /* PARAMETERS FOR THE BOND MODEL */
              real_t k_bond,                 // Normal stiffness coefficient
              real_t gamma_n_bond,           // Normal damping coefficient
              real_t k_t_bond,               // Stiffness coefficient for sticking/sliding
              real_t gamma_t_bond,           // Damping coefficient for sticking/sliding
              real_t k_r_bond,               // Stiffness coefficient for rolling
              real_t gamma_r_bond,           // Damping coefficient for rolling
              real_t k_o_bond,               // Stiffness for torsion
              real_t gamma_o_bond,           // Damping coefficient for torsion
              real_t critical_separation,

              /* PARAMETERS FOR THE HAMAKER MODEL */
              real_t A,
              real_t h0,

              /* GENERAL PARAMETERS */
              std::vector<field_value_t> x0,         // Initial positions
              size_t n_part,            // Number of particles in the system
              real_t r_part,            // Radius of a particle
              real_t mass,              // Mass
              real_t inertia,           // Moment of inertia
              real_t dt,                // Time step for spring update (same as integration time step for 1st order schemes)
              field_value_t field_zero, // Zero-valued field_value_t
              real_t real_zero) :

    contact_model {
        n_part, k_contact, gamma_n_contact, k_t_contact, gamma_t_contact, mu_s_contact, phi_d_contact,
        k_r_contact, gamma_r_contact, mu_r_contact, phi_r_contact,
        k_o_contact, gamma_o_contact, mu_o_contact, phi_o_contact,
        r_part, mass, inertia, dt, field_zero, real_zero
    }, sinter_model {
        n_part, std::move(x0), k_bond, gamma_n_bond, k_t_bond, gamma_t_bond,
        k_r_bond, gamma_r_bond, k_o_bond, gamma_o_bond,
        r_part, mass, inertia, dt, field_zero, real_zero, critical_separation, contact_model
    }, hamaker_model {
        A, h0, r_part, mass, field_zero, real_zero
    } {}

    std::pair<field_value_t, field_value_t> operator () (size_t i, size_t j,
            std::vector<field_value_t> const & x,
             std::vector<field_value_t> const & v,
             std::vector<field_value_t> const & theta,
             std::vector<field_value_t> const & omega,
             real_t t) {


        auto result = sinter_model(i, j, x, v, theta, omega, t);
        result += hamaker_model(i, j, x, v, theta, omega, t);
        return result;
    }

    [[nodiscard]]
    std::vector<bool> const & get_bonded_contacts() const {
        return sinter_model.bonded_contacts;
    }

    [[deprecated]] // TODO: Later needs to be replaced with high level functions that handle neck breakage
    std::vector<bool> & get_bonded_contacts() {
        return sinter_model.bonded_contacts;
    }

private:

    /*void load_mackowski_aggregate(std::string const & path) {
        std::ifstream ifs(path);

        if (!ifs.good()) {
            std::cerr << "Unable to read the aggregate file: " << path << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string line;

        while (getline(ifs, line)) {
            if (!line.empty()) {
                std::istringstream oss(line);
                double _, x, y, z;
                oss >> _ >> x >> y >> z;
                x0.emplace_back(x * r_part, y * r_part, z * r_part);
            }
        }
    }*/

    contact_model_t contact_model;
    sinter_model_t sinter_model;
    hamaker_model_t hamaker_model;
};

#endif //SOOT_AFM_AGGREGATE_H
