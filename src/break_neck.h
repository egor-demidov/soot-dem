/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_BREAK_NECK_H
#define SOOT_AFM_BREAK_NECK_H

#include <vector>
#include <iostream>

#include "aggregate.h"

void break_random_neck(std::vector<bool> & bonded_contacts, size_t n_part);

template <typename field_value_t, typename real_t>
void break_strained_necks(aggregate<field_value_t, real_t> & aggregate_model,
                          std::vector<field_value_t> const & x,
                          real_t k_n_bond,
                          real_t k_t_bond,
                          real_t k_r_bond,
                          real_t k_o_bond,
                          real_t r_part) {

    const auto n_part = x.size();

    for (size_t i = 0; i < n_part - 1; i ++) {
        for (size_t j = i + 1; j < n_part; j ++) {
            bool bond = aggregate_model.get_bonded_contacts()[i * n_part + j];

            if (!bond)
                continue;

            real_t xi_n = (x[i] - x[j]).norm() - 2.0 * r_part;
            auto [xi_t, xi_r, xi_o] = aggregate_model.get_sinter_model().get_contact_springs()[i * n_part + j];

            real_t total_strain_energy = k_t_bond * xi_t.dot(xi_t)
                                        + k_r_bond * xi_r.dot(xi_r)
                                        + k_o_bond * xi_o.dot(xi_o);

            // Normal component only contributes in case of tension
            if (xi_n > 0.0)
                total_strain_energy += k_n_bond * xi_n * xi_n;

            constexpr double e_crit = 1e-20;

            if (total_strain_energy < e_crit)
                continue;

            std::cout << "Breaking neck" << std::endl;

            aggregate_model.get_bonded_contacts()[i * n_part + j] = false;
            aggregate_model.get_bonded_contacts()[j * n_part + i] = false;
        }
    }
}

#endif //SOOT_AFM_BREAK_NECK_H
