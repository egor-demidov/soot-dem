// Created by egor on 2/28/24.

#ifndef SOOT_AFM_CENTRAL_FORCE_H
#define SOOT_AFM_CENTRAL_FORCE_H

#include <vector>

template <typename field_value_t, typename real_t>
struct central_force_functor {
    central_force_functor(field_value_t const & center_of_mass, real_t magnitude, real_t mass, field_value_t zero_field) :
        magnitude{magnitude}, mass{mass}, zero_field{zero_field}, center_of_mass{center_of_mass} {}

    void update_center_of_mass(field_value_t const & new_center_of_mass) {
        center_of_mass = new_center_of_mass;
    }

    std::pair<field_value_t, field_value_t> operator () (size_t i,
            std::vector<field_value_t> const & x,
         std::vector<field_value_t> const & v [[maybe_unused]],
         std::vector<field_value_t> const & theta [[maybe_unused]],
         std::vector<field_value_t> const & omega [[maybe_unused]],
         real_t t [[maybe_unused]]) const {

        field_value_t n = (center_of_mass - x[i]).normalized();
        field_value_t f = magnitude * n;

        return std::make_pair(f / mass, zero_field);
    }

private:
    const real_t magnitude, mass;
    const field_value_t zero_field;
    field_value_t center_of_mass;
};

#endif //SOOT_AFM_CENTRAL_FORCE_H
