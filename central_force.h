//
// Created by egor on 2/26/24.
//

#ifndef SOOT_AFM_CENTRAL_FORCE_H
#define SOOT_AFM_CENTRAL_FORCE_H

template <typename field_value_t, typename real_t>
struct central_force_functor {
    central_force_functor(std::vector<field_value_t> const & x, real_t k, real_t mass, field_value_t zero_field) :
        k{k}, mass{mass}, zero_field{zero_field} {

        update_center_of_mass(x);
    }

    std::pair<field_value_t, field_value_t> operator () (size_t i,
            std::vector<field_value_t> const & x,
             std::vector<field_value_t> const & v [[maybe_unused]],
             std::vector<field_value_t> const & theta [[maybe_unused]],
             std::vector<field_value_t> const & omega [[maybe_unused]],
             real_t t [[maybe_unused]]) const {

        return std::make_pair(k * (center_of_mass - x[i]) / mass, zero_field);
    }

    void update_center_of_mass(std::vector<field_value_t> const & x) {
        center_of_mass = zero_field;
        for (auto const & pt : x) {
            center_of_mass += pt;
        }

        center_of_mass /= real_t(x.size());
    }

private:
    field_value_t center_of_mass;
    const real_t k, mass;
    const field_value_t zero_field;
};

#endif //SOOT_AFM_CENTRAL_FORCE_H
