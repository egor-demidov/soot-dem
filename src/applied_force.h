//
// Created by egor on 9/24/24.
//

#ifndef APPLIED_FORCE_H
#define APPLIED_FORCE_H

// Unary coating functor can be used to model capillary force between a static plane and particles
template <typename field_value_t, typename real_t>
struct unary_applied_force_functor {
    unary_applied_force_functor(
        size_t part_i_index,
        size_t part_j_index,
        real_t magnitude,
        real_t mass,
        field_value_t field_zero
    )
        : part_i_index{part_i_index}
        , part_j_index{part_j_index}
        , magnitude{magnitude}
        , mass{mass}
        , field_zero{field_zero}
    {}

    std::pair<field_value_t, field_value_t> operator () (
        size_t i,
        std::vector<field_value_t> const & x,
        std::vector<field_value_t> const & v,
        std::vector<field_value_t> const & theta,
        std::vector<field_value_t> const & omega,
        real_t t) const {

        field_value_t f = field_zero;

        if (i == part_i_index) [[unlikely]] {
            field_value_t dir = (x[part_i_index] - x[part_j_index]).normalized();
            f = dir * magnitude;
        } else if (i == part_j_index) [[unlikely]] {
            field_value_t dir = (x[part_j_index] - x[part_i_index]).normalized();
            f = dir * magnitude;
        }

        return std::make_pair(f/mass, field_zero);
    }

private:
    const size_t part_i_index, part_j_index;
    const real_t magnitude, mass;
    const field_value_t field_zero;
};

#endif //APPLIED_FORCE_H
