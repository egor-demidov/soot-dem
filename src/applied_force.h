//
// Created by egor on 9/24/24.
//

#ifndef APPLIED_FORCE_H
#define APPLIED_FORCE_H

// Unary coating functor can be used to model capillary force between a static plane and particles

template <typename real_t>
struct ConstantAppliedForce {
public:
    explicit ConstantAppliedForce(real_t magnitude)
        : magnitude_(magnitude)
    {}

    real_t operator () (real_t t [[maybe_unused]]) const noexcept {
        return magnitude_;
    }

private:
    const real_t magnitude_;
};

template <typename field_value_t, typename real_t, typename f_of_t = ConstantAppliedForce<real_t>>
struct unary_applied_force_functor {
    unary_applied_force_functor(
        size_t part_i_index,
        size_t part_j_index,
        real_t mass,
        field_value_t field_zero,
        f_of_t applied_force_magnitude
    )
        : part_i_index{part_i_index}
        , part_j_index{part_j_index}
        , mass{mass}
        , field_zero{field_zero}
        , applied_force_magnitude_(applied_force_magnitude)
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
            f = dir * applied_force_magnitude_(t);
        } else if (i == part_j_index) [[unlikely]] {
            field_value_t dir = (x[part_j_index] - x[part_i_index]).normalized();
            f = dir * applied_force_magnitude_(t);
        }

        return std::make_pair(f/mass, field_zero);
    }

private:
    const size_t part_i_index, part_j_index;
    const real_t mass;
    const field_value_t field_zero;
    f_of_t applied_force_magnitude_;
};

#endif //APPLIED_FORCE_H
