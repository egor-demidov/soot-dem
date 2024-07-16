//
// Created by egor on 3/3/24.
//

#ifndef SOOT_AFM_COATING_FORCE_H
#define SOOT_AFM_COATING_FORCE_H

template <typename field_value_t, typename real_t>
struct binary_coating_functor {
    binary_coating_functor(real_t d_cutoff, real_t magnitude, real_t drop_rate, real_t mass, field_value_t field_zero) :
    d_cutoff{d_cutoff}, magnitude{magnitude}, drop_rate{drop_rate}, mass{mass}, field_zero{field_zero} {}

    std::pair<field_value_t, field_value_t> operator () (size_t i, size_t j,
            std::vector<field_value_t> const & x,
             std::vector<field_value_t> const & v [[maybe_unused]],
             std::vector<field_value_t> const & theta [[maybe_unused]],
             std::vector<field_value_t> const & omega [[maybe_unused]],
             real_t t [[maybe_unused]]) const {

        field_value_t distance = x[j] - x[i];
        real_t distance_magnitude = distance.norm();
        field_value_t n = distance / distance_magnitude;

        field_value_t f = magnitude * (0.5 - 0.5 * tanh(drop_rate * (distance_magnitude - d_cutoff))) * n;

        return std::make_pair(f / mass, field_zero);
    }

private:
    const real_t d_cutoff, magnitude, drop_rate, mass;
    const field_value_t field_zero;
};


// Unary coating functor can be used to model capillary force between a static plane and particles
template <typename field_value_t, typename real_t>
struct surface_coating_functor {
    surface_coating_functor(real_t d_cutoff, real_t magnitude, real_t drop_rate, real_t mass, field_value_t field_zero) :
            d_cutoff{d_cutoff}, magnitude{magnitude}, drop_rate{drop_rate}, mass{mass}, field_zero{field_zero} {}

    std::pair<field_value_t, field_value_t> operator () (size_t i, field_value_t const & x_facet, field_value_t const & v_facet [[maybe_unused]],
                                                         std::vector<field_value_t> const & x, std::vector<field_value_t> const & v [[maybe_unused]],
                                                         std::vector<field_value_t> const & theta [[maybe_unused]], std::vector<field_value_t> const & omega [[maybe_unused]], real_t t [[maybe_unused]]) const {

        field_value_t distance = (x_facet - x[i]); // Distance vector
        real_t distance_magnitude = distance.norm();
        field_value_t n = distance / distance_magnitude;

        field_value_t f = magnitude * (0.5 - 0.5 * tanh(drop_rate * (distance_magnitude - d_cutoff))) * n;


        return std::make_pair(f/mass, field_zero);
    }

private:
    const real_t d_cutoff, magnitude, drop_rate, mass;
    const field_value_t field_zero;
};

#endif //SOOT_AFM_COATING_FORCE_H
