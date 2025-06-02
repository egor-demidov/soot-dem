/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

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

template <typename field_value_t, typename real_t>
struct sphere_coating_functor {
    sphere_coating_functor(real_t magnitude, real_t mass, field_value_t field_zero, real_t particle_radius, std::vector<field_value_t> const & initial_position, real_t total_time)
        : magnitude{magnitude}, mass{mass}, r_part{particle_radius}, t_tot{total_time}, field_zero{field_zero}, initial_position(initial_position)
    {
        // Compute center of mass
        field_value_t sum = field_value_t::Zero();
        for (const auto& xi : initial_position) {
            sum += xi;
        }
        com_current = sum / static_cast<real_t>(initial_position.size());

        // Compute initial radius (max distance from COM)
        initial_radius = 0.0;
        for (const auto& xi : initial_position) {
            real_t dist = (xi - com_current).norm();
            if (dist > initial_radius) {
                initial_radius = dist;
            }
        }

        // Compute shrinking velocity
        shrink_rate = initial_radius / t_tot;
    }

    std::pair<field_value_t, field_value_t> operator () (
        size_t i, size_t j,
        std::vector<field_value_t> const & x,
        std::vector<field_value_t> const & v [[maybe_unused]],
        std::vector<field_value_t> const & theta [[maybe_unused]],
        std::vector<field_value_t> const & omega [[maybe_unused]],
        real_t t) const
    {
        // Compute current sphere radius and make sure it doesn't go below 0
        real_t current_radius = std::max(initial_radius - shrink_rate * current_time, 0.0);

        // Check if both particles are on the shrinking sphere's surface
        real_t dist_i = (x[i] - com_current).norm();
        real_t dist_j = (x[j] - com_current).norm();

        // tolerance = particle radius
        real_t tolerance = r_part;

        bool i_on_surface = std::abs(dist_i - current_radius) < tolerance;
        bool j_on_surface = std::abs(dist_j - current_radius) < tolerance;

        if (i_on_surface && j_on_surface) {
            field_value_t distance = x[j] - x[i];
            real_t distance_magnitude = distance.norm();

            field_value_t n = distance / distance_magnitude;
            field_value_t f = magnitude * n;

            return std::make_pair(f / mass, field_zero);
        }

        return std::make_pair(field_zero, field_zero);
    }

    void set_time(real_t t) const {
        current_time = t;
    }

    // update center of mass
    void updateCOM(std::vector<field_value_t> particlePos) {
        field_value_t sum = field_value_t::Zero();
        for (const auto& xi : particlePos) {
            sum += xi;
        }
        com_current = sum / static_cast<real_t>(particlePos.size());
    }

    real_t get_radius() {
        return std::max(initial_radius - shrink_rate * current_time, 0.0);
    }

    field_value_t get_COM() {
        return com_current;
    }

private:
    const real_t magnitude, mass, r_part, t_tot;
    const field_value_t field_zero;
    real_t initial_radius;
    real_t shrink_rate;
    std::vector<field_value_t> initial_position;
    mutable real_t current_time = 0.0;
    mutable field_value_t com_current;
};

#endif //SOOT_AFM_COATING_FORCE_H
