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
    sphere_coating_functor(real_t magnitude, real_t mass, field_value_t field_zero, real_t particle_radius, std::vector<field_value_t> const & initial_position, real_t drag_coefficient, real_t total_time)
        : magnitude{magnitude}, mass{mass}, r_part{particle_radius}, t_tot{total_time}, drag_coefficient(drag_coefficient), field_zero{field_zero}, initial_position(initial_position)
    {
        // Compute initial center of mass
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

        //move sphere

        //Eigen::Vector3d vec;
        //vec << -0.1, 1.0, 0.0;
        //vec *= initial_radius*1.45;

        //com_current += vec;

        // Compute shrinking velocity
        shrink_rate = initial_radius / t_tot;
    }

    std::pair<field_value_t, field_value_t> operator () (
        size_t i,
        std::vector<field_value_t> const & x,
        std::vector<field_value_t> const & v [[maybe_unused]],
        std::vector<field_value_t> const & theta [[maybe_unused]],
        std::vector<field_value_t> const & omega [[maybe_unused]],
        real_t t) const
    {
        // Compute current sphere radius and make sure it doesn't go below 0
        //real_t current_radius = std::max(initial_radius - shrink_rate * current_time, 0.0);

        // Check if particle is on the shrinking sphere's surface
        if (interface_particles_active[i]) {

            field_value_t total_force = field_zero;

            //drag force
            field_value_t difference = com_current - x[i];
            field_value_t r = difference / difference.norm();
            field_value_t relative_velocity = shrink_rate * r - (v[i].dot(r))*r;

            total_force += relative_velocity * drag_coefficient;

            // find force acting on particle from all other particles in sphere
            // for (long j : interface_particles) {
            //     if (j != i) {
            //         field_value_t distance = x[j] - x[i];
            //         real_t distance_magnitude = distance.norm();

            //         field_value_t n = distance / distance_magnitude;
            //         field_value_t f = magnitude * n;

            //         total_force += f;
            //     }
            // }

            return std::make_pair(total_force / mass, field_zero);
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

    // check which particles are on the sphere and add them to the vector; also activate them in the bool vector
    void update_interface_particles(std::vector<field_value_t> const & particle_positions) {
        interface_particles.clear();
        interface_particles_active.assign(particle_positions.size(), false);
        real_t current_radius = get_radius();

        for (int i = 0; i < particle_positions.size(); i++) {
            if (abs((particle_positions[i] - com_current).norm() - current_radius) < r_part) {
                interface_particles.push_back(i);
                interface_particles_active[i] = true;
            }
        }
    }

    real_t get_radius() {
        return std::max(initial_radius - shrink_rate * current_time, 0.0);
    }

    field_value_t get_COM() {
        return com_current;
    }

    std::vector<long> const & get_interface_particles() const {
        return interface_particles;
    }

private:
    const real_t magnitude, mass, r_part, t_tot, drag_coefficient;
    const field_value_t field_zero;
    real_t initial_radius;
    real_t shrink_rate;
    std::vector<field_value_t> initial_position;
    std::vector<long> interface_particles;
    std::vector<bool> interface_particles_active;
    mutable real_t current_time = 0.0;
    mutable field_value_t com_current;
};

template <typename field_value_t, typename real_t>
struct hemisphere_coating_model_t {
    hemisphere_coating_model_t(real_t magnitude, real_t mass, field_value_t field_zero, real_t particle_radius, std::vector<field_value_t> const & initial_position, real_t drag_coefficient, real_t total_time)
        : magnitude{magnitude}, mass{mass}, r_part{particle_radius}, t_tot{total_time}, drag_coefficient(drag_coefficient), field_zero{field_zero}, initial_position(initial_position)
    {
        // Compute initial center of mass
        field_value_t sum = field_value_t::Zero();
        for (const auto& xi : initial_position) {
            sum += xi;
        }
        com_current = sum / static_cast<real_t>(initial_position.size());
        com_current[2] = 0;

        // Compute initial radius (max distance from COM)
        initial_radius = 0.0;
        for (const auto& xi : initial_position) {
            real_t dist = (xi - com_current).norm();
            if (dist > initial_radius) {
                initial_radius = dist;
            }
        }

        //move sphere

        //Eigen::Vector3d vec;
        //vec << -0.1, 1.0, 0.0;
        //vec *= initial_radius*1.45;

        //com_current += vec;

        // Compute shrinking velocity
        shrink_rate = initial_radius / t_tot;
    }

    std::pair<field_value_t, field_value_t> operator () (
        size_t i,
        std::vector<field_value_t> const & x,
        std::vector<field_value_t> const & v [[maybe_unused]],
        std::vector<field_value_t> const & theta [[maybe_unused]],
        std::vector<field_value_t> const & omega [[maybe_unused]],
        real_t t) const
    {
        // Compute current sphere radius and make sure it doesn't go below 0
        //real_t current_radius = std::max(initial_radius - shrink_rate * current_time, 0.0);

        // Check if particle is on the shrinking sphere's surface
        if (interface_particles_active[i]) {

            field_value_t total_force = field_zero;

            field_value_t difference = com_current - x[i];
            field_value_t r = difference / difference.norm();
            field_value_t relative_velocity = shrink_rate * r - (v[i].dot(r))*r;

            total_force += relative_velocity * drag_coefficient;

            // find force acting on particle from all other particles in sphere
            // for (long j : interface_particles) {
            //     if (j != i) {
            //         field_value_t distance = x[j] - x[i];
            //         real_t distance_magnitude = distance.norm();

            //         field_value_t n = distance / distance_magnitude;
            //         field_value_t f = magnitude * n;

            //         total_force += f;
            //     }
            // }

            return std::make_pair(total_force / mass, field_zero);
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
        com_current[2] = 0;
    }

    // check which particles are on the sphere and add them to the vector; also activate them in the bool vector
    void update_interface_particles(std::vector<field_value_t> const & particle_positions) {
        interface_particles.clear();
        interface_particles_active.assign(particle_positions.size(), false);
        real_t current_radius = get_radius();

        for (int i = 0; i < particle_positions.size(); i++) {
            if (abs((particle_positions[i] - com_current).norm() - current_radius) < r_part) {
                interface_particles.push_back(i);
                interface_particles_active[i] = true;
            }
        }
    }

    real_t get_radius() {
        return std::max(initial_radius - shrink_rate * current_time, 0.0);
    }

    field_value_t get_COM() {
        return com_current;
    }

    std::vector<long> const & get_interface_particles() const {
        return interface_particles;
    }

private:
    const real_t magnitude, mass, r_part, t_tot, drag_coefficient;
    const field_value_t field_zero;
    real_t initial_radius;
    real_t shrink_rate;
    std::vector<field_value_t> initial_position;
    std::vector<long> interface_particles;
    std::vector<bool> interface_particles_active;
    mutable real_t current_time = 0.0;
    mutable field_value_t com_current;
};

#endif //SOOT_AFM_COATING_FORCE_H
