/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_FIXED_STEP_HANDLER_H
#define SOOT_AFM_FIXED_STEP_HANDLER_H

#include <vector>

template <typename field_container_t,
        typename field_value_t,
        template <
        typename _field_container_t,
        typename _field_value_t>
        typename base_step_handler_t>
struct fixed_step_handler {
    fixed_step_handler(std::vector<bool> fixed_particles, base_step_handler_t<field_container_t, field_value_t> base_step_handler) :
            fixed_particles(std::move(fixed_particles)), base_step_handler(std::move(base_step_handler)) {}

    // This method increments the specified value in the x buffer
    void increment_x(size_t n,                                                                              // index of the value to increment
                     field_value_t const & dx,                                                              // value of the position increment
                     typename field_container_t::iterator x_begin_itr,                                      // iterator pointing to the start of the x buffer
                     typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],               // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                     typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],           // iterator pointing to the start of the theta buffer
                     typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],           // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) const {   // iterator pointing to the start of the alpha buffer

        if (!fixed_particles[n])
            base_step_handler.increment_x(n, dx, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

    // This method increments the specified value in the v buffer
    void increment_v(size_t n,                                                                              // index of the value to increment
                     field_value_t const & dv,                                                              // value of the velocity increment
                     typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],               // iterator pointing to the start of the x buffer
                     typename field_container_t::iterator v_begin_itr,                                      // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                     typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],           // iterator pointing to the start of the theta buffer
                     typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],           // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) const {   // iterator pointing to the start of the alpha buffer

        if (!fixed_particles[n])
            base_step_handler.increment_v(n, dv, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

    // This method increments the specified value in the theta buffer
    void increment_theta(size_t n,                                                                              // index of the value to increment
                         field_value_t const & dtheta,                                                          // value of the angle increment
                         typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],               // iterator pointing to the start of the x buffer
                         typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],               // iterator pointing to the start of the v buffer
                         typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                         typename field_container_t::iterator theta_begin_itr,                                  // iterator pointing to the start of the theta buffer
                         typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],           // iterator pointing to the start of the omega buffer
                         typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) const {   // iterator pointing to the start of the alpha buffer

        base_step_handler.increment_theta(n, dtheta, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

    // This method increments the specified value in the omega buffer
    void increment_omega(size_t n,                                                                              // index of the value to increment
                         field_value_t const & domega,                                                          // value of the angular velocity increment
                         typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],               // iterator pointing to the start of the x buffer
                         typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],               // iterator pointing to the start of the v buffer
                         typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                         typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],           // iterator pointing to the start of the theta buffer
                         typename field_container_t::iterator omega_begin_itr,                                  // iterator pointing to the start of the omega buffer
                         typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) const {   // iterator pointing to the start of the alpha buffer

        base_step_handler.increment_omega(n, domega, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

private:
    const std::vector<bool> fixed_particles;
    base_step_handler_t<field_container_t, field_value_t> base_step_handler;
};

#endif //SOOT_AFM_FIXED_STEP_HANDLER_H
