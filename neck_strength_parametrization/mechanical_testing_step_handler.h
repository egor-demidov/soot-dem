//
// Created by mail on 12/23/2025.
//

#ifndef MECHANICAL_TESTING_STEP_HANDLER_H
#define MECHANICAL_TESTING_STEP_HANDLER_H

#include <libgran/granular_system/granular_system_neighbor_list.h>

// A custom step handler that prevents velocity updates of primary particles
template <typename field_container_t, typename field_value_t>
struct mechanical_testing_step_handler : public rotational_step_handler<field_container_t, field_value_t> {
public:
    // This method increments the specified value in the x buffer
    void increment_v(long n,
                     field_value_t const & dx [[maybe_unused]],
                     typename field_container_t::iterator x_begin_itr [[maybe_unused]],
                     typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],
                     typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],
                     typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],
                     typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],
                     typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) const {}
};

#endif //MECHANICAL_TESTING_STEP_HANDLER_H
