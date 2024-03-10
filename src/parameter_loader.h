//
// Created by egor on 3/6/24.
//

#ifndef SOOT_AFM_PARAMETER_LOADER_H
#define SOOT_AFM_PARAMETER_LOADER_H

#include <map>
#include <string>
#include <iostream>

struct parameter_store_t {
    std::string simulation_type; // Used to validate that the correct input file is being used
    std::map<std::string, double> reals;
    std::map<std::string, long> integers;
    std::map<std::string, std::string> strings;
};

std::ostream & operator << (std::ostream & os, parameter_store_t const & parameter_store);

double get_real_parameter(parameter_store_t const & parameter_store, std::string const & id);

long get_integer_parameter(parameter_store_t const & parameter_store, std::string const & id);

std::string const & get_string_parameter(parameter_store_t const & parameter_store, std::string const & id);

parameter_store_t load_parameters(std::string const & parameter_file_path);

#endif //SOOT_AFM_PARAMETER_LOADER_H
