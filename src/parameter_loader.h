/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_PARAMETER_LOADER_H
#define SOOT_AFM_PARAMETER_LOADER_H

#include <map>
#include <string>
#include <iostream>
#include <filesystem>

struct parameter_store_t {
    std::string simulation_type; // Used to validate that the correct input file is being used
    std::map<std::string, double> reals;
    std::map<std::string, long> integers;
    std::map<std::string, std::string> strings;
    std::map<std::string, std::filesystem::path> paths;
};

std::ostream & operator << (std::ostream & os, parameter_store_t const & parameter_store);

double get_real_parameter(parameter_store_t const & parameter_store, std::string const & id);

long get_integer_parameter(parameter_store_t const & parameter_store, std::string const & id);

std::string const & get_string_parameter(parameter_store_t const & parameter_store, std::string const & id);

std::filesystem::path const & get_path_parameter(parameter_store_t const & parameter_store, std::string const & id);

parameter_store_t load_parameters(std::filesystem::path const & parameter_file_path);

bool has_string_parameter(parameter_store_t const & parameter_store, std::string const & id);

bool has_integer_parameter(parameter_store_t const & parameter_store, std::string const & id);

bool has_real_parameter(parameter_store_t const & parameter_store, std::string const & id);

bool has_path_parameter(parameter_store_t const & parameter_store, std::string const & id);

#endif //SOOT_AFM_PARAMETER_LOADER_H
