/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include <iostream>

#include <tinyxml2/tinyxml2.h>

#include "parameter_loader.h"

std::ostream & operator << (std::ostream & os, parameter_store_t const & parameter_store) {
    os << "REAL PARAMETERS:\n";
    for (const auto & real : parameter_store.reals) {
        os << "\t" << real.first << ": " << real.second << "\n";
    }
    os << "INTEGER PARAMETERS:\n";
    for (const auto & integer : parameter_store.integers) {
        os << "\t" << integer.first << ": " << integer.second << "\n";
    }
    os << "STRING PARAMETERS:\n";
    for (const auto & string : parameter_store.strings) {
        os << "\t" << string.first << ": " << string.second << "\n";
    }
    os << "PATH PARAMETERS:\n";
    for (const auto & path : parameter_store.paths) {
        os << "\t" << path.first << ": " << path.second << "\n";
    }
    return os;
}

bool has_string_parameter(parameter_store_t const & parameter_store, std::string const & id) {
    return parameter_store.strings.find(id) != parameter_store.strings.end();
}

bool has_integer_parameter(parameter_store_t const & parameter_store, std::string const & id) {
    return parameter_store.integers.find(id) != parameter_store.integers.end();
}

bool has_real_parameter(parameter_store_t const & parameter_store, std::string const & id) {
    return parameter_store.reals.find(id) != parameter_store.reals.end();
}

bool has_path_parameter(parameter_store_t const & parameter_store, std::string const & id) {
    return parameter_store.paths.find(id) != parameter_store.paths.end();
}

double get_real_parameter(parameter_store_t const & parameter_store, std::string const & id) {
    auto parameter_itr = parameter_store.reals.find(id);

    if (parameter_itr == parameter_store.reals.end()) {
        std::cerr << "Real parameter `" << id << "` is required but was not specified in the parameter file" << std::endl;
        exit(EXIT_FAILURE);
    }

    return parameter_itr->second;
}

long get_integer_parameter(parameter_store_t const & parameter_store, std::string const & id) {
    auto parameter_itr = parameter_store.integers.find(id);

    if (parameter_itr == parameter_store.integers.end()) {
        std::cerr << "Integer parameter `" << id << "` is required but was not specified in the parameter file" << std::endl;
        exit(EXIT_FAILURE);
    }

    return parameter_itr->second;
}

std::string const & get_string_parameter(parameter_store_t const & parameter_store, std::string const & id) {
    auto parameter_itr = parameter_store.strings.find(id);

    if (parameter_itr == parameter_store.strings.end()) {
        std::cerr << "String parameter `" << id << "` is required but was not specified in the parameter file" << std::endl;
        exit(EXIT_FAILURE);
    }

    return parameter_itr->second;
}

std::filesystem::path const & get_path_parameter(parameter_store_t const & parameter_store, std::string const & id) {
    auto parameter_itr = parameter_store.paths.find(id);

    if (parameter_itr == parameter_store.paths.end()) {
        std::cerr << "Path parameter `" << id << "` is required but was not specified in the parameter file" << std::endl;
        exit(EXIT_FAILURE);
    }

    return parameter_itr->second;
}

void set_real(parameter_store_t & parameter_store, std::string const & id, double value) {
    auto parameter_itr = parameter_store.reals.find(id);

    if (parameter_itr != parameter_store.reals.end()) {
        std::cerr << "Multiple definition of real parameter `" << id << "`" << std::endl;
        exit(EXIT_FAILURE);
    }

    parameter_store.reals[id] = value;
}

void set_integer(parameter_store_t & parameter_store, std::string const & id, long value) {
    auto parameter_itr = parameter_store.integers.find(id);

    if (parameter_itr != parameter_store.integers.end()) {
        std::cerr << "Multiple definition of integer parameter `" << id << "`" << std::endl;
        exit(EXIT_FAILURE);
    }

    parameter_store.integers[id] = value;
}

void set_string(parameter_store_t & parameter_store, std::string const & id, std::string const & value) {
    auto parameter_itr = parameter_store.strings.find(id);

    if (parameter_itr != parameter_store.strings.end()) {
        std::cerr << "Multiple definition of string parameter `" << id << "`" << std::endl;
        exit(EXIT_FAILURE);
    }

    parameter_store.strings[id] = value;
}

void set_path(parameter_store_t & parameter_store, std::string const & id, std::filesystem::path const & value) {
    auto parameter_itr = parameter_store.paths.find(id);

    if (parameter_itr != parameter_store.paths.end()) {
        std::cerr << "Multiple definition of path parameter `" << id << "`" << std::endl;
        exit(EXIT_FAILURE);
    }

    parameter_store.paths[id] = value;
}

void parse_parameters(parameter_store_t & parameter_store, tinyxml2::XMLElement * root_element, std::filesystem::path const & reference_path);

void load_include_file(parameter_store_t & parameter_store, std::filesystem::path const & include_file_path,
    std::filesystem::path const & reference_path) {

    auto effective_path = reference_path / include_file_path;

    tinyxml2::XMLDocument doc;
    doc.LoadFile(effective_path.string().c_str());

    if (doc.Error()) {
        std::cerr << "Unable to read the parameter file at " << effective_path << std::endl;
        exit(EXIT_FAILURE);
    }

    auto root = doc.FirstChildElement("include_file");

    if (root == nullptr) {
        std::cerr << "The root element of the include file must be `include_file`" << std::endl;
        exit(EXIT_FAILURE);
    }

    parse_parameters(parameter_store, root, effective_path.parent_path());
}

void parse_parameters(parameter_store_t & parameter_store, tinyxml2::XMLElement * root_element, std::filesystem::path const & reference_path) {
    // Iterate over parameter declarations
    for (auto element = root_element->FirstChildElement("let"); element != nullptr; element = element->NextSiblingElement("let")) {
        auto id = element->FindAttribute("id");

        if (id == nullptr) {
            std::cerr << "Every `let` element must contain an `id` attribute" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string id_string(id->Value());

        auto type = element->FindAttribute("type");

        if (type == nullptr) {
            std::cerr << "Every `let` element must contain a `type` attribute" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string type_string(type->Value());
        auto text = element->GetText();
        if (type_string == "real") {
            set_real(parameter_store, id_string, std::stod(text));
        } else if (type_string == "integer") {
            set_integer(parameter_store, id_string, std::stol(text));
        } else if (type_string == "string") {
            set_string(parameter_store, id_string, text);
        } else if (type_string == "path") {
            set_path(parameter_store, id_string, reference_path / std::filesystem::path(text));
        } else {
            std::cerr << "Unrecognized parameter type `" << type_string << "`" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    for (auto element = root_element->FirstChildElement("include"); element != nullptr; element = element->NextSiblingElement("include")) {
        auto text = element->GetText();
        load_include_file(parameter_store, text, reference_path);
    }
}

parameter_store_t load_parameters(std::filesystem::path const & parameter_file_path) {
    parameter_store_t parameter_store;

    tinyxml2::XMLDocument doc;
    doc.LoadFile(parameter_file_path.string().c_str());

    if (doc.Error()) {
        std::cerr << "Unable to read the parameter file at " << parameter_file_path << std::endl;
        exit(EXIT_FAILURE);
    }

    auto root = doc.FirstChildElement("simulation");

    if (root == nullptr) {
        std::cerr << "The root element of the parameter file must be `simulation`" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto simulation_type = root->FindAttribute("type");

    if (simulation_type == nullptr) {
        std::cerr << "The `simulation` element must contain a `type` attribute" << std::endl;
        exit(EXIT_FAILURE);
    }

    parameter_store.simulation_type = std::string(simulation_type->Value());

    parse_parameters(parameter_store, root, parameter_file_path.parent_path());

    return parameter_store;
}
