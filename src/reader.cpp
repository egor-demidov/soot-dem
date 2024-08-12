/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include <iostream>
#include <fstream>
#include <sstream>

#include <tinyxml2/tinyxml2.h>

#include "reader.h"

std::vector<Eigen::Vector3d> load_flage_aggregate(std::filesystem::path const & path, double r_part) {
    tinyxml2::XMLDocument doc;
    doc.LoadFile(path.string().c_str());

    if (doc.Error()) {
        std::cerr << "Unable to open file " << path << std::endl;
        return {};
    }

    auto root = doc.FirstChildElement("geometry_complex");
    auto aggregate = root->FirstChildElement("data_particle");

    std::vector<Eigen::Vector3d> out;

    // Iterate over monomers
    for (auto monomer = aggregate->FirstChildElement("particle"); monomer != nullptr; monomer = monomer->NextSiblingElement("particle")) {
        auto position_string = monomer->FirstChildElement("position")->GetText();
        std::stringstream ss(position_string);
        double x, y, z;
        std::string _;
        ss >> x >> _ >> y >> _ >> z;
        out.emplace_back(Eigen::Vector3d{x, y, z} * r_part);
    }

    return out;
}

std::vector<Eigen::Vector3d> load_mackowski_aggregate(std::filesystem::path const & path, double r_part) {
    std::ifstream ifs(path);

    if (!ifs.good()) {
        std::cerr << "Unable to read the aggregate file: " << path << std::endl;
        return {};
    }

    std::string line;
    std::vector<Eigen::Vector3d> x0;

    while (getline(ifs, line)) {
        if (!line.empty()) {
            std::istringstream oss(line);
            double _, x, y, z;
            oss >> _ >> x >> y >> z;
            x0.emplace_back(x * r_part, y * r_part, z * r_part);
        }
    }

    return x0;
}

std::vector<Eigen::Vector3d> load_vtk_aggregate(std::filesystem::path const & path, double r_part) {
    std::ifstream ifs(path);

    if (!ifs.good()) {
        std::cerr << "Unable to read the aggregate file: " << path << std::endl;
        return {};
    }

    std::string _;
    size_t num_pts;

    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    ifs >> _ >> num_pts >> _;

    std::vector<Eigen::Vector3d> out;
    out.resize(num_pts);

    double x, y, z;
    for (size_t i = 0; i < num_pts; i ++) {
        ifs >> x >> y >> z;
        out[i] = Eigen::Vector3d {x, y, z} * r_part;
    }

    return out;
}

std::vector<bool> load_necks(std::filesystem::path const & path, size_t n_part) {
    std::ifstream ifs(path);

    if (!ifs.good()) {
        std::cerr << "Unable to read the aggregate file: " << path << std::endl;
        return {};
    }

    std::string _;
    size_t num_pts;

    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    ifs >> _ >> num_pts >> _;
    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);
    std::getline(ifs, _);

    std::vector<bool> bonded_contacts(n_part*n_part, false);
    for (size_t i = 0; i < num_pts * 2; i ++) {
        int index;
        ifs >> index;
        bonded_contacts[index] = true;
    }

    return bonded_contacts;
}
