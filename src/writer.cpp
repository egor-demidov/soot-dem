/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include <iostream>
#include <fstream>

#include "writer.h"

// Overloaded dump_necks with energy write-outs
bool dump_necks(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                std::vector<bool> const & bonded_contacts, double r_part, double k_n_bond,
                double k_t_bond, double k_r_bond, double k_o_bond,
                std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> const & contact_springs,
                size_t n_part) {

    std::stringstream out_file_name;
    out_file_name << dir << "/necks_" << count << ".vtk";
    std::ofstream ofs(out_file_name.str());

    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        return false;
    }

    size_t neck_count = std::count(bonded_contacts.begin(), bonded_contacts.end(), true) / 2u;

    ofs << "# vtk DataFile Version 4.0" << "\n";
    ofs << "Generated by libgran" << "\n";
    ofs << "ASCII" << "\n";
    ofs << "DATASET POLYDATA" << "\n";
    ofs << "POINTS " << neck_count << " FLOAT" << "\n";

    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if (!bonded_contacts[i * x.size() + j])
                continue;

            // This is a bonded contact
            Eigen::Vector3d position = (x[j] + x[i]) / 2.0 / r_part; // Compute the position

            ofs << position[0] << " " << position[1] << " " << position[2] << " ";
        }
    }

    ofs << "\n" << "\n";
    ofs << "POINT_DATA " << neck_count << "\n";
    ofs << "FIELD FieldData 5" << "\n";
    ofs << "normals 3 " << neck_count << " double" << "\n";
    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if (!bonded_contacts[i * x.size() + j])
                continue;

            // This is a bonded contact
            Eigen::Vector3d orientation = (x[j] - x[i]).normalized(); // Compute the orientation vector

            ofs << orientation[0] << " " << orientation[1] << " " << orientation[2] << " ";
        }
    }
    ofs << "\n";
    ofs << "connections 2 " << neck_count << " int" << "\n";  // Tuples of 2 to avoid errors with ParaView
    for (size_t i = 0; i < bonded_contacts.size(); i ++) {
        if (bonded_contacts[i])
            ofs << i << " ";
    }
    ofs << "energies 1 " << neck_count << " double" << "\n";
    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if (!bonded_contacts[i * x.size() + j])
                continue;

            double xi_n = (x[i] - x[j]).norm() - 2.0 * r_part;
            auto [xi_t, xi_r, xi_o] = contact_springs[i * n_part + j];

            double total_strain_energy = k_t_bond * xi_t.dot(xi_t)
                                         + k_r_bond * xi_r.dot(xi_r)
                                         + k_o_bond * xi_o.dot(xi_o);

            // Normal component only contributes in case of tension
            if (xi_n > 0.0)
                total_strain_energy += k_n_bond * xi_n * xi_n;

            ofs << total_strain_energy << " ";
        }
    }
    ofs << "energies_normal 1 " << neck_count << " double" << "\n";
    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if (!bonded_contacts[i * x.size() + j])
                continue;

            double xi_n = (x[i] - x[j]).norm() - 2.0 * r_part;

            ofs << k_n_bond * xi_n * xi_n << " ";
        }
    }
    ofs << "energies_tangential_rolling_torsional 1 " << neck_count << " double" << "\n";
    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if (!bonded_contacts[i * x.size() + j])
                continue;

            auto [xi_t, xi_r, xi_o] = contact_springs[i * n_part + j];

            ofs << k_t_bond * xi_t.dot(xi_t) + k_r_bond * xi_r.dot(xi_r) + k_o_bond * xi_o.dot(xi_o) << " ";
        }
    }
    ofs << "\n\n";
    return true;
}

bool dump_necks(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                std::vector<bool> const & bonded_contacts, double r_part) {

    std::stringstream out_file_name;
    out_file_name << dir << "/necks_" << count << ".vtk";
    std::ofstream ofs(out_file_name.str());

    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        return false;
    }

    size_t neck_count = std::count(bonded_contacts.begin(), bonded_contacts.end(), true) / 2u;

    ofs << "# vtk DataFile Version 4.0" << "\n";
    ofs << "Generated by libgran" << "\n";
    ofs << "ASCII" << "\n";
    ofs << "DATASET POLYDATA" << "\n";
    ofs << "POINTS " << neck_count << " FLOAT" << "\n";

    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if (!bonded_contacts[i * x.size() + j])
                continue;

            // This is a bonded contact
            Eigen::Vector3d position = (x[j] + x[i]) / 2.0 / r_part; // Compute the position

            ofs << position[0] << " " << position[1] << " " << position[2] << " ";
        }
    }

    ofs << "\n" << "\n";
    ofs << "POINT_DATA " << neck_count << "\n";
    ofs << "FIELD FieldData 2" << "\n";
    ofs << "normals 3 " << neck_count << " double" << "\n";
    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if (!bonded_contacts[i * x.size() + j])
                continue;

            // This is a bonded contact
            Eigen::Vector3d orientation = (x[j] - x[i]).normalized(); // Compute the orientation vector

            ofs << orientation[0] << " " << orientation[1] << " " << orientation[2] << " ";
        }
    }
    ofs << "\n";
    ofs << "connections 2 " << neck_count << " int" << "\n";  // Tuples of 2 to avoid errors with ParaView
    for (size_t i = 0; i < bonded_contacts.size(); i ++) {
        if (bonded_contacts[i])
            ofs << i << " ";
    }
    ofs << "\n\n";
    return true;
}

// Overload to write neck strengths
bool dump_necks(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                std::vector<bool> const & bonded_contacts, double r_part,
                std::vector<double> const & neck_strengths) {

    std::stringstream out_file_name;
    out_file_name << dir << "/necks_" << count << ".vtk";
    std::ofstream ofs(out_file_name.str());

    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        return false;
    }

    size_t neck_count = std::count(bonded_contacts.begin(), bonded_contacts.end(), true) / 2u;

    ofs << "# vtk DataFile Version 4.0" << "\n";
    ofs << "Generated by libgran" << "\n";
    ofs << "ASCII" << "\n";
    ofs << "DATASET POLYDATA" << "\n";
    ofs << "POINTS " << neck_count << " FLOAT" << "\n";

    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if (!bonded_contacts[i * x.size() + j])
                continue;

            // This is a bonded contact
            Eigen::Vector3d position = (x[j] + x[i]) / 2.0 / r_part; // Compute the position

            ofs << position[0] << " " << position[1] << " " << position[2] << " ";
        }
    }

    ofs << "\n" << "\n";
    ofs << "POINT_DATA " << neck_count << "\n";
    ofs << "FIELD FieldData 3" << "\n";
    ofs << "normals 3 " << neck_count << " double" << "\n";
    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if (!bonded_contacts[i * x.size() + j])
                continue;

            // This is a bonded contact
            Eigen::Vector3d orientation = (x[j] - x[i]).normalized(); // Compute the orientation vector

            ofs << orientation[0] << " " << orientation[1] << " " << orientation[2] << " ";
        }
    }
    ofs << "\n";
    ofs << "connections 2 " << neck_count << " int" << "\n";  // Tuples of 2 to avoid errors with ParaView
    for (size_t i = 0; i < bonded_contacts.size(); i ++) {
        if (bonded_contacts[i])
            ofs << i << " ";
    }
    ofs << "\n";

    ofs << "strengths 1 " << neck_count << " double" << "\n";
    for (int i = 0; i < neck_strengths.size(); i++) {
        ofs << neck_strengths[i] << " ";
    }
    ofs << "\n\n";
    return true;
}

// Overload to write particle positions without velocities and orientations
bool dump_particles(std::string const & name, std::vector<Eigen::Vector3d> const & x, double r_part) {
    std::stringstream out_file_name;
    out_file_name << name << ".vtk";
    std::ofstream ofs(out_file_name.str());

    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        return false;
    }

    ofs << "# vtk DataFile Version 4.0" << "\n";
    ofs << "Generated by libFractalCommon" << "\n";
    ofs << "ASCII" << "\n";
    ofs << "DATASET POLYDATA" << "\n";
    ofs << "POINTS " << x.size() << " FLOAT" << "\n";

    // Coordinates are normalized to avoid issues with ParaView and small particles
    for (auto const & p : x) {
        ofs << p[0] / r_part << " " << p[1] / r_part << " " << p[2] / r_part << " ";
    }
    ofs << "\n";
    return true;
}

bool dump_particles(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                    std::vector<Eigen::Vector3d> const & theta,
                    std::vector<Eigen::Vector3d> const & v,
                    std::vector<Eigen::Vector3d> const & a,
                    std::vector<Eigen::Vector3d> const & omega,
                    std::vector<Eigen::Vector3d> const & alpha,
                    double r_part) {
    std::stringstream out_file_name;
    out_file_name << dir << "/particles_" << count << ".vtk";
    std::ofstream ofs(out_file_name.str());

    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        return false;
    }

    ofs << "# vtk DataFile Version 4.0" << "\n";
    ofs << "Generated by libFractalCommon" << "\n";
    ofs << "ASCII" << "\n";
    ofs << "DATASET POLYDATA" << "\n";
    ofs << "POINTS " << x.size() << " FLOAT" << "\n";

    // Coordinates are normalized to avoid issues with ParaView and small particles
    for (auto const & p : x) {
        ofs << p[0] / r_part << " " << p[1] / r_part << " " << p[2] / r_part << " ";
    }

    ofs << "\n\n";
    ofs << "POINT_DATA " << x.size() << "\n";
    ofs << "FIELD FieldData 7" << "\n";
    ofs << "v 1 " << x.size() << " double\n";
    for (auto const & p : v) {
        ofs << p.norm() << " ";
    }
    ofs << "\n";
    ofs << "a 1 " << x.size() << " double\n";
    for (auto const & p : a) {
        ofs << p.norm() << " ";
    }
    ofs << "\n";
    ofs << "omega 1 " << x.size() << " double\n";
    for (auto const & p : omega) {
        ofs << p.norm() << " ";
    }
    ofs << "\n";
    ofs << "alpha 1 " << x.size() << " double\n";
    for (auto const & p : alpha) {
        ofs << p.norm() << " ";
    }
    ofs << "\n";
    ofs << "thetax 3 " << x.size() << " double\n";
    for (auto const & t : theta) {
        Eigen::Matrix3d m;
        m = Eigen::AngleAxis(t[0], Eigen::Vector3d::UnitX())
            * Eigen::AngleAxis(t[1], Eigen::Vector3d::UnitY())
            * Eigen::AngleAxis(t[2], Eigen::Vector3d::UnitZ());

        auto unit = Eigen::Vector3d::UnitX();
        auto orient = m*unit;

        ofs << orient[0] << " " << orient[1] << " " << orient[2] << " ";
    }
    ofs << "\n" << "thetay 3 " << x.size() << " double" << "\n";
    for (auto const & t : theta) {
        Eigen::Matrix3d m;
        m = Eigen::AngleAxis(t[0], Eigen::Vector3d::UnitX())
            * Eigen::AngleAxis(t[1], Eigen::Vector3d::UnitY())
            * Eigen::AngleAxis(t[2], Eigen::Vector3d::UnitZ());

        auto unit = Eigen::Vector3d::UnitY();
        auto orient = m*unit;

        ofs << orient[0] << " " << orient[1] << " " << orient[2] << " ";
    }
    ofs << "\n" << "thetaz 3 " << x.size() << " double" << "\n";
    for (auto const & t : theta) {
        Eigen::Matrix3d m;
        m = Eigen::AngleAxis(t[0], Eigen::Vector3d::UnitX())
            * Eigen::AngleAxis(t[1], Eigen::Vector3d::UnitY())
            * Eigen::AngleAxis(t[2], Eigen::Vector3d::UnitZ());

        auto unit = Eigen::Vector3d::UnitZ();
        auto orient = m*unit;

        ofs << orient[0] << " " << orient[1] << " " << orient[2] << " ";
    }
    return true;
}

bool dump_particles(std::string const & dir, size_t count, std::vector<Eigen::Vector3d> const & x,
                    std::vector<Eigen::Vector3d> const & v,
                    std::vector<Eigen::Vector3d> const & a,
                    std::vector<Eigen::Vector3d> const & omega,
                    std::vector<Eigen::Vector3d> const & alpha,
                    double r_part) {
    std::stringstream out_file_name;
    out_file_name << dir << "/particles_" << count << ".vtk";
    std::ofstream ofs(out_file_name.str());

    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        return false;
    }

    ofs << "# vtk DataFile Version 4.0" << "\n";
    ofs << "Generated by libFractalCommon" << "\n";
    ofs << "ASCII" << "\n";
    ofs << "DATASET POLYDATA" << "\n";
    ofs << "POINTS " << x.size() << " FLOAT" << "\n";

    // Coordinates are normalized to avoid issues with ParaView and small particles
    for (auto const & p : x) {
        ofs << p[0] / r_part << " " << p[1] / r_part << " " << p[2] / r_part << " ";
    }

    ofs << "\n\n";
    ofs << "POINT_DATA " << x.size() << "\n";
    ofs << "FIELD FieldData 4" << "\n";
    ofs << "v 1 " << x.size() << " double\n";
    for (auto const & p : v) {
        ofs << p.norm() << " ";
    }
    ofs << "\n";
    ofs << "a 1 " << x.size() << " double\n";
    for (auto const & p : a) {
        ofs << p.norm() << " ";
    }
    ofs << "\n";
    ofs << "omega 1 " << x.size() << " double\n";
    for (auto const & p : omega) {
        ofs << p.norm() << " ";
    }
    ofs << "\n";
    ofs << "alpha 1 " << x.size() << " double\n";
    for (auto const & p : alpha) {
        ofs << p.norm() << " ";
    }
    ofs << "\n";
    return true;
}

bool dump_sphere(std::string const & dir, size_t count,
                 Eigen::Vector3d const & sphereCenter,
                 double const sphereRadius, double r_part) {
    std::stringstream out_file_name;
    out_file_name << dir << "/sphere_" << count << ".vtk";
    std::ofstream ofs(out_file_name.str());

    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        return false;
    }

    ofs << "# vtk DataFile Version 4.0" << "\n";
    ofs << "Generated by libFractalCommon" << "\n";
    ofs << "ASCII" << "\n";
    ofs << "DATASET POLYDATA" << "\n";
    ofs << "POINTS " << 1 << " FLOAT" << "\n";

    // Coordinates are normalized to avoid issues with ParaView and small particles
    ofs << sphereCenter[0] / r_part << " " << sphereCenter[1] / r_part << " " << sphereCenter[2] / r_part << "\n";

    ofs << "\n\n";
    ofs << "POINT_DATA " << 1 << "\n";
    ofs << "FIELD FieldData 1" << "\n";
    ofs << "r 1 " << 1 << " double\n";
    // sphere radius is normalized
    ofs << sphereRadius / r_part << " ";
    return true;
}
