//
// Created by mail on 12/23/2025.
//

#include <iostream>

#include "force_writer.h"

ForceWriter::ForceWriter(std::filesystem::path dump_file_name) : dump_file_name(std::move(dump_file_name)), ofs(dump_file_name) {
    if (!ofs.good()) {
        std::cerr << "Unable to create a dump force dump at " << this->dump_file_name << std::endl;
        return;
    }

    // Write header
    ofs << "n, dx, dy, dz, xi_n, xi_tx, xi_ty, xi_tz, xi_rx, xi_ry, xi_rz, xi_ox, xi_oy, xi_oz, fx, fy, fz, taux, tauy, tauz\n";
}

void ForceWriter::add_line( long step,
                            Eigen::Vector3d const & dr,
                            double xi_n,
                            Eigen::Vector3d const & xi_t,
                            Eigen::Vector3d const & xi_r,
                            Eigen::Vector3d const & xi_o,
                            Eigen::Vector3d const & f,
                            Eigen::Vector3d const & tau) {
    ofs << step
        << ", " << dr[0] << ", " << dr[1] << ", " << dr[2]
        << ", " << xi_n
        << ", " << xi_t[0] << ", " << xi_t[1] << ", " << xi_t[2]
        << ", " << xi_r[0] << ", " << xi_r[1] << ", " << xi_r[2]
        << ", " << xi_o[0] << ", " << xi_o[1] << ", " << xi_o[2]
        << ", " << f[0] << ", " << f[1] << ", " << f[2]
        << ", " << tau[0] << ", " << tau[1] << ", " << tau[2] << "\n";
}

