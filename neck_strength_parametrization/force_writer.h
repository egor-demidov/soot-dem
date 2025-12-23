//
// Created by mail on 12/23/2025.
//

#ifndef FORCE_WRITER_H
#define FORCE_WRITER_H

#include <filesystem>
#include <fstream>
#include <vector>

#include <Eigen/Eigen>

class ForceWriter {
public:
    explicit ForceWriter(std::filesystem::path dump_file_name);

    void add_line( long step,
                            Eigen::Vector3d const & dr,
                            double xi_n,
                            Eigen::Vector3d const & xi_t,
                            Eigen::Vector3d const & xi_r,
                            Eigen::Vector3d const & xi_o,
                            Eigen::Vector3d const & f,
                            Eigen::Vector3d const & tau);

private:
    std::ofstream ofs;
    const std::filesystem::path dump_file_name;
};

#endif //FORCE_WRITER_H
