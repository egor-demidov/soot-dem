//
// Created by egor on 10/27/26.
//

#ifndef SOOT_AFM_CONTACT_LINE_OUTPUT_H
#define SOOT_AFM_CONTACT_LINE_OUTPUT_H

#include "surface_contact_line.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

inline void writeContactLineVtp(
    const SurfaceLineResult& line,
    const std::filesystem::path& filename
) {
    std::ofstream out(filename);

    if (!out)
        throw std::runtime_error(
            "Cannot create " + filename.string()
        );

    const std::size_t count = line.segments.size();

    out << std::setprecision(17)
        << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"PolyData\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">\n"
        << "<PolyData><Piece NumberOfPoints=\"" << 2 * count
        << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << count
        << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
        << "<Points><DataArray type=\"Float64\" "
        << "NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto& segment : line.segments)
    {
        for (const auto& point : segment)
        {
            out << point.x() << ' '
                << point.y() << ' '
                << point.z() << '\n';
        }
    }

    out << "</DataArray></Points>\n"
        << "<Lines>\n"
        << "<DataArray type=\"Int64\" Name=\"connectivity\" "
        << "format=\"ascii\">\n";

    for (std::size_t i = 0; i < count; ++i)
        out << 2 * i << ' ' << 2 * i + 1 << '\n';

    out << "</DataArray>\n"
        << "<DataArray type=\"Int64\" Name=\"offsets\" "
        << "format=\"ascii\">\n";

    for (std::size_t i = 0; i < count; ++i)
        out << 2 * (i + 1) << '\n';

    out << "</DataArray></Lines>\n"
        << "</Piece></PolyData></VTKFile>\n";

    out.close();

    if (!out)
        throw std::runtime_error(
            "Failed writing " + filename.string()
        );
}

#endif //SOOT_AFM_CONTACT_LINE_OUTPUT_H
