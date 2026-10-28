#pragma once
#include "../deps/tinyxml2/tinyxml2.h"
#include "reader.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Physical input geometry. All internal coordinates and distances are SI.
struct CapillarySetup {
    using Vec = Eigen::Vector3d;
    struct Drop {
        Vec center;
        double radius;
    };
    std::vector<Vec> positions, velocities;
    std::vector<Drop> drops;
    std::string liquidMode;
    double alpha = 0.1, threshold = 0.5;
    Vec lower, upper;
    std::array<unsigned, 3> cells{};
    double radius = 0, dx = 0;
    unsigned samples = 4, meshFrequency = 1000, particleFrequency = 1000, vtkFrequency = 1000,
             demFrequency = 1000, logFrequency = 500;
    static double number(const tinyxml2::XMLElement *e, const char *key) {
        double v;
        if (!e || e->QueryDoubleAttribute(key, &v) != tinyxml2::XML_SUCCESS || !std::isfinite(v))
            throw std::runtime_error(std::string("Missing or invalid XML attribute: ") + key);
        return v;
    }
    static unsigned frequency(const tinyxml2::XMLElement *e, const char *key) {
        double v = number(e, key);
        if (v < 0 || v > 1000000000 || std::floor(v) != v)
            throw std::runtime_error(std::string("Invalid frequency: ") + key);
        return static_cast<unsigned>(v);
    }
    static Vec vector(const tinyxml2::XMLElement *e) {
        return {number(e, "x"), number(e, "y"), number(e, "z")};
    }
    static std::string attribute(const tinyxml2::XMLElement *e, const char *key) {
        if (!e || !e->Attribute(key))
            throw std::runtime_error(std::string("Missing XML attribute: ") + key);
        return e->Attribute(key);
    }
    static double unit(const std::string &s) {
        if (s == "m")
            return 1;
        if (s == "nm")
            return 1e-9;
        if (s == "um")
            return 1e-6;
        throw std::runtime_error("length_unit must be m, nm, or um");
    }
    static CapillarySetup load(const std::filesystem::path &xml, double r, unsigned resolution) {
        CapillarySetup c;
        c.radius = r;
        if (!(std::isfinite(r) && r > 0) || resolution < 4)
            throw std::runtime_error(
                "Positive radius and at least 4 cells per particle diameter required");
        c.dx = 2 * r / resolution;
        tinyxml2::XMLDocument d;
        if (d.LoadFile(xml.string().c_str()) != tinyxml2::XML_SUCCESS)
            throw std::runtime_error("Cannot read setup XML");
        auto *root = d.FirstChildElement("simulation");
        auto *init = root ? root->FirstChildElement("initialization") : nullptr;
        auto *aggregate = init ? init->FirstChildElement("aggregate") : nullptr;
        if (attribute(aggregate, "format") != "vtk")
            throw std::runtime_error("Aggregate format must be vtk");
        if (aggregate->Attribute("length_unit"))
            throw std::runtime_error(
                "VTK aggregate coordinates are in particle-radius units; remove length_unit");
        const auto path = xml.parent_path() / attribute(aggregate, "file");
        // The existing soot-dem loader reads POINTS only and multiplies by r_part.
        // Preflight its supported legacy ASCII layout before invoking it.
        std::ifstream in(path);
        std::string line;
        if (!in || !std::getline(in, line) || line.rfind("# vtk DataFile Version", 0) != 0)
            throw std::runtime_error("Invalid legacy VTK aggregate header: " + path.string());
        if (!std::getline(in, line))
            throw std::runtime_error("Missing VTK title");
        if (!std::getline(in, line) || line.find("ASCII") == std::string::npos)
            throw std::runtime_error("VTK aggregate must be legacy ASCII");
        if (!std::getline(in, line) || line.find("DATASET POLYDATA") == std::string::npos)
            throw std::runtime_error("VTK aggregate must use DATASET POLYDATA");
        std::string token, type;
        long long pointCount = 0;
        if (!(in >> token >> pointCount >> type) || token != "POINTS" || pointCount <= 0 ||
            pointCount > 10000000)
            throw std::runtime_error("Invalid VTK POINTS pointCount");
        for (long long i = 0; i < pointCount; ++i)
            for (int axis = 0; axis < 3; ++axis) {
                double value;
                if (!(in >> value) || !std::isfinite(value) || !std::isfinite(value * r))
                    throw std::runtime_error("Truncated or nonfinite VTK POINTS data");
            }
        c.positions = load_vtk_aggregate(path, r);
        if (c.positions.size() != static_cast<std::size_t>(pointCount))
            throw std::runtime_error("VTK loader returned an unexpected particle pointCount");
        c.velocities.assign(c.positions.size(), Vec::Zero());
        if (c.positions.empty())
            throw std::runtime_error("Aggregate contains no particles");
        auto *liquid = init->FirstChildElement("liquid");
        if (liquid && liquid->NextSiblingElement("liquid"))
            throw std::runtime_error("Specify exactly one liquid element");
        c.liquidMode = attribute(liquid, "mode");
        if (c.liquidMode == "droplets") {
            const double s = unit(attribute(liquid, "length_unit"));
            for (auto *drop = liquid->FirstChildElement("droplet"); drop;
                 drop = drop->NextSiblingElement("droplet")) {
                double radius = s * number(drop, "radius");
                if (!(radius > 0))
                    throw std::runtime_error("Droplet radius must be positive");
                c.drops.push_back({s * vector(drop), radius});
            }
            if (c.drops.empty())
                throw std::runtime_error("No droplets specified");
        } else if (c.liquidMode == "levelset") {
            c.alpha = number(liquid, "alpha");
            c.threshold = number(liquid, "threshold");
            if (c.alpha <= 0 || c.threshold <= 0 || c.threshold >= double(c.positions.size()))
                throw std::runtime_error("Require alpha>0 and 0<threshold<particle count");
        } else
            throw std::runtime_error("liquid mode must be droplets or levelset");
        auto *lattice = root->FirstChildElement("lattice");
        c.samples = frequency(lattice, "initialization_samples_per_axis");
        if (c.samples < 1 || c.samples > 16)
            throw std::runtime_error("Initialization samples per axis must be 1..16");
        auto *output = root->FirstChildElement("output");
        c.meshFrequency = frequency(output, "surface_every_lb_steps");
        c.particleFrequency = frequency(output, "particle_every_lb_steps");
        c.vtkFrequency = frequency(output, "fluid_every_lb_steps");
        c.demFrequency = frequency(output, "dem_every_dem_steps");
        c.logFrequency = frequency(output, "log_every_lb_steps");
        auto *margin = root->FirstChildElement("domain_margin");
        double ms = unit(attribute(margin, "length_unit"));
        Vec m = ms * vector(margin);
        if ((m.array() < 0).any())
            throw std::runtime_error("Domain margins must be nonnegative");
        // Keep at least three lattice cells of boundary clearance for remapping/stencils.
        m = m.cwiseMax(Vec::Constant(3 * c.dx));
        double extent = r;
        if (c.liquidMode == "levelset") {
            // Outside every expanded particle AABB each term is <= threshold/N.
            // erfc is algebraically identical to 1+erf(-z), without tail cancellation.
            double target = c.threshold / double(c.positions.size());
            double lo = 0, hi = 1;
            while (.5 * std::erfc(hi) > target)
                hi *= 2;
            for (int k = 0; k < 80; ++k) {
                double mid = (lo + hi) / 2;
                if (.5 * std::erfc(mid) > target)
                    lo = mid;
                else
                    hi = mid;
            }
            extent = r * std::max(1.0, 1 + c.alpha * hi);
            if (!std::isfinite(extent))
                throw std::runtime_error("Levelset extent overflow");
        }
        c.lower = c.positions[0] - Vec::Constant(extent);
        c.upper = c.positions[0] + Vec::Constant(extent);
        for (const auto &p : c.positions) {
            c.lower = c.lower.cwiseMin(p - Vec::Constant(extent));
            c.upper = c.upper.cwiseMax(p + Vec::Constant(extent));
        }
        for (const auto &drop : c.drops) {
            c.lower = c.lower.cwiseMin(drop.center - Vec::Constant(drop.radius));
            c.upper = c.upper.cwiseMax(drop.center + Vec::Constant(drop.radius));
        }
        c.lower -= m;
        c.upper += m;
        for (int a = 0; a < 3; ++a) {
            c.lower[a] = std::floor(c.lower[a] / c.dx) * c.dx;
            double n = std::ceil((c.upper[a] - c.lower[a]) / c.dx);
            if (!std::isfinite(n) || n < 7 || n > 100000)
                throw std::runtime_error("Invalid automatic domain size");
            c.cells[a] = static_cast<unsigned>(n);
            c.upper[a] = c.lower[a] + n * c.dx;
        }
        double count = double(c.cells[0]) * c.cells[1] * c.cells[2];
        if (count > 100000000)
            throw std::runtime_error(
                "Automatic domain exceeds 100 million cells; reduce resolution/margins");
        return c;
    }
    bool solid(const Vec &p) const {
        for (const auto &x : positions)
            if ((p - x).squaredNorm() <= radius * radius)
                return true;
        return false;
    }
    double psi(const Vec &p) const {
        double sum = 0;
        for (const auto &x : positions)
            sum += .5 * std::erfc(((p - x).norm() / radius - 1) / alpha);
        return sum;
    }
    bool liquid(const Vec &p) const {
        bool inside = false;
        if (liquidMode == "levelset")
            inside = psi(p) >= threshold;
        else
            for (const auto &d : drops)
                if ((p - d.center).squaredNorm() <= d.radius * d.radius) {
                    inside = true;
                    break;
                }
        return inside && !solid(p);
    }
    double fill(const Vec &center) const {
        unsigned wet = 0;
        for (unsigned z = 0; z < samples; ++z)
            for (unsigned y = 0; y < samples; ++y)
                for (unsigned x = 0; x < samples; ++x)
                    if (liquid(center + dx * Vec((x + .5) / samples - .5, (y + .5) / samples - .5,
                                                 (z + .5) / samples - .5)))
                        ++wet;
        return double(wet) / (samples * samples * samples);
    }
};
