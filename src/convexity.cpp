//
// Created by egor on 6/12/25.
//

#include "convexity.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>

using Point = boost::geometry::model::d2::point_xy<double>;
using Polygon = boost::geometry::model::polygon<Point>;

double compute_convexity() {
    Polygon poly1;
    boost::geometry::append(poly1.outer(), Point(0, 0));
    boost::geometry::append(poly1.outer(), Point(0, 1));
    boost::geometry::append(poly1.outer(), Point(1, 1));
    boost::geometry::append(poly1.outer(), Point(1, 0));
    boost::geometry::append(poly1.outer(), Point(0, 0));

    std::cout << boost::geometry::area(poly1) << std::endl;

    Polygon poly2;

    for (long n = 0; n < 50; n ++) {
        double theta = -2.0 * M_PI * static_cast<double>(n) / static_cast<double>(50-1);
        double x = 1.0 * cos(theta);
        double y = 1.0 * sin(theta);
        boost::geometry::append(poly2.outer(), Point(x, y));
    }

    std::cout << boost::geometry::area(poly2) / M_PI << std::endl;

    return 0.0;
}

static constexpr long N_THETA_PTS = 25l;

double compute_z_convexity(std::vector<Eigen::Vector3d> const & pos, double r_part) {

    std::stack<Polygon> projections_of_monomers;
    for (long i = 0; i < pos.size(); i ++) {

        Polygon poly;
        for (long n = 0; n < N_THETA_PTS; n ++) {
            double theta = -2.0 * M_PI * static_cast<double>(n) / static_cast<double>(N_THETA_PTS-1);
            double x = pos[i][0] + r_part * cos(theta);
            double y = pos[i][1] + r_part * sin(theta);
            boost::geometry::append(poly.outer(), Point(x, y));
        }

        projections_of_monomers.push(poly);
    }

    std::vector<Polygon> projected_union = {projections_of_monomers.top()};
    projections_of_monomers.pop();

    while (!projections_of_monomers.empty()) {
        boost::geometry::union_(projected_union[0], projections_of_monomers.top(), projected_union);
        projections_of_monomers.pop();
    }

}
