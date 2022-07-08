// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "defs.hpp"

namespace actsvg {

/** The view convert 3D point collections as given from
 * the geometry output, and convert them into 2D contours that can be displayed.
 *
 * @note Given the mathemtatical negative orientation of the x/y system in SVG,
 * it will flip the y coordinate accordingly.
 *
 */
namespace views {

/// A contour is a 2-dimensional object
using contour = std::vector<point2>;

/// Single view per module/surface
struct x_y {
    /// Make it screen obvious
    std::array<std::string, 2> _axis_names = {"x", "y"};

    /** A single planar view operator, asuming a contour in
     *  x/y plane
     *
     * @tparam point3_container is a 3D-point container, where elements
     * of a single p3 object can be accessed via [] operators
     *
     * @param vertices_ the vertices of the surface, they are
     * assumed to be in global/display coorindates
     *
     * @return a 2D contour array
     */
    template <typename point3_container>
    contour operator()(const point3_container &vertices_) const {
        contour c;
        c.reserve(vertices_.size());
        for (const auto &v : vertices_) {
            // flip y coordinate */
            c.push_back({static_cast<scalar>(v[0]), static_cast<scalar>(v[1])});
        }
        return c;
    }
};

/// z_r projection view
struct z_r {

    /// Make it screen obvious
    std::array<std::string, 2> _axis_names = {"z", "r"};

    /** A r-z view operator
     *
     * @tparam point3_container is a 3D-point container, where elements
     * of a single p3 object can be accessed via [] operators
     *
     * @param vertices_ the vertices of the surface
     *
     * @return a 2D contour array
     */
    template <typename point3_container>
    contour operator()(const point3_container &vertices_) const {
        // Fill the contour accordingly
        contour c;
        c.reserve(vertices_.size());
        for (const auto &v : vertices_) {
            scalar r = std::sqrt(v[0] * v[0] + v[1] * v[1]);
            c.push_back({static_cast<scalar>(v[2]), r});
        }
        return c;
    }
};

// z-phi projection view
struct z_phi {
    /// Switch the phi protection on (wrapping in phi detection)
    bool _protect_phi = true;

    /// Make it screen obvious
    std::array<std::string, 2> _axis_names = {"z", "phi"};

    /** A z-phi view operator
     *
     * @tparam point3_container is a 3D-point container, where elements
     * of a single p3 object can be accessed via [] operators
     *
     * @param vertices_ the vertices of the surface
     *
     *
     * @return a 2D contour array
     */
    template <typename point3_container>
    contour operator()(const point3_container &vertices_) const {

        // Fill the contour accordingly
        contour c;
        c.reserve(vertices_.size());
        std::vector<scalar> phi_values;
        phi_values.reserve(vertices_.size());

        for (const auto &v : vertices_) {
            scalar phi = std::atan2(v[1], v[0]);
            phi_values.push_back(phi);
            c.push_back({static_cast<scalar>(v[2]), phi});
        }
        // Run the phi detection and protection
        if (_protect_phi) {
            auto min_phi =
                std::min_element(phi_values.begin(), phi_values.end());
            auto max_phi =
                std::max_element(phi_values.begin(), phi_values.end());

            if ((*min_phi) < 0. and (*max_phi) > 0. and
                ((*max_phi) - (*min_phi)) > M_PI) {
                for (auto &cv : c) {
                    if (cv[1] < 0.) {
                        cv[1] += 2 * M_PI;
                    }
                }
            }
        }
        return c;
    }
};

// z_rphi projection view
struct z_rphi {
    scalar fixed_r = std::numeric_limits<scalar>::quiet_NaN();

    /// Make it screen obvious
    std::array<std::string, 2> _axis_names = {"z", "r · phi"};

    /** A z-rphi view operator
     *
     * @tparam point3_container is a 3D-point container, where elements
     * of a single p3 object can be accessed via [] operators
     *
     * @param vertices_ the vertices of the surface
     *
     * @return a 2D contour array
     */
    template <typename point3_container>
    contour operator()(const point3_container &vertices_) const {
        // Fill the contour accordingly
        contour c;
        c.reserve(vertices_.size());
        for (const auto &v : vertices_) {
            scalar r = fixed_r;
            if (std::isnan(r)) {
                r = std::sqrt(v[0] * v[0] + v[1] * v[1]);
            }
            scalar phi = std::atan2(v[1], v[0]);
            c.push_back({static_cast<scalar>(v[2]), r * phi});
        }
        return c;
    }
};

}  // namespace views

}  // namespace actsvg
