// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
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

/// A scene description
struct scene {
    /// A camera position
    point2 _camera = {0., 0.};
    /// A view angle
    point2 _view = {0., -1.};
    /// A strict view, i.e. only see along the view direction
    bool _strict = false;
    /// What's actually shown
    std::array<std::array<scalar, 2>, 2> _range = {
        std::array<scalar, 2>{std::numeric_limits<scalar>::lowest(),
                              std::numeric_limits<scalar>::max()},
        std::array<scalar, 2>{std::numeric_limits<scalar>::lowest(),
                              std::numeric_limits<scalar>::max()}};
};

/// Single view per module/surface
struct x_y {

    /// Scene setting
    scene _scene = scene();

    /// Make it screen obvious
    std::array<std::string, 2> _axis_names = {"x", "y"};

    /** A planar view operator, asuming a contour in x/y plane
     *
     * @tparam point3 is a 3D-point type assuming [] operators exit
     *
     * @param point_ the point global coordinates
     *
     * @return a 2D point
     */
    template <typename point3_type>
    point2 point(const point3_type &point_) const {
        return point2{static_cast<scalar>(point_[0]),
                      static_cast<scalar>(point_[1])};
    }

    /** A planar view operator, asuming a contour in x/y plane
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
    contour path(const point3_container &vertices_) const {
        contour c;
        c.reserve(vertices_.size());
        for (const auto &v : vertices_) {
            // flip y coordinate */
            c.push_back(point(v));
        }
        return c;
    }

    /** A planar view operator, asuming a contour in x/y plane
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
        return path(vertices_);
    }
};

/// z_r projection view
struct z_r {

    /// Scene setting
    scene _scene = scene();

    /// Make it screen obvious
    std::array<std::string, 2> _axis_names = {"z", "r"};

    /** A r-z view operator for single objects
     *
     * @tparam point3 is a 3D-point type assuming [] operators exit
     *
     * @param point_ the point global coordinates
     *
     * @return a 2D point
     */
    template <typename point3_type>
    point2 point(const point3_type &point_) const {
        scalar r = std::sqrt(point_[0] * point_[0] + point_[1] * point_[1]);
        return point2{static_cast<scalar>(point_[2]), r};
    }

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
    contour path(const point3_container &vertices_) const {
        // Fill the contour accordingly
        contour c;
        c.reserve(vertices_.size());
        for (const auto &v : vertices_) {
            c.push_back(point(v));
        }
        return c;
    }

    /** A planar view operator
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
        return path(vertices_);
    }
};

// z-phi projection view
struct z_phi {

    /// Scene setting
    scene _scene = scene();

    /// Switch the phi protection on (wrapping in phi detection)
    bool _protect_phi = true;

    /// Make it screen obvious
    std::array<std::string, 2> _axis_names = {"z", "phi"};

    /** A z-phi view operator
     *
     * @tparam point3 is a 3D-point type assuming [] operators exit
     *
     * @param point_ the point global coordinates
     *
     * @return a 2D point
     */
    template <typename point3_type>
    point2 point(const point3_type &point_) const {
        scalar phi = std::atan2(point_[1], point_[0]);
        return point2{static_cast<scalar>(point_[2]), phi};
    }

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
    contour path(const point3_container &vertices_) const {

        // Fill the contour accordingly
        contour c;
        c.reserve(vertices_.size());
        std::vector<scalar> phi_values;
        phi_values.reserve(vertices_.size());

        for (const auto &v : vertices_) {
            auto p2 = point(v);
            phi_values.push_back(v[1u]);
            c.push_back(p2);
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
                        cv[1] += 2_scalar * pi;
                    }
                }
            }
        }
        return c;
    }

    /** A planar view operator
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
        return path(vertices_);
    }
};

// z_rphi projection view
struct z_rphi {

    /// Scene setting
    scene _scene = scene();

    /// A fixed radius is needed for this view
    scalar _fixed_r = std::numeric_limits<scalar>::quiet_NaN();

    /// Make it screen obvious
    std::array<std::string, 2> _axis_names = {"z", "r · phi"};

    /** A z-rphi view operator
     *
     * @tparam point3 is a 3D-point type assuming [] operators exit
     *
     * @param point_ the point global coordinates
     *
     * @return a 2D point
     */
    template <typename point3_type>
    point2 point(const point3_type &point_) const {
        scalar r = _fixed_r;
        if (std::isnan(r)) {
            r = std::sqrt(point_[0] * point_[0] + point_[1] * point_[1]);
        }
        scalar phi = std::atan2(point_[1], point_[0]);
        return point2{static_cast<float>(point_[2]), r * phi};
    }

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
    contour path(const point3_container &vertices_) const {
        // Fill the contour accordingly
        contour c;
        c.reserve(vertices_.size());
        for (const auto &v : vertices_) {
            c.push_back(point(v));
        }
        return c;
    }

    /** A planar view operator
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
        return path(vertices_);
    }
};

}  // namespace views

}  // namespace actsvg
