// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "defs.hpp"
#include "generators.hpp"
#include "style.hpp"
#include "svg.hpp"
#include "utils.hpp"

namespace actsvg {

namespace detail {
/** Helper method to estimate ranges of an object
 *
 * @param _o_ is the svg object in question (to be adapted)
 * @param vertices_ are the input vertices
 **/
static inline void adapt_range(svg::object &_o_,
                               const std::vector<point2> &vertices_) {
    for (const auto &v : vertices_) {
        _o_._x_range = {std::min(_o_._x_range[0], v[0]),
                        std::max(_o_._x_range[1], v[0])};
        _o_._y_range = {std::min(_o_._y_range[0], v[1]),
                        std::max(_o_._y_range[1], v[1])};

        scalar r = std::sqrt(v[0] * v[0] + v[1] * v[1]);
        scalar phi = std::atan2(-v[1], v[0]);
        _o_._r_range = {std::min(_o_._r_range[0], r),
                        std::max(_o_._r_range[1], r)};
        _o_._phi_range = {std::min(_o_._phi_range[0], phi),
                          std::max(_o_._phi_range[1], phi)};
    }
}

}  // namespace detail

/** The draw namespace encapsulates the CORRECT
 * left handed x-y system from the AWKWARD SVG draw system, and applies
 * the transform to it if necessary
 *
 * That is its main purpose (togehter with given allowing to
 * give an identifier an connect objects.
 *
 * */

namespace draw {

/** Method to draw a simple line
 *
 * @note will perform the y switch
 *
 * @param id_ is the identification tag of this line
 * @param start_ is the start point of the line
 * @param end_ is the end point of the line
 * @param stroke_ are the stroke parameters
 *
 * @param transform_ is an optional transform of the object
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the line
 */
static inline svg::object line(
    const std::string &id_, const point2 &start_, const point2 &end_,
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {
    svg::object l;
    l._tag = "line";
    l._id = id_;

    // Apply the transform & scale
    // - the line needs the scale direcly applied
    scalar tx = transform_._tr[0];
    scalar ty = transform_._tr[1];
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];

    scalar st_x = sx * (start_[0] + tx);
    scalar st_y = sx * (-start_[1] - ty);
    scalar en_x = sy * (end_[0] + ty);
    scalar en_y = sy * (-end_[1] - ty);

    l._barycenter =
        utils::barycenter<std::array<scalar, 2>>({{st_x, st_y}, {en_x, en_y}});

    // Draw the line, remember the sign flip
    l._attribute_map["x1"] = utils::to_string(st_x);
    l._attribute_map["y1"] = utils::to_string(st_y);
    l._attribute_map["x2"] = utils::to_string(en_x);
    l._attribute_map["y2"] = utils::to_string(en_y);

    // Adapt the range of this object
    detail::adapt_range(l, {{st_x, st_y}, {en_x, en_y}});

    // Remember the stroke attributes and add them
    l._stroke = stroke_;
    return l;
}

/** Method to draw an arc
 *
 * @param id_ is the identification tag of this line
 * @param r_ the radius
 * @param start_ is the start point of the line
 * @param end_ is the end point of the line
 * @param stroke_ are the stroke parameters
 *
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the line
 */
static inline svg::object arc(
    const std::string &id_, scalar r_, const point2 &start_, const point2 &end_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {
    svg::object a;
    a._tag = "path";
    a._id = id_;

    // Apply the transform & scale
    // - the arc needs the scale directly attached
    scalar tx = transform_._tr[0];
    scalar ty = transform_._tr[1];
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];

    scalar x_min = sx * (start_[0] + tx);
    scalar y_min = -sy * (start_[1] + ty);
    scalar x_max = sx * (end_[0] + tx);
    scalar y_max = -sy * (end_[1] + ty);

    std::string arc_string =
        "M " + utils::to_string(x_min) + " " + utils::to_string(y_min);
    arc_string +=
        " A " + utils::to_string(sx * r_) + " " + utils::to_string(sy * r_);
    arc_string += " 0 0 0 ";
    arc_string += utils::to_string(x_max) + " " + utils::to_string(y_max);
    a._attribute_map["d"] = arc_string;

    // Adapt the range of this object
    detail::adapt_range(a, {{x_min, y_min}, {x_max, y_max}});

    /// @todo add sagitta
    a._barycenter = utils::barycenter<std::array<scalar, 2>>(
        {{x_min, y_min}, {x_max, y_max}});

    // Remember the stroke attributes and add them
    a._stroke = stroke_;
    a._fill = fill_;
    return a;
}

/** Draw a circle object
 *  - will translate into ellipse to allow for a scale
 *
 * @param id_ is the identification
 * @param p_ the position
 * @param r_ the radius
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the circle
 */
static inline svg::object circle(
    const std::string &id_, const point2 &p_, scalar r_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {
    // Create the object, tag it, id it (if given)
    svg::object e;
    e._tag = "ellipse";
    e._id = id_;

    // Apply the transform & scale

    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];
    scalar cx = sx * (p_[0] + transform_._tr[0]);
    scalar cy = sy * (-p_[1] - transform_._tr[1]);

    // Fill the points attributes
    e._attribute_map["cx"] = utils::to_string(cx);
    e._attribute_map["cy"] = utils::to_string(cy);
    e._attribute_map["rx"] = utils::to_string(r_ * sx);
    e._attribute_map["ry"] = utils::to_string(r_ * sy);

    // Adapt the range of this object
    detail::adapt_range(
        e, {{cx - sx * r_, cy - sy * r_}, {cx + sx * r_, cy + sy * r_}});

    e._barycenter = {cx, cy};

    // Attach fill, stroke & transform attributes and apply
    e._fill = fill_;
    e._stroke = stroke_;

    // The svg object is now set up
    return e;
}

/** Draw an ellipse object
 *
 * @param id_ is the identification
 * @param p_ the position
 * @param rs_ the radii
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the ellipse
 */
static inline svg::object ellipse(
    const std::string &id_, const point2 &p_, const std::array<scalar, 2> &rs_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {

    // Create the object, tag it, id it (if given)
    svg::object e;
    e._tag = "ellipse";
    e._id = id_;

    // Apply the transform & scale
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];
    scalar cx = sx * (p_[0] + transform_._tr[0]);
    scalar cy = sy * (-p_[1] - transform_._tr[1]);

    // Fill the points attributes
    e._attribute_map["cx"] = utils::to_string(cx);
    e._attribute_map["cy"] = utils::to_string(cy);
    e._attribute_map["rx"] = utils::to_string(rs_[0] * sx);
    e._attribute_map["ry"] = utils::to_string(rs_[1] * sy);

    // Adapt the range of this object
    detail::adapt_range(e, {{cx - sx * rs_[0], cy - sy * rs_[1]},
                            {cx + sx * rs_[0], cy + sy * rs_[1]}});

    // Barycenter
    e._barycenter = {cx, cy};

    // Attach fill, stroke & transform attributes and apply
    e._fill = fill_;
    e._stroke = stroke_;
    // The svg object is now set up
    return e;
}

/** Draw a polygon object
 *
 * @param id_ is the identification
 * @param polygon_ the polygon points
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the polygon
 */
static inline svg::object polygon(
    const std::string &id_, const std::vector<point2> &polygon_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform())

{
    // Create the object, tag it, id it (if given)
    svg::object p;
    p._tag = "polygon";
    p._id = id_;
    // Apply the scale
    scalar tx = transform_._tr[0];
    scalar ty = transform_._tr[1];
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];
    // Write attributes and measure object size, length
    std::string svg_points_string;
    std::vector<point2> display_vertices;
    display_vertices.reserve(polygon_.size());
    for (auto v : polygon_) {
        scalar alpha = transform_._rot[0];
        if (alpha != 0.) {
            scalar alpha_rad = static_cast<scalar>(alpha / 180. * M_PI);
            v = utils::rotate(v, alpha_rad);
        }

        // Add display scaling
        v[0] *= sx;
        v[1] *= sy;
        // Add scaled * translation
        v[0] += sx * tx;
        v[1] += sy * ty;
        v[1] *= -1;
        // Per vertex range estimation
        detail::adapt_range(p, {v});

        // Convert to string attributes, y-switch
        svg_points_string += utils::to_string(v[0]);
        svg_points_string += ",";
        svg_points_string += utils::to_string(v[1]);
        svg_points_string += " ";
    }
    // Barycenter
    p._barycenter = utils::barycenter(display_vertices);
    // Fill the points attributes
    p._attribute_map["points"] = svg_points_string;
    // Attach fill, stroke & transform attributes and apply
    p._fill = fill_;
    p._stroke = stroke_;
    // The svg object is now set up
    return p;
}

/** Draw a text object - unconnected
 *
 * @param id_ is the text object id
 * @param p_ is the text position
 * @param text_ is the actual text to be drawn
 * @param font_ is the font sytle specification
 * @param transform_ defines the text transform
 *
 * @return an svg object for the text
 *
 **/
static inline svg::object text(
    const std::string &id_, const point2 &p_,
    const std::vector<std::string> &text_,
    const style::font &font_ = style::font(),
    const style::transform &transform_ = style::transform()) {
    // Create the object, tag it, id it (if given)
    svg::object t;
    t._tag = "text";
    t._id = id_;
    // Apply the scale
    scalar x = p_[0];
    scalar y = p_[1];

    x *= transform_._scale[0];
    y *= transform_._scale[1];

    t._fill = font_._fc;

    // Fill the field
    t._field = text_;
    t._attribute_map["x"] = utils::to_string(x);
    t._attribute_map["y"] = utils::to_string(-y);
    t._attribute_map["font-family"] = font_._family;

    t._field_span = font_._size * font_._line_spacing;

    size_t l = 0;
    for (const auto &tl : text_) {
        l = l > tl.size() ? l : tl.size();
    }

    scalar fs = font_._size;

    detail::adapt_range(
        t, {{x, y - l}, {static_cast<scalar>(x + 0.7 * fs * l), y + l}});

    t._barycenter = utils::barycenter<std::array<scalar, 2>>(
        {{x, y - l}, {static_cast<scalar>(x + 0.7 * fs * l), y + l}});

    return t;
}

/** Draw a text object - unconnected
 *
 * @param id_ is the text object id
 * @param p_ is the text position
 * @param text_ is the actual text to be drawn
 * @param font_ is the font sytle specification
 * @param transform_ defines the text transform
 * @param object_ is the connected object
 * @param highlight_ are the hightlighting options
 *
 * @return an svg object with highlight connection
 *
 **/
static inline svg::object connected_text(
    const std::string &id_, const point2 &p_,
    const std::vector<std::string> &text_, const style::font &font_,
    const style::transform &transform_, const svg::object &object_,
    const std::vector<std::string> &highlight_ = {"mouseover", "mouseout"}) {
    auto t = text(id_, p_, text_, font_, transform_);

    t._attribute_map["display"] = "none";

    svg::object on;
    on._tag = "animate";
    on._attribute_map["fill"] = "freeze";
    on._attribute_map["attributeName"] = "display";
    on._attribute_map["from"] = "none";
    on._attribute_map["to"] = "block";
    on._attribute_map["begin"] = object_._id + __d + highlight_[1];

    svg::object off;

    off._tag = "animate";
    off._attribute_map["fill"] = "freeze";
    off._attribute_map["attributeName"] = "display";
    off._attribute_map["to"] = "none";
    off._attribute_map["from"] = "block";
    off._attribute_map["begin"] = object_._id + __d + highlight_[0];

    // Store the animation
    t._sub_objects.push_back(on);
    t._sub_objects.push_back(off);
    return t;
}

/** Draw a tiled cartesian grid - ready for connecting
 *
 * @param id_ the grid identification
 * @param l0_edges_ are the edges in l0
 * @param l1_edges_ are the edges in l1
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @return a simple cartesian grid
 */
static inline svg::object cartesian_grid(
    const std::string &id_, const std::vector<scalar> &l0_edges_,
    const std::vector<scalar> &l1_edges_,
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {
    // The grid group object
    svg::object grid;
    grid._tag = "g";
    grid._id = id_;

    scalar l0_min = l0_edges_[0];
    scalar l0_max = l0_edges_[l0_edges_.size() - 1];

    scalar l1_min = l1_edges_[0];
    scalar l1_max = l1_edges_[l1_edges_.size() - 1];

    for (auto [i0, l0] : utils::enumerate(l0_edges_)) {
        grid.add_object(line(id_ + "_l0_" + std::to_string(i0), {l0, l1_min},
                             {l0, l1_max}, stroke_, transform_));
    }

    for (auto [i1, l1] : utils::enumerate(l1_edges_)) {
        grid.add_object(line(id_ + "_l1_" + std::to_string(i1), {l0_min, l1},
                             {l0_max, l1}, stroke_, transform_));
    }

    return grid;
}

/** Draw a tiled cartesian grid - ready for connecting
 *
 * @param id_ the grid identification
 * @param l0_edges_ are the edges in l0
 * @param l1_edges_ are the edges in l1
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note - the single objects of the grid can be found by
 * their unique name "id_+_X_Y" when X is the bin in the
 * first local and j the bin in the second local coordinate
 *
 * @return a tiled grid with internal objects
 */
static inline svg::object tiled_cartesian_grid(
    const std::string &id_, const std::vector<scalar> &l0_edges_,
    const std::vector<scalar> &l1_edges_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {
    // The grid group object
    svg::object grid;
    grid._tag = "g";
    grid._id = id_;

    // The list of grid sectors
    for (size_t il0 = 1; il0 < l0_edges_.size(); ++il0) {
        // Grid svg object
        std::string gs = id_ + "_";
        gs += std::to_string(il0 - 1);
        gs += "_";
        for (size_t il1 = 1; il1 < l1_edges_.size(); ++il1) {
            std::array<scalar, 2u> llc = {l0_edges_[il0 - 1],
                                          l1_edges_[il1 - 1]};
            std::array<scalar, 2u> lrc = {l0_edges_[il0], l1_edges_[il1 - 1]};
            std::array<scalar, 2u> rrc = {l0_edges_[il0], l1_edges_[il1]};
            std::array<scalar, 2u> rlc = {l0_edges_[il0 - 1], l1_edges_[il1]};

            std::vector<std::array<scalar, 2u>> tile = {llc, lrc, rrc, rlc};

            auto grid_tile = polygon(gs + std::to_string(il1 - 1), tile, fill_,
                                     stroke_, transform_);
            grid.add_object(grid_tile);
        }
    }
    return grid;
}

/** Draw a simple fan grid
 *
 * @param id_ the grid identification
 * @param x_low_edges_ are the edges in x at low y
 * @param x_high_edges_ are the edges in x at high y
 * @param y_edges_ are the edges in phi
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @return a simple faned grid structure
 */
static inline svg::object fan_grid(
    const std::string &id_, const std::vector<scalar> &x_low_edges_,
    const std::vector<scalar> &x_high_edges_,
    const std::vector<scalar> &y_edges_,
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) noexcept(false) {
    // The list of grid sectors
    svg::object grid;
    grid._tag = "g";
    grid._id = id_;

    if (x_low_edges_.size() != x_high_edges_.size()) {
        throw std::invalid_argument(
            "fan_grid: mismatch in low/high edge numbers");
    }

    scalar x_low_min = x_low_edges_[0];
    scalar x_low_max = x_low_edges_[x_low_edges_.size() - 1];

    scalar x_high_min = x_high_edges_[0];
    scalar x_high_max = x_high_edges_[x_high_edges_.size() - 1];

    scalar y_min = y_edges_[0];
    scalar y_max = y_edges_[y_edges_.size() - 1];

    // Calculate slopes
    scalar x_low_slope = (x_high_min - x_low_min) / (y_max - y_min);
    scalar x_high_slope = (x_high_max - x_low_max) / (y_max - y_min);

    for (const auto [ix, x] : utils::enumerate(x_low_edges_)) {
        scalar x_min = x;
        scalar x_max = x_high_edges_[ix];
        grid.add_object(line(id_ + "_l0_" + std::to_string(ix), {x_min, y_min},
                             {x_max, y_max}, stroke_, transform_));
    }

    for (const auto &[iy, y] : utils::enumerate(y_edges_)) {
        scalar x_min = x_low_min + (y - y_min) * x_low_slope;
        scalar x_max = x_low_max + (y - y_min) * x_high_slope;
        grid.add_object(line(id_ + "_l1_" + std::to_string(iy), {x_min, y},
                             {x_max, y}, stroke_, transform_));
    }

    return grid;
}

/** Draw a simple fan grid
 *
 * @param id_ the grid identification
 * @param x_low_edges_ are the edges in x at low y
 * @param x_high_edges_ are the edges in x at high y
 * @param y_edges_ are the edges in phi
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note - the single objects of the grid can be found by
 * their unique name "id_+_X_Y" when X is the bin in the
 * first local and j the bin in the second local coordinate
 *
 * @return a tiled grid with contained individual cells
 */
static inline svg::object tiled_fan_grid(
    const std::string &id_, const std::vector<scalar> &x_low_edges_,
    const std::vector<scalar> &x_high_edges_,
    const std::vector<scalar> &y_edges_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) noexcept(false) {
    svg::object grid;
    grid._tag = "g";
    grid._id = id_;

    // The list of grid sectors
    std::vector<svg::object> grid_tiles;

    if (x_low_edges_.size() != x_high_edges_.size()) {
        throw std::invalid_argument(
            "fan_grid: mismatch in low/high edge numbers");
    }

    scalar y_min = y_edges_[0];
    scalar y_max = y_edges_[y_edges_.size() - 1];

    std::vector<scalar> slopes;
    for (auto [ix, xl] : utils::enumerate(x_low_edges_)) {
        scalar xh = x_high_edges_[ix];
        slopes.push_back((xh - xl) / (y_max - y_min));
    }

    for (auto [iy, yl] : utils::enumerate(y_edges_)) {
        if (iy + 1 == y_edges_.size()) {
            continue;
        }
        scalar yh = y_edges_[iy + 1];
        for (auto [ix, xll] : utils::enumerate(x_low_edges_)) {
            if (ix + 1 == x_low_edges_.size()) {
                continue;
            }
            scalar xlr = x_low_edges_[ix + 1];

            scalar x_l_slope = slopes[ix];
            scalar x_r_slope = slopes[ix + 1];

            // the corrected position given the y values
            scalar xll_c = xll + x_l_slope * (yl - y_min);
            scalar xlr_c = xlr + x_r_slope * (yl - y_min);

            scalar xhl_c = xll + x_l_slope * (yh - y_min);
            scalar xhr_c = xlr + x_r_slope * (yh - y_min);

            std::string g_n =
                id_ + "_" + std::to_string(ix) + "_" + std::to_string(iy);
            grid.add_object(draw::polygon(
                g_n, {{xll_c, yl}, {xlr_c, yl}, {xhr_c, yh}, {xhl_c, yh}},
                fill_, stroke_, transform_));
        }
    }

    return grid;
}

/** Draw a simple polar grid
 *
 * @param id_ the grid identification
 * @param r_edges_ are the edges in r
 * @param phi_edges_ are the edges in phi
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @return a simple polar grid object
 */
static inline svg::object polar_grid(
    const std::string &id_, const std::vector<scalar> &r_edges_,
    const std::vector<scalar> &phi_edges_,
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {
    // The list of grid sectors
    svg::object grid;
    grid._tag = "g";
    grid._id = id_;

    scalar r_min = r_edges_[0];
    scalar r_max = r_edges_[r_edges_.size() - 1];

    scalar phi_min = phi_edges_[0];
    scalar phi_max = phi_edges_[phi_edges_.size() - 1];

    scalar cos_phi_min = std::cos(phi_min);
    scalar sin_phi_min = std::sin(phi_min);
    scalar cos_phi_max = std::cos(phi_max);
    scalar sin_phi_max = std::sin(phi_max);

    style::fill fill;

    for (const auto [ir, r] : utils::enumerate(r_edges_)) {
        grid.add_object(draw::arc(id_ + "_r_" + std::to_string(ir), r,
                                  {r * cos_phi_min, r * sin_phi_min},
                                  {r * cos_phi_max, r * sin_phi_max}, fill,
                                  stroke_, transform_));
    }

    for (const auto &[iphi, phi] : utils::enumerate(phi_edges_)) {
        scalar cos_phi = std::cos(phi);
        scalar sin_phi = std::sin(phi);
        grid.add_object(draw::line(id_ + "_l_" + std::to_string(iphi),
                                   {r_min * cos_phi, r_min * sin_phi},
                                   {r_max * cos_phi, r_max * sin_phi}, stroke_,
                                   transform_));
    }
    return grid;
}

/** Draw a connected polar grid
 *
 * @param id_ the grid identification
 * @param r_edges_ are the edges in r
 * @param phi_edges_ are the edges in phi
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note - the single objects of the grid can be found by
 * their unique name "id_+_X_Y" when X is the bin in the
 * first local and j the bin in the second local coordinate
 *
 * @return a tiled polar grid in individual objects
 */
static inline svg::object tiled_polar_grid(
    const std::string &id_, const std::vector<scalar> &r_edges_,
    const std::vector<scalar> &phi_edges_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {
    svg::object grid;
    grid._tag = "g";
    grid._id = id_;

    for (size_t ir = 1; ir < r_edges_.size(); ++ir) {
        // Grid svg object
        std::string gs = id_ + "_";
        gs += std::to_string(ir - 1);
        gs += "_";
        for (size_t iphi = 1; iphi < phi_edges_.size(); ++iphi) {
            auto sector_contour = generators::sector_contour(
                r_edges_[ir - 1], r_edges_[ir], phi_edges_[iphi - 1],
                phi_edges_[iphi]);

            auto grid_sector =
                polygon(gs + std::to_string(iphi - 1), sector_contour, fill_,
                        stroke_, transform_);
            grid.add_object(grid_sector);
        }
    }
    return grid;
}

/** Marker definition
 *
 *  Arrorws types are: <, <<, <|, |<, |<<, |<|, o, x, *
 * @param id_ is the marker identification
 * @param at_ is the poistion of the marker
 * @param marker_ is the marker style
 * @param rot_ is the rotation in [pi,phi)]
 *
 * @return an svg object for the marker group
 **/
static inline svg::object marker(const std::string &id_, const point2 &at_,
                                 const style::marker &marker_,
                                 scalar rot_ = 0.) {
    svg::object marker_group;
    marker_group._tag = "g";

    std::vector<point2> arrow_head;
    auto size = marker_._size;

    // Offset due to measureing
    scalar m_offset = 0.;

    // It's a measure type
    if (marker_._type.substr(0u, 1u) == "|") {
        auto measure_line = line(
            id_ + "_line", {at_[0], static_cast<scalar>(at_[1] - 2 * size)},
            {at_[0], static_cast<scalar>(at_[1] + 2 * size)}, marker_._stroke);
        marker_group.add_object(measure_line);
        m_offset = -size;
    }
    // Still an arrow to draw
    if (marker_._type.find("<") != std::string::npos) {
        arrow_head = {{at_[0] - size + m_offset, at_[1] - size},
                      {at_[0] + size + m_offset, at_[1]},
                      {at_[0] - size + m_offset, at_[1] + size}};

        // Modify the arrow-type marker
        if (marker_._type.find("<<") != std::string::npos) {
            // Filled arrow-head
            arrow_head.push_back(
                {static_cast<scalar>(at_[0] - 0.25 * size + m_offset), at_[1]});
        } else if (marker_._type.substr(1u, marker_._type.size()).find("|") ==
                   std::string::npos) {
            arrow_head.push_back(
                {static_cast<scalar>(at_[0] + 1 * size + m_offset), at_[1]});
        }
    } else if (marker_._type.find("o") != std::string::npos) {
        // A dot marker
        svg::object dot =
            circle(id_, at_, 0.5 * size, marker_._fill, marker_._stroke);
        marker_group.add_object(dot);
    } else if (marker_._type.find("x") != std::string::npos) {
        scalar a_x = at_[0];
        scalar a_y = at_[1];
        scalar h_s = 0.5 * size;
        marker_group.add_object(
            line(id_ + "_ml0", {a_x - h_s, a_y - h_s}, {a_x + h_s, a_y + h_s}));
        marker_group.add_object(
            line(id_ + "_ml1", {a_x - h_s, a_y + h_s}, {a_x + h_s, a_y - h_s}));
    }

    // Plot the arrow if not empty
    if (not arrow_head.empty()) {
        auto arrow = polygon(id_, arrow_head, marker_._fill, marker_._stroke);
        marker_group.add_object(arrow);
    }

    // We need to rotate the marker group
    if (rot_ != 0.) {
        style::transform group_transform = style::transform{};
        group_transform._rot = {static_cast<scalar>(rot_ * 180. / M_PI), at_[0],
                                at_[1]};
        marker_group._transform = group_transform;
    }

    return marker_group;
}

/** Draw a measure in z-y
 *
 * @param id_ is the identification tag of this object
 * @param start_ is the start point of the line
 * @param end_ is the end point of the line
 * @param stroke_ are the stroke parameters
 * @param start_marker_ are the marker parmeters at start
 * @param end_marker_ are the marker parmeters at start
 * @param font_ are the font parameters
 * @param label_ is the label associated
 * @param label_pos_ is the label position
 *
 * @return an svg object for the measurexs
 */
static inline svg::object measure(
    const std::string &id_, const point2 &start_, const point2 &end_,
    const style::stroke &stroke_ = style::stroke(),
    const style::marker &start_marker_ = style::marker({"|<"}),
    const style::marker &end_marker_ = style::marker({"|<"}),
    const style::font &font_ = style::font(), const std::string &label_ = "",
    const point2 &label_pos_ = {0., 0.}) {
    // Measure group here we go
    svg::object measure_group;
    measure_group._tag = "g";
    measure_group._id = id_;

    auto mline = line(id_ + "_line", start_, end_, stroke_);
    measure_group.add_object(mline);

    // Calculate the rotation
    scalar theta = std::atan2(end_[1] - start_[1], end_[0] - start_[0]);
    if (std::abs(end_[1] - start_[1]) <
        std::numeric_limits<scalar>::epsilon()) {
        theta = 0.;
    }

    measure_group.add_object(marker(id_ + "_start_tag", {start_[0], start_[1]},
                                    start_marker_,
                                    M_PI + static_cast<scalar>(theta)));
    measure_group.add_object(marker(id_ + "_end_tag", {end_[0], end_[1]},
                                    end_marker_, static_cast<scalar>(theta)));

    if (not label_.empty()) {
        auto ltext = text(id_ + "_label", label_pos_, {label_}, font_);
        measure_group.add_object(ltext);
    }
    return measure_group;
}

/** Draw an arc measure
 *
 * @param id_ is the identification tag of this object
 * @param r_ the radius
 * @param start_ is the start point of the line, per definition with smaller phi
 * @param end_ is the end point of the line, defines the marker
 * @param stroke_ are the stroke parameters
 * @param start_marker_ are the marker parmeters
 * @param end_marker_ are the marker parmeters
 * @param font_ are the font parameters
 * @param label_ is the label associated
 * @param label_pos_ is the label position
 *
 * @return an svg object for the measures
 */
static inline svg::object arc_measure(
    const std::string &id_, scalar r_, const point2 &start_, const point2 &end_,
    const style::stroke &stroke_ = style::stroke(),
    const style::marker &start_marker_ = style::marker(),
    const style::marker &end_marker_ = style::marker({"|<"}),
    const style::font &font_ = style::font(), const std::string &label_ = "",
    const point2 &label_pos_ = {0., 0.}) {
    // Measure group here we go
    svg::object measure_group;
    measure_group._tag = "g";
    measure_group._id = id_;

    measure_group.add_object(
        arc((id_ + "_arc"), r_, start_, end_, style::fill(), stroke_));

    // Arrow is at end point
    if (not start_marker_._type.empty() and
        start_marker_._type != std::string("none")) {
        scalar theta_start = atan2(start_[1], start_[0]);
        measure_group.add_object(
            marker(id_ + "_start_tag", {start_[0], start_[1]}, start_marker_,
                   static_cast<scalar>(theta_start - 0.5 * M_PI)));
    }

    scalar theta_end = atan2(end_[1], end_[0]);
    measure_group.add_object(
        marker(id_ + "_end_tag", {end_[0], end_[1]}, end_marker_,
               static_cast<scalar>(theta_end + 0.5 * M_PI)));

    if (not label_.empty()) {
        auto ltext = text(id_ + "_label", label_pos_, {label_}, font_);
        measure_group.add_object(ltext);
    }

    return measure_group;
}

/** Draw an x-y axes system
 *
 * @param id_ is the id tag of this object
 * @param x_range_ is the x range of the axes to be drawn
 * @param y_range_ is the y range of the axes to be drawn
 * @param stroke_ are the stroke parameters
 * @param x_label_ is the x label of the axis system
 * @param y_label_ is the y label of the axis system
 * @param font_ are the font parameters
 * @param markers_ are the 4 markers on each axis end
 *
 * @return an svg object representing the axes
 */
static inline svg::object x_y_axes(
    const std::string &id_, const std::array<scalar, 2> &x_range_,
    const std::array<scalar, 2> &y_range_,
    const style::stroke &stroke_ = style::stroke(),
    const std::string &x_label_ = "", const std::string &y_label_ = "",
    const style::font &font_ = style::font(),
    const style::axis_markers<2u> &markers_ = {__standard_axis_markers,
                                               __standard_axis_markers}) {
    svg::object axes;
    axes._tag = "g";
    axes._id = id_;

    auto x =
        line(id_ + "_x_axis", {x_range_[0], 0.}, {x_range_[1], 0.}, stroke_);
    auto y =
        line(id_ + "_y_axis", {0., y_range_[0]}, {0., y_range_[1]}, stroke_);

    axes.add_object(x);
    axes.add_object(y);

    /** Helper method to add marker heads
     *
     * @param p_ is the position of the marker
     * @param b0_ and @param b1_ are the accessors into the marker styles
     * @param rot_ is the rotation parameter
     * @param mid_ is the marker identification
     *
     * */
    auto add_marker = [&](const point2 &p_, unsigned int b0_, unsigned int b1_,
                          scalar rot_, const std::string &mid_) -> void {
        auto lmarker = markers_[b0_][b1_];
        if (not lmarker._type.empty()) {
            axes.add_object(marker(mid_, {p_[0], p_[1]}, lmarker, rot_));
        }
    };

    // Add the markers to the arrows
    add_marker({x_range_[0], 0.}, 0, 0, M_PI, id_ + "_neg_x_head");
    add_marker({x_range_[1], 0.}, 0, 1, 0., id_ + "pos_x_head");
    add_marker({0., y_range_[0]}, 1, 0, -0.5 * M_PI, id_ + "neg_y_head");
    add_marker({0., y_range_[1]}, 1, 1, 0.5 * M_PI, id_ + "pos_y_head");

    // Add the labels: x
    if (not x_label_.empty()) {
        scalar size = markers_[0][1]._size;
        auto xlab = text(id_ + "_x_label", {x_range_[1] + 2 * size, size},
                         {x_label_}, font_);
        axes.add_object(xlab);
    }
    // Add the labels: y
    if (not y_label_.empty()) {
        scalar size = markers_[1][1]._size;
        auto ylab = text(id_ + "_y_label", {-size, y_range_[1] + 2 * size},
                         {y_label_}, font_);
        axes.add_object(ylab);
    }

    return axes;
}

/** Place a copy with different attributes via xling
 *
 * @param id_ the identification of this surface
 * @param ro_ the reference object
 * @param f_ the fill of the new object
 * @param s_ the stroke of the new object
 * @param t_ the transform of the new object
 *
 * @return the create object from reference
 */
static inline svg::object from_template(
    const std::string &id_, const svg::object &ro_,
    const style::fill &f_ = style::fill(),
    const style::stroke &s_ = style::stroke(),
    const style::transform t_ = style::transform()) {
    // Create new svg object
    svg::object nsvg;
    nsvg._sterile = true;
    nsvg._id = id_;
    nsvg._tag = "g";

    // Refer to as the linker object
    svg::object use_obj;
    use_obj._tag = "use";
    use_obj._id = id_ + "_use";
    use_obj._attribute_map["xlink:href"] = "#" + ro_._id;
    use_obj._definitions.push_back(ro_);

    // Barycenter is shifted
    use_obj._barycenter = {(ro_._barycenter[0] + t_._tr[0]),
                           (ro_._barycenter[1] + t_._tr[1])};

    // Set the fill attributes
    use_obj._fill = f_;
    use_obj._stroke = s_;
    use_obj._transform = t_;

    // Add it to the main object
    nsvg.add_object(use_obj);

    return nsvg;
}

}  // namespace draw

}  // namespace actsvg
