// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <string>
#include <vector>

#include "defs.hpp"
#include "generators.hpp"
#include "style.hpp"
#include "svg.hpp"
#include "utils.hpp"

namespace actsvg {

namespace draw {
/** Draw a polygon object
 *
 * @param polygon_ the polygon points
 * @param id_ is the identification
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 */
static inline svg::object polygon(
    const std::vector<point2> &polygon_, const std::string &id_ = "",
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
    for (auto [x, y] : polygon_) {
        // Add to the real barycenter (without display scaling)
        p._real_barycenter[0] += x;
        p._real_barycenter[1] += y;
        // Add display scaling
        x *= sx;
        y *= sy;
        // Record min/max for the view point
        p._x_range = {std::min(p._x_range[0], x + tx),
                      std::max(p._x_range[1], x + tx)};
        p._y_range = {std::min(p._y_range[0], y + ty),
                      std::max(p._y_range[1], y + ty)};
        // Convert to string attributes
        svg_points_string += std::to_string(x);
        svg_points_string += ",";
        svg_points_string += std::to_string(y);
        svg_points_string += " ";
    }
    // Re-normalize the barycenter
    p._real_barycenter[0] /= polygon_.size();
    p._real_barycenter[1] /= polygon_.size();
    // Fill the points attributes
    p._attribute_map["points"] = svg_points_string;
    // Attach fill, stroke & transform attributes and apply
    p._fill = fill_;
    p._stroke = stroke_;
    p._transform = transform_;
    fill_.attach_attributes(p);
    stroke_.attach_attributes(p);
    transform_.attach_attributes(p);
    // The svg object is now set up
    return p;
}

/** Draw a text object - unconnected
 *
 * @note will perform the y switxh
 *
 * @param p_ is the text position
 * @param tid_ is the text object id
 * @param text_ is the actual text to be drawn
 * @param font_ is the font sytle specification
 * @param transform_ defines the text transform
 **/
static inline svg::object text(
    const point2 &p_, const std::string &tid_,
    const std::vector<std::string> &text_,
    const style::font &font_ = style::font(),
    const style::transform &transform_ = style::transform()) {

    // Create the object, tag it, id it (if given)
    svg::object t;
    t._tag = "text";
    t._id = tid_;
    // Apply the scale
    scalar x = p_[0];
    scalar y = p_[1];
    x *= transform_._scale[0];
    y *= transform_._scale[1];
    // Fill the field
    t._field = text_;
    t._attribute_map["x"] = std::to_string(x);
    t._attribute_map["y"] = std::to_string(-y);
    // barycenter
    t._real_barycenter = {x, -y};

    font_.attach_attributes(t);
    transform_.attach_attributes(t);
    return t;
}

/** Draw a text object - unconnected
 *
 * @param p_ is the text position
 * @param tid_ is the text object id
 * @param text_ is the actual text to be drawn
 * @param font_ is the font sytle specification
 * @param transform_ defines the text transform
 * @param object_ is the connected object
 * @param highlight_ are the hightlighting options
 **/
static inline svg::object connected_text(
    const point2 &p_, const std::string &tid_,
    const std::vector<std::string> &text_, const style::font &font_,
    const style::transform &transform_, const svg::object &object_,
    const std::vector<std::string> &highlight_ = {"mouseover", "mouseout"}) {

    auto t = text(p_, tid_, text_, font_, transform_);

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

/** Draw a connected r-phi-grid in x_y
 *
 * @param r_edges_ are the edges in r
 * @param phi_edges_ are the edges in phi
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 */
static inline std::vector<svg::object> r_phi_grid(
    const std::vector<scalar> &r_edges_, const std::vector<scalar> &phi_edges_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {
    // The list of grid sectors
    std::vector<svg::object> grid_sectors;
    for (size_t ir = 1; ir < r_edges_.size(); ++ir) {
        // Grid svg object
        std::string gs = "g_r";
        gs += std::to_string(ir - 1);
        gs += "_phi";
        for (size_t iphi = 1; iphi < phi_edges_.size(); ++iphi) {
            auto sector_contour = generators::sector_contour(
                r_edges_[ir - 1], r_edges_[ir], phi_edges_[iphi - 1],
                phi_edges_[iphi]);

            auto grid_sector =
                polygon(sector_contour, gs + std::to_string(iphi - 1), fill_,
                        stroke_, transform_);
            scalar r = r_edges_[ir - 1], r_edges_[ir];
            scalar phi = phi_edges_[iphi - 1], phi_edges_[iphi];
            grid_sector._real_barycenter = {r * std::cos(phi),
                                            r * std::sin(phi)};
            grid_sectors.push_back(grid_sector);
        }
    }
    return grid_sectors;
}

/** Draw a connected z-phi-grid in x_y
 *
 * @param z_edges_ are the edges in r
 * @param phi_edges_ are the edges in phi
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 */
static inline std::vector<svg::object> z_phi_grid(
    const std::vector<scalar> &z_edges_, const std::vector<scalar> &phi_edges_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) {

    // The list of grid sectors
    std::vector<svg::object> grid_tiles;
    for (size_t iz = 1; iz < z_edges_.size(); ++iz) {
        // Grid svg object
        std::string gs = "g_z";
        gs += std::to_string(iz - 1);
        gs += "_phi";
        for (size_t iphi = 1; iphi < phi_edges_.size(); ++iphi) {
            std::array<scalar, 2u> llc = {z_edges_[iz - 1],
                                          phi_edges_[iphi - 1]};
            std::array<scalar, 2u> lrc = {z_edges_[iz], phi_edges_[iphi - 1]};
            std::array<scalar, 2u> rrc = {z_edges_[iz], phi_edges_[iphi]};
            std::array<scalar, 2u> rlc = {z_edges_[iz - 1], phi_edges_[iphi]};

            std::vector<std::array<scalar, 2u>> tile = {llc, lrc, rrc, rlc};

            auto grid_tile = polygon(tile, gs + std::to_string(iphi - 1), fill_,
                                     stroke_, transform_);
            grid_tiles.push_back(grid_tile);
        }
    }
    return grid_tiles;
}

/** Method to draw a simple line
 * @note will perform the y switch
 * 
 * @param start_ is the start point of the line
 * @param end_ is the end point of the line
 * @param stroke_ are the stroke parameters
 *
 * @return an svg object for the line
 */
static svg::object line(const point2 &start_, const point2 &end_,
                        const style::stroke &stroke_ = style::stroke()) {
    svg::object l;
    l._tag = "line";
    // Draw the line, remember the sign flip
    l._attribute_map["x1"] = std::to_string(start_[0]);
    l._attribute_map["y1"] = std::to_string(-start_[1]);
    l._attribute_map["x2"] = std::to_string(end_[0]);
    l._attribute_map["y2"] = std::to_string(-end_[1]);
    // We have the range
    l._x_range = {std::min(start_[0], end_[0])};
    l._y_range = {std::min(start_[1], end_[1])};
    // Remember the stroke attributes and add them
    l._stroke = stroke_;
    stroke_.attach_attributes(l);
    return l;
}

/** Marker definition
 *
 *  Arrorws types are: <, <<, <|, |<, |<<, |<|, o, *
 *
 * @param at_ is the poistion of the marker
 * @param marker_ is the marker style
 * @param m_id_ is the marker identification
 * 
 * @return an svg object for the marker
 **/
static inline svg::object marker(const point2 &at_, const style::marker &marker_,
                          const std::string &m_id_ = "x0") {

    svg::object marker_group;
    marker_group._tag = "g";

    std::vector<point2> arrow_head;
    auto size = marker_._size;

    // offset due to measureing
    scalar m_offset = 0.;

    // It's a measure type
    if (marker_._type.substr(0u, 1u) == "|") {
        auto measure_line = line(
            {at_[0], static_cast<scalar>(at_[1] - 2 * size)},
            {at_[0], static_cast<scalar>(at_[1] + 2 * size)}, marker_._stroke);
        marker_._transform.attach_attributes(measure_line);
        marker_group.add_object(measure_line);
        m_offset = -size;
    }
    // Still an arrow to draw
    if (marker_._type.find("<") != std::string::npos) {
        arrow_head = {{at_[0] - size + m_offset, at_[1] - size},
                      {at_[0] + size + m_offset, at_[1]},
                      {at_[0] - size + m_offset, at_[1] + size}};

        // Modify the arrow-type marker
        if (marker_._type == "<<") {
            arrow_head.push_back(
                {static_cast<scalar>(at_[0] - 0.25 * size + m_offset), at_[1]});
        } else if (marker_._type.substr(1u, marker_._type.size()).find("|") ==
                   std::string::npos) {
            arrow_head.push_back(
                {static_cast<scalar>(at_[0] + 1 * size + m_offset), at_[1]});
        }
    }
    // Plot the arrow if not empty
    if (not arrow_head.empty()) {
        auto arrow = polygon(arrow_head, m_id_, marker_._fill, marker_._stroke,
                             marker_._transform);
        marker_group.add_object(arrow);
    }

    return marker_group;
}

/** Draw a measure in z-y
 *
 * @param start_ is the start point of the line
 * @param end_ is the end point of the line
 * @param stroke_ are the stroke parameters
 * @param marker_ are the marker parmeters
 * @param label_ is the label associated
 * @param font_ are the font parameters
 * @param side_x_ is the x offset of the label
 * @param side_y_ is the y offset of the label
 * 
 * @return an svg object for the measurexs
 */
static inline svg::object measure(const point2 &start_, const point2 &end_,
                           const style::stroke &stroke_ = style::stroke(),
                           const style::marker &marker_ = style::marker({"|<"}),
                           const std::string &label_ = "",
                           const style::font &font_ = style::font(),
                           int side_x_ = 1, int side_y_ = 1) {

    // Measure group here we go
    svg::object measure_group;
    measure_group._tag = "g";

    auto mline = line(start_, end_, stroke_);
    measure_group.add_object(mline);

    // Calculate the rotation
    scalar theta = std::atan2(end_[1] - start_[1], -end_[0] + start_[0]);
    scalar theta_deg = theta * 180 / M_PI;

    style::marker lmarker = marker_;
    lmarker._transform = style::transform(
        {end_[0], -end_[1], static_cast<scalar>(theta_deg + 180.)});
    measure_group.add_object(marker({0., 0.}, lmarker));

    style::marker rmarker = marker_;
    rmarker._transform = style::transform({start_[0], -start_[1], theta_deg});
    measure_group.add_object(marker({0., 0.}, rmarker));

    if (not label_.empty()) {
        scalar size = marker_._size;

        scalar x_off = side_x_ * 2 * std::abs(std::sin(theta)) * size;
        scalar y_off = -side_y_ * 2 * std::abs(std::cos(theta)) * size;

        scalar xl = 0.5 * (start_[0] + end_[0]) + x_off;
        scalar yl = 0.5 * (start_[1] + end_[1]) - y_off;
        auto ltext = text({xl, yl}, "t1", {label_}, font_);
        measure_group.add_object(ltext);
    }

    return measure_group;
}

/** Draw an x-y axes system
 * 
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
static inline svg::object x_y_axes(const std::array<scalar, 2> &x_range_,
                            const std::array<scalar, 2> &y_range_,
                            const style::stroke &stroke_ = style::stroke(),
                            const std::string &x_label_ = "",
                            const std::string &y_label_ = "",
                            const style::font &font_ = style::font(),
                            const style::axis_markers<2u> &markers_ = {
                                __standard_axis_markers,
                                __standard_axis_markers}) {

    svg::object axes;
    axes._tag = "g";

    auto x = line({x_range_[0], 0.}, {x_range_[1], 0.}, stroke_);
    auto y = line({0., y_range_[0]}, {0., y_range_[1]}, stroke_);

    axes.add_object(x);
    axes.add_object(y);

    /** Helper method to add marker heads
     *
     * @param p_ is the position of the marker
     * @param b0_ and @param b1_ are the accessors into the marker styles
     * @param rot_ is the rotation parameter
     * @param mid_ is the marker idengification
     *
     * */
    auto add_marker = [&](const point2 &p_, unsigned int b0_, unsigned int b1_,
                          scalar rot_, const std::string &mid_) -> void {
        auto lmarker = markers_[b0_][b1_];
        if (lmarker._type != "none") {
            lmarker._transform = style::transform({p_[0], -p_[1], -rot_});
            axes.add_object(marker({0., 0.}, lmarker, mid_));
        }
    };

    // Add the markers to the arrows
    add_marker({x_range_[0], 0.}, 0, 0, 180., "neg_x_head");
    add_marker({x_range_[1], 0.}, 0, 1, 0., "pos_x_head");
    add_marker({0., y_range_[0]}, 1, 0, -90., "neg_y_head");
    add_marker({0., y_range_[1]}, 1, 1, 90., "pos_y_head");

    // Add the labels: x
    if (not x_label_.empty()) {
        scalar size = markers_[0][1]._size;
        auto xlab =
            text({x_range_[1] + 2 * size, size}, "t1", {x_label_}, font_);
        axes.add_object(xlab);
    }
    // Add the labels: y
    if (not y_label_.empty()) {
        scalar size = markers_[1][1]._size;
        auto ylab =
            text({-size, y_range_[1] + 2 * size}, "t1", {y_label_}, font_);
        axes.add_object(ylab);
    }

    return axes;
}

}  // namespace draw

}  // namespace actsvg
