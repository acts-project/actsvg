// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"
#include "actsvg/styles/defaults.hpp"
#include "actsvg/display/geometry.hpp"

namespace actsvg {

using namespace defaults;

namespace display {

/** Helper method to calculate the center
 *
 * @param vs are the vertices that build up this module
 *
 * @return the string
 **/
template <typename point3_container>
std::string center_string(const point3_container& vs) {
    std::string c_str = "center = (";
    scalar c_x = 0;
    scalar c_y = 0;
    scalar c_z = 0;
    for (auto& v : vs) {
        c_x += v[0];
        c_y += v[1];
        c_z += v[2];
    }
    c_x /= vs.size();
    c_y /= vs.size();
    c_z /= vs.size();
    return c_str + std::to_string(c_x) + "," + std::to_string(c_y) + ", " +
           std::to_string(c_z) + ")";
}

/** Helper method to prepare axis for a view point
 *
 * @param first_ is the first axis
 * @param second_ is the second axis
 * @param sx_ is the first axis scale
 * @param sy_ is the second axis scale
 * @param ax_ is the first axis addon
 * @param ay_ is the second axis addon
 *
 * @return the marker size as 1 percent of the range
 **/
static inline void prepare_axes(std::array<scalar, 2>& first_,
                                std::array<scalar, 2>& second_, scalar sx_,
                                scalar sy_, scalar ax_ = 0., scalar ay_ = 0.) {
    // Add some extra space for the axis
    first_[0] *= sx_;
    first_[1] *= sx_;
    first_[0] -= ax_;
    first_[1] += ax_;

    second_[0] *= sy_;
    second_[1] *= sy_;
    second_[0] -= ay_;
    second_[1] += ay_;
}

/** Helper method to scale the axis accordingly
 *
 * @param cc_ is an iterable container of contours
 *
 * @return a view range
 **/
template <typename contour_container>
static std::array<std::array<scalar, 2>, 2> view_range(
    const contour_container& cc_) {

    std::array<scalar, 2> x_range = __e_object._x_range;
    std::array<scalar, 2> y_range = __e_object._y_range;

    for (auto& c : cc_) {
        for (auto& v : c) {
            x_range[0] = std::min(v[0], x_range[0]);
            x_range[1] = std::max(v[0], x_range[1]);
            y_range[0] = std::min(v[1], y_range[0]);
            y_range[1] = std::max(v[1], y_range[1]);
        }
    }
    return {x_range, y_range};
}

/** Helper method to process the modules, estimate the scale and the axes
 *
 * @param v_ volume of the detector
 * @param view_ the view used for this
 * @param sh_ the sheet size
 *
 * @returns the modules, a scale transform & the axes
 **/
template <typename volume_type, typename view_type>
std::tuple<std::vector<svg::object>, style::transform,
           std::array<std::array<scalar, 2>, 2> >
process_modules(const volume_type& v_, const view_type& view_,
                const std::array<scalar, 2>& sh_ = {600., 600.}) {

    using surface_type = typename volume_type::surface_type;

    // Axis range & pre-loop
    std::vector<views::contour> contours;
    contours.reserve(v_._surfaces.size());
    for (auto [is, s] : utils::enumerate(v_._surfaces)) {
        auto surface_contour = view_(s._vertices);
        contours.push_back(surface_contour);
    }

    // Get the scaling right
    auto axes = display::view_range(contours);
    scalar s_x = sh_[0] / (axes[0][1] - axes[0][0]);
    scalar s_y = sh_[1] / (axes[1][1] - axes[1][0]);

    // Create the scale transform
    style::transform scale_transform;
    scale_transform._scale = {s_x, s_y};

    // Draw the modules and estimate axis ranges
    std::vector<svg::object> modules;
    modules.reserve(contours.size());
    for (auto [ic, c] : utils::enumerate(contours)) {
        surface_type draw_surface = v_._surfaces[ic];
        draw_surface._transform._scale = {s_x, s_y};
        auto surface_module = display::surface(draw_surface._name, draw_surface, view_);
        modules.push_back(surface_module);
    }

    // Prepare the axis for the view range
    prepare_axes(axes[0], axes[1], s_x, s_y, 30., 30.);

    return {modules, scale_transform, axes};
}

/** Helper method to connect the surface sheets to the
 * surfaces of the layer_sheets
 *
 * @tparam volume_type the type of volume (templated on point3_container)
 *
 * @param v_ the input volume
 * @param templates_ the given module templtes
 * @param o_ the object to which they are attached
 * @param yt_ is the position of the title text in y
 *
 **/
template <typename volume_type>
void connect_surface_sheets(const volume_type& v_,
                            std::vector<svg::object>& templates_,
                            svg::object& o_, scalar yt_ = 200.) {
    // Now create an item per surface
    for (auto [is, s] : utils::enumerate(v_._surfaces)) {
        std::string sid = s._name;

        // Template copy
        size_t it = v_._templates[is];
        if (it >= v_._templates.size()) {
            continue;
        }
        svg::object s_sheet_s = templates_[v_._templates[is]];

        s_sheet_s._attribute_map["display"] = "none";

        auto surface_info = draw::text("info_" + s._name, {0, yt_}, s._info);
        s_sheet_s.add_object(surface_info);

        // Object information to appear
        svg::object on;
        on._tag = "animate";
        on._attribute_map["fill"] = "freeze";
        on._attribute_map["attributeName"] = "display";
        on._attribute_map["from"] = "none";
        on._attribute_map["to"] = "block";
        on._attribute_map["begin"] = sid + __d + "mouseout";

        svg::object off;

        off._tag = "animate";
        off._attribute_map["fill"] = "freeze";
        off._attribute_map["attributeName"] = "display";
        off._attribute_map["to"] = "none";
        off._attribute_map["from"] = "block";
        off._attribute_map["begin"] = sid + __d + "mouseover";

        // Store the animation
        s_sheet_s._sub_objects.push_back(off);
        s_sheet_s._sub_objects.push_back(on);

        o_.add_object(s_sheet_s);
    }
}

}  // namespace display

}  // namespace actsvg
