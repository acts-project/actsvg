// This file is part of the actsvg packge.
//
// Copyright (C) 2022-2024 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"
#include "actsvg/styles/defaults.hpp"

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
void prepare_axes(std::array<scalar, 2>& first_, std::array<scalar, 2>& second_,
                  scalar sx_, scalar sy_, scalar ax_ = 0., scalar ay_ = 0.);

/** Helper method to get the contours
 *
 * @param s_ is the surface
 * @param fs_ draw in focus mode
 *
 * @return a vector of contours;
 **/
template <typename surface_type, typename view_type = views::x_y>
views::contour range_contour(const surface_type& s_,
                             bool fs_ = false) noexcept(false) {

    view_type view;

    using point3 = typename surface_type::point3_type;

    // Add the surface
    views::contour contour;
    if (s_._template_object.is_defined()) {
        // Use a bounding box trick
        // Get the contour from the template as bounding box
        point2 bbl = {s_._template_object._x_range[0],
                      s_._template_object._y_range[0]};
        point2 bbr = {s_._template_object._x_range[1],
                      s_._template_object._y_range[1]};
        if (not fs_) {
            bbl[0] += s_._transform._tr[0];
            bbl[1] += s_._transform._tr[1];
            bbr[0] += s_._transform._tr[0];
            bbr[1] += s_._transform._tr[1];
            scalar alpha = s_._transform._rot[0];
            if (alpha != 0.) {
                scalar alpha_rad = static_cast<scalar>(alpha / M_PI * 180.);
                bbl = utils::rotate(bbl, alpha_rad);
                bbr = utils::rotate(bbr, alpha_rad);
            }
        }
        contour = {bbl, bbr};
    } else if (not s_._vertices.empty()) {
        // Plain contours are present
        contour = view.path(s_._vertices);
        if (not fs_) {
            std::for_each(contour.begin(), contour.end(), [&](auto& v) {
                scalar alpha = s_._transform._rot[0];
                scalar alpha_rad = static_cast<scalar>(alpha / M_PI * 180.);
                v = utils::rotate(v, alpha_rad);
                v[0] += s_._transform._tr[0];
                v[1] += s_._transform._tr[1];
            });
        }
    } else if (s_._radii[1] != 0.) {
        // Create a contour
        scalar ri = s_._radii[0];
        scalar ro = s_._radii[1];
        scalar phi_low = s_._opening[0];
        scalar phi_high = s_._opening[1];
        scalar phi = 0.5_scalar * (phi_low + phi_high);
        scalar cos_phi_low = std::cos(phi_low);
        scalar sin_phi_low = std::sin(phi_low);
        scalar cos_phi_high = std::cos(phi_high);
        scalar sin_phi_high = std::cos(phi_high);
        /// @todo add transform
        point3 A = {ri * cos_phi_low, ri * sin_phi_low, 0.};
        point3 B = {ri * std::cos(phi), ri * std::sin(phi), 0.};
        point3 C = {ri * cos_phi_high, ri * sin_phi_high, 0.};
        point3 D = {ro * cos_phi_high, ro * sin_phi_high, 0.};
        point3 E = {ro * std::cos(phi), ro * std::sin(phi), 0.};
        point3 F = {ro * cos_phi_low, ro * sin_phi_low, 0.};

        std::vector<point3> vertices_disc = {A, B, C, D, E, F};
        contour = view.path(vertices_disc);
    } else {
        throw std::invalid_argument(
            "surface_sheet_xy(...) - could not estimate range.");
    }
    return contour;
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
 * @param es_ is the equal scale
 * @param hl_ switch highlighting on/off
 * @param dt_ draw template is available
 *
 * @returns the modules (per surface batch), a scale transform & the axes
 **/
template <typename volume_type, typename view_type>
std::tuple<std::vector<std::vector<svg::object>>, style::transform,
           std::array<std::array<scalar, 2>, 2>>
process_modules(const volume_type& v_, const view_type& view_,
                const std::array<scalar, 2>& sh_ = {600., 600.},
                bool es_ = false, bool hl_ = true, bool dt_ = true) {

    using surface_type = typename volume_type::surface_type;

    // Axis range & pre-loop, create contours
    std::vector<std::vector<views::contour>> all_contours;
    for (const auto& surfaces : v_._surfaces) {
        std::vector<views::contour> contours;
        contours.reserve(surfaces.size());
        for (const auto& s : surfaces) {
            surface_type range_surface = s;
            if (not dt_) {
                range_surface._template_object = svg::object{};
            }
            contours.push_back(
                range_contour<surface_type, view_type>(range_surface));
        }
        all_contours.push_back(contours);
    }

    // Concatenate all contours into one for the range display
    std::vector<views::contour> flattened_contours;
    for (const auto& c : all_contours) {
        flattened_contours.insert(flattened_contours.end(), c.begin(), c.end());
    }

    // Get the scaling right - same scaling for all potential views
    auto axes = display::view_range(flattened_contours);

    scalar s_x = sh_[0] / (axes[0][1] - axes[0][0]);
    scalar s_y = sh_[1] / (axes[1][1] - axes[1][0]);

    if (es_) {
        // Harmonize the view window with equal scales
        s_x = s_x < s_y ? s_x : s_y;
        s_y = s_y < s_x ? s_y : s_x;
    }

    // Create the scale transform
    style::transform scale_transform;
    scale_transform._scale = {s_x, s_y};

    // Draw the modules and estimate axis ranges
    std::vector<std::vector<svg::object>> all_modules;
    for (auto [ics, contours] : utils::enumerate(all_contours)) {
        std::vector<svg::object> modules;
        modules.reserve(contours.size());
        for (auto [ic, c] : utils::enumerate(contours)) {
            surface_type draw_surface = v_._surfaces[ics][ic];
            draw_surface._transform._scale = {s_x, s_y};
            if (not hl_) {
                draw_surface._fill._fc._highlight = {};
            }
            if (not dt_) {
                draw_surface._template_object = svg::object{};
            }

            auto surface_module =
                display::surface(draw_surface._name, draw_surface, view_,
                                 {true, false, true, false});
            modules.push_back(surface_module);
        }
        all_modules.push_back(modules);
    }

    // Prepare the axis for the view range
    prepare_axes(axes[0], axes[1], s_x, s_y, 30., 30.);

    return {all_modules, scale_transform, axes};
}

/** Helper method to connect the surface sheets to the
 * surfaces of the layer_sheets
 *
 * @tparam volume_type the type of volume (templated on point3_container)
 *
 * @param v_ the input volume
 * @param templates_ the given module templtes
 * @param o_ the object to which they are attached
 *
 **/
template <typename volume_type>
void connect_surface_sheets(const volume_type& v_,
                            std::vector<std::vector<svg::object>>& templates_,
                            svg::object& o_) {
    // Now create an item per surface
    for (auto [ib, surface_batch] : utils::enumerate(v_._surfaces)) {
        for (auto [is, s] : utils::enumerate(surface_batch)) {
            std::string sid = s._name;

            svg::object& s_sheet_s = templates_[ib][is];
            s_sheet_s._attribute_map["display"] = "none";

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
}

}  // namespace display

}  // namespace actsvg
