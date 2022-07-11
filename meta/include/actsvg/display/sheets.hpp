// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <exception>
#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/display/helpers.hpp"
#include "actsvg/proto/cluster.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"

namespace actsvg {

using namespace defaults;

namespace display {

enum layer_type { e_endcap = 0, e_barrel = 1 };

enum sheet_type { e_module_info = 0, e_grid_info = 1 };

/** Make a surface sheet
 *
 * @tparam volume3_container is the type of the 3D point container
 *
 * @param id_ is the identifcation tag
 * @param s_ is the surface to be displayed
 * @param sh_ is the sheet size for displaying
 * @param kr_ is a directive to keep the ratio
 *
 * @return a surface sheet svg object
 **/
template <typename point3_container>
svg::object surface_sheet_xy(const std::string& id_,
                             const proto::surface<point3_container>& s_,
                             const std::array<scalar, 2>& sh_ = {400., 400.},
                             bool kr_ = true) noexcept(false) {
    svg::object so;
    so._tag = "g";
    so._id = id_;

    using point3 = typename point3_container::value_type;

    views::x_y x_y_view;

    // Add the surface
    std::vector<views::contour> contours;
    if (not s_._vertices.empty()) {
        auto surface_contour = x_y_view(s_._vertices);
        contours = {surface_contour};
    } else if (s_._radii[1] != 0.) {
        scalar ri = s_._radii[0];
        scalar ro = s_._radii[1];
        scalar phi_low = s_._opening[0];
        scalar phi_high = s_._opening[1];
        scalar phi = 0.5 * (phi_low + phi_high);
        scalar cos_phi_low = std::cos(phi_low);
        scalar sin_phi_low = std::sin(phi_low);
        scalar cos_phi_high = std::cos(phi_high);
        scalar sin_phi_high = std::cos(phi_high);
        point3 A = {ri * cos_phi_low, ri * sin_phi_low, 0.};
        point3 B = {ri * std::cos(phi), ri * std::sin(phi), 0.};
        point3 C = {ri * cos_phi_high, ri * sin_phi_high, 0.};
        point3 D = {ro * cos_phi_high, ro * sin_phi_high, 0.};
        point3 E = {ro * std::cos(phi), ro * std::sin(phi), 0.};
        point3 F = {ro * cos_phi_low, ro * sin_phi_low, 0.};
        std::vector<point3> vertices_disc = {A, B, C, D, E, F};
        auto surface_contour = x_y_view(vertices_disc);
        contours = {surface_contour};
    } else {
        throw std::invalid_argument(
            "surface_sheet_xy(...) - could not estimate range.");
    }

    auto [x_axis, y_axis] = display::view_range(contours);

    // Stitch when disc
    if (s_._type == proto::surface<point3_container>::e_disc) {
        x_axis[0] = 0.;
        y_axis[0] = 0.;
    }

    scalar s_x = sh_[0] / (x_axis[1] - x_axis[0]);
    scalar s_y = sh_[1] / (y_axis[1] - y_axis[0]);

    // Harmonize the view window
    if (kr_) {
        s_x = s_x < s_y ? s_x : s_y;
        s_y = s_y < s_x ? s_y : s_x;
    }

    // Create the scale transform
    style::transform scale_transform;
    scale_transform._scale = {s_x, s_y};

    // Copy in order to modify the transform
    proto::surface<point3_container> draw_surface = s_;
    draw_surface._transform._scale = {s_x, s_y};
    auto surface = display::surface(s_._name, draw_surface, x_y_view);
    so.add_object(surface);

    display::prepare_axes(x_axis, y_axis, s_x, s_y, 30., 30.);
    auto axis_font = __a_font;
    axis_font._size = 10;

    so.add_object(draw::x_y_axes(id_ + "_axes_xy", x_axis, y_axis, __a_stroke,
                                 "x", "y", axis_font));

    // The measures
    // - Trapezoid
    if (s_._type == proto::surface<point3_container>::e_trapez and
        s_._measures.size() == 3u) {

        scalar hlx_min = s_._measures[0] * s_x;
        scalar hlx_max = s_._measures[1] * s_x;
        scalar hly = s_._measures[2] * s_y;
        scalar ms = 2 * __m_marker._size;

        scalar hly_x = hlx_max + 2 * __m_marker._size;

        auto measure_hlx_min = draw::measure(
            "hlx_min", {0, -hly - ms}, {hlx_min, -hly - ms}, __m_stroke,
            __m_marker, "h_x_min = " + std::to_string(s_._measures[0]),
            __m_font, 1, -2);
        auto measure_hlx_max = draw::measure(
            "hlx_max", {0, hly + ms}, {hlx_max, hly + ms}, __m_stroke,
            __m_marker, "h_x_max = " + std::to_string(s_._measures[1]),
            __m_font, 1, 1);
        auto measure_hly = draw::measure(
            "hly", {hly_x, 0}, {hly_x, hly}, __m_stroke, __m_marker,
            "h_y = " + std::to_string(s_._measures[2]), __m_font, 1, 1);
        so.add_object(measure_hlx_min);
        so.add_object(measure_hlx_max);
        so.add_object(measure_hly);
    } else if (s_._type == proto::surface<point3_container>::e_rectangle and
               s_._measures.size() == 2u) {

        scalar hlx = s_._measures[0] * s_x;
        scalar hly = s_._measures[1] * s_y;
        scalar ms = 2 * __m_marker._size;

        auto measure_hlx = draw::measure(
            "hlx", {0, hly + ms}, {hlx, hly + ms}, __m_stroke, __m_marker,
            "h_x = " + std::to_string(s_._measures[0]), __m_font, 1, 1);
        auto measure_hly = draw::measure(
            "hly", {hlx + ms, 0}, {hlx + ms, hly}, __m_stroke, __m_marker,
            "h_y = " + std::to_string(s_._measures[1]), __m_font, 1, 1);
        so.add_object(measure_hlx);
        so.add_object(measure_hly);
    }

    return so;
}

/** Make a summary sheet for an endcap type volume
 *
 * @tparam volume3_container is the type of the 3D point container
 *
 * @param id_ is the identification tag of this sheet
 * @param v_ is the volume to be displayed
 * @param sh_ is the sheet size for displaying
 * @param t_ is the sheet type
 * @param s_sh_ is the surface sheet sub display size
 **/
template <typename point3_container, typename view_type, layer_type lT>
svg::object sheet(const std::string& id_,
                  const proto::volume<point3_container>& v_,
                  const std::array<scalar, 2>& sh_ = {600., 600.},
                  sheet_type t_ = e_module_info,
                  const std::array<scalar, 2>& s_sh_ = {200., 400.}) {

    svg::object eo;
    eo._tag = "g";
    eo._id = id_;

    view_type view;

    auto [modules, scale_transform, axes] = process_modules(v_, view, sh_);
    auto x_axis = axes[0];
    auto y_axis = axes[1];

    std::vector<svg::object> extra_objects;

    // Draw the template module surfaces
    if (t_ == e_module_info and not v_._template_surfaces.empty()) {
        // The templates
        std::vector<svg::object> templates;
        for (auto s : v_._template_surfaces) {
            // Use a local copy of the surface to modify color
            style::fill s_fill;
            s_fill._fc._rgb = s._fill._fc._hl_rgb;
            s_fill._fc._opacity = s._fill._fc._opacity;
            s._fill = s_fill;
            // The template sheet
            auto s_sheet =
                display::surface_sheet_xy(id_ + "_surface_sheet", s, s_sh_);
            style::transform(
                {{static_cast<scalar>(0.5 * sh_[0] + 0.5 * s_sh_[0] + 100), 0.,
                  0.}})
                .attach_attributes(s_sheet);
            templates.push_back(s_sheet);
        }
        // Connect the surface sheets
        if (v_._surfaces.size() == v_._templates.size()) {
            connect_surface_sheets(v_, templates, eo, 0.5 * sh_[1]);
        }
    } else if (t_ == e_grid_info and not v_._surface_grid._edges_0.empty() and
               not v_._surface_grid._edges_1.empty()) {

        // Draw the grid with the appropriate scale transform
        if (lT == e_endcap) {
            extra_objects =
                draw::tiled_polar_grid(id_, v_._surface_grid._edges_0,
                                       v_._surface_grid._edges_1, __g_fill,
                                       __g_stroke, scale_transform)
                    ._sub_objects;
        } else if (lT == e_barrel) {
            extra_objects =
                draw::tiled_cartesian_grid(id_, v_._surface_grid._edges_0,
                                           v_._surface_grid._edges_1, __g_fill,
                                           __g_stroke, scale_transform)
                    ._sub_objects;
        }
        // Connect grid and surfaces
        connectors::connect_objects(extra_objects, modules,
                                    v_._surface_grid._associations);
    }

    // Add the modules & eventual extra objects
    eo.add_objects(modules);
    eo.add_objects(extra_objects);

    // Add the axes on top
    auto axis_font = __a_font;
    axis_font._size = 10;
    eo.add_object(draw::x_y_axes("xy", x_axis, y_axis, __a_stroke,
                                 view._axis_names[0], view._axis_names[1],
                                 axis_font));

    //  Add the title text
    auto title_font = __t_font;
    title_font._size = 0.03 * sh_[0];
    auto title = draw::text("sheet_title",
                            {static_cast<scalar>(-0.55 * sh_[0]),
                             static_cast<scalar>(0.55 * sh_[1])},
                            {v_._name}, title_font);
    eo.add_object(title);

    return eo;
}

/** Make a summary sheet for an endcap type volume
 *
 * @tparam volume3_container is the type of the 3D point container
 *
 * @param id_ is the idenfication tag of this sheet
 * @param v_ is the volume to be displayed
 * @param sh_ is the sheet size for displaying
 * @param t_ is the sheet type
 * @param s_sh_ is the surface sheet sub display size
 **/
template <typename point3_container>
svg::object endcap_sheet(const std::string& id_,
                         const proto::volume<point3_container>& v_,
                         const std::array<scalar, 2>& sh_ = {600., 600.},
                         sheet_type t_ = e_module_info,
                         const std::array<scalar, 2>& s_sh_ = {200., 400.}) {

    return sheet<point3_container, views::x_y, e_endcap>(id_, v_, sh_, t_,
                                                         s_sh_);
}

/** Make a summary sheet for an barrel type volume
 *
 * @tparam volume3_container is the type of the 3D point container
 *
 * @param id_ is the idenfication tag of this sheet
 * @param v_ is the volume to be displayed
 * @param sh_ is the sheet size for displaying
 * @param t_ is the sheet type
 * @param s_sh_ is the surface sheet sub display size
 **/
template <typename point3_container>
svg::object barrel_sheet(const std::string& id_,
                         const proto::volume<point3_container>& v_,
                         const std::array<scalar, 2>& sh_ = {600., 600.},
                         sheet_type t_ = e_module_info,
                         const std::array<scalar, 2>& s_sh_ = {200., 400.}) {

    return sheet<point3_container, views::z_phi, e_barrel>(id_, v_, sh_, t_,
                                                           s_sh_);
}

}  // namespace display
}  // namespace actsvg
