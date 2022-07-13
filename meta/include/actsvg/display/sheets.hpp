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

    // Stitch when disc - and make x and y axis similarly long
    bool full = false;
    if (s_._type == proto::surface<point3_container>::e_disc or
        s_._type == proto::surface<point3_container>::e_annulus) {
        x_axis[0] = 0.;
        y_axis[0] = 0.;
        if (s_._opening == std::array<scalar, 2>({-M_PI, M_PI})) {
            y_axis[1] = x_axis[1];
            full = true;
        }
    }

    scalar s_x = sh_[0] / (x_axis[1] - x_axis[0]);
    scalar s_y = sh_[1] / (y_axis[1] - y_axis[0]);

    // Harmonize the view window with equal scales
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
        scalar ms = 2. * __m_marker._size;
        // Measure names
        std::string h_x_min = "h_x_min";
        std::string h_x_max = "h_x_max";
        std::string h_y = "h_y";

        scalar hly_x = hlx_max + 2 * __m_marker._size;

        auto measure_hlx_min = draw::measure(
            id_ + "_hlx_min", {0, -hly - ms}, {hlx_min, -hly - ms}, __m_stroke,
            __m_marker, __m_marker, __m_font,
            h_x_min + " = " + utils::to_string(s_._measures[0]),
            {static_cast<scalar>(0.5 * hlx_min),
             static_cast<scalar>((-hly - ms - 2.2 * __m_font._size))});
        auto measure_hlx_max = draw::measure(
            id_ + "_hlx_max", {0, hly + ms}, {hlx_max, hly + ms}, __m_stroke,
            __m_marker, __m_marker, __m_font,
            h_x_max + " = " + utils::to_string(s_._measures[1]),
            {static_cast<scalar>(0.5 * hlx_max),
             static_cast<scalar>(1.2 * __m_font._size + (hly + ms))});
        auto measure_hly =
            draw::measure(id_ + "_hly", {hly_x, 0}, {hly_x, hly}, __m_stroke,
                          __m_marker, __m_marker, __m_font,
                          h_y + " = " + utils::to_string(s_._measures[2]),
                          {static_cast<scalar>(__m_font._size + hly_x),
                           static_cast<scalar>(0.5 * hly)});
        so.add_object(measure_hlx_min);
        so.add_object(measure_hlx_max);
        so.add_object(measure_hly);
    } else if (s_._type == proto::surface<point3_container>::e_rectangle and
               s_._measures.size() == 2u) {

        scalar hlx = s_._measures[0] * s_x;
        scalar hly = s_._measures[1] * s_y;
        scalar ms = 2. * __m_marker._size;
        // Measure names
        std::string h_x = "h_x";
        std::string h_y = "h_y";

        auto measure_hlx = draw::measure(
            id_ + "_hlx", {0, hly + ms}, {hlx, hly + ms}, __m_stroke,
            __m_marker, __m_marker, __m_font,
            h_x + " = " + utils::to_string(s_._measures[0]),
            {static_cast<scalar>(0.5 * hlx),
             static_cast<scalar>((hly + ms + 1.2 * __m_font._size))});
        auto measure_hly =
            draw::measure(id_ + "_hly", {hlx + ms, 0}, {hlx + ms, hly},
                          __m_stroke, __m_marker, __m_marker, __m_font,
                          h_y + " = " + utils::to_string(s_._measures[1]),
                          {static_cast<scalar>(hlx + 1.2 * __m_font._size),
                           static_cast<scalar>(0.5 * hly)});
        so.add_object(measure_hlx);
        so.add_object(measure_hly);
    } else if (s_._type == proto::surface<point3_container>::e_polygon and
               s_._measures.size() == 2 * s_._vertices.size()) {

        point2 gcenter = {0., 0.};

        for (unsigned int iv = 0; iv < s_._vertices.size(); ++iv) {
            auto v = draw::marker(
                id_ + "_vertex_" + std::to_string(iv),
                {static_cast<scalar>(s_x * s_._measures[2 * iv]),
                 static_cast<scalar>(s_y * s_._measures[2 * iv + 1u])},
                style::marker({"o"}));

            gcenter[0] += s_x * s_._measures[2 * iv];
            gcenter[1] += s_y * s_._measures[2 * iv + 1u];

            so.add_object(v);
        }

        gcenter[0] /= s_._vertices.size();
        gcenter[1] /= s_._vertices.size();

        for (unsigned int iv = 0; iv < s_._vertices.size(); ++iv) {
            scalar x = s_x * s_._measures[2 * iv];
            scalar y = s_x * s_._measures[2 * iv + 1];

            scalar dx = x - gcenter[0];
            scalar dy = y - gcenter[1];

            scalar dnorm = std::sqrt(dx * dx + dy * dy);
            dx /= dnorm;
            dy /= dnorm;

            std::string label_v = "v" + std::to_string(iv) + " = ";
            label_v += utils::to_string(std::array<scalar, 2>{
                s_._measures[2 * iv], s_._measures[2 * iv + 1]});

            scalar offx = dx > 0 ? 2. * __m_font._size * dx
                                 : 0.5 * label_v.size() * __m_font._size * dx;
            scalar offy = 2 * __m_font._size * dy;

            so.add_object(draw::text(id_ + "_label_v" + std::to_string(iv),
                                     {x + offx, y + offy}, {label_v}));
        }

    } else if (s_._type == proto::surface<point3_container>::e_disc and
               not s_._measures.empty()) {

        std::string dphi = "h_phi";
        std::string aphi = "avg_phi";

        // Where to set the labels and how to label them
        std::vector<scalar> r_label;
        scalar phi_span = 2 * M_PI;

        if (full) {
            r_label = {M_PI * 0.25, -M_PI * 0.25};
        } else {
            phi_span = s_._opening[1] - s_._opening[0];
            r_label = {static_cast<scalar>(s_._opening[0] - 0.1),
                       static_cast<scalar>(s_._opening[1] + 0.1)};
        }
        // Start/end parameters of the labels
        scalar r_max = 0.;
        for (auto [ir, r] : utils::enumerate(s_._measures)) {

            std::array<scalar, 2> xs = {
                static_cast<scalar>(s_x * r * std::cos(r_label[0])),
                static_cast<scalar>(s_x * r * std::cos(r_label[1]))};
            std::array<scalar, 2> ys = {
                static_cast<scalar>(s_y * r * std::sin(r_label[0])),
                static_cast<scalar>(s_y * r * std::sin(r_label[1]))};

            // Radial labelling
            if (ir < 2 and
                std::abs(r) > std::numeric_limits<scalar>::epsilon()) {
                // Create a measurement helper
                if (not full) {
                    auto helper_r = draw::arc(id_ + "_arc_helper", s_x * r,
                                              {xs[0], ys[0]}, {xs[1], ys[1]},
                                              style::fill(), __m_stroke_guide);
                    so.add_object(helper_r);
                }

                // Record the maximum radius
                r_max = r > r_max ? r : r_max;
                std::string r_o = ir == 0 ? "r" : "R";
                auto measure_r = draw::measure(
                    id_ + "_r", {0., 0.}, {xs[ir], ys[ir]}, __m_stroke,
                    style::marker(), __m_marker, __m_font,
                    r_o + " = " + utils::to_string(r),
                    {static_cast<scalar>(__m_font._size + xs[ir]),
                     static_cast<scalar>(__m_font._size + ys[ir])});
                so.add_object(measure_r);
            }
            // Phi labelling
            if (ir == 2 and not full) {

                // Place it outside
                scalar lr = s_x * r_max + 2. * __m_marker._size;

                point2 start = {
                    static_cast<scalar>(lr * std::cos(s_._opening[0])),
                    static_cast<scalar>(lr * std::sin(s_._opening[0]))};

                // Medium phi arc and line
                scalar mphi = 0.5 * (s_._opening[0] + s_._opening[1]);
                point2 end = {static_cast<scalar>(lr * std::cos(mphi)),
                              static_cast<scalar>(lr * std::sin(mphi))};

                scalar mmphi = 0.5 * (s_._opening[0] + mphi);
                point2 mend = {
                    static_cast<scalar>(__m_font._size + lr * std::cos(mmphi)),
                    static_cast<scalar>(__m_font._size + lr * std::sin(mmphi))};

                auto maesure_arc = draw::arc_measure(
                    id_ + "_arc", lr, start, end, __m_stroke, __m_marker,
                    __m_marker, __m_font,
                    dphi + " = " + utils::to_string(phi_span), mend);
                so.add_object(maesure_arc);

                auto medium_phi_line = draw::line(id_ + "medium_phi", {0., 0.},
                                                  end, __m_stroke_guide);
                so.add_object(medium_phi_line);

                // Measure to the medium phi
                scalar r_avg_phi = 0.2 * s_x * r_max;
                std::array<scalar, 2> r_avg_start = {r_avg_phi, 0.};
                std::array<scalar, 2> r_avg_end = {
                    static_cast<scalar>(r_avg_phi * std::cos(mphi)),
                    static_cast<scalar>(r_avg_phi * std::sin(mphi))};

                std::array<scalar, 2> r_avg_mend = {
                    static_cast<scalar>(__m_font._size +
                                        r_avg_phi * std::cos(0.5 * mphi)),
                    static_cast<scalar>(__m_font._size +
                                        r_avg_phi * std::sin(0.5 * mphi))};

                auto measure_avg_phi = draw::arc_measure(
                    id_ + "_avg_phi", r_avg_phi, r_avg_start, r_avg_end,
                    __m_stroke, style::marker(), __m_marker, __m_font,
                    aphi + " = " + utils::to_string(mphi), r_avg_mend);
                so.add_object(measure_avg_phi);
            }
        }
    } else if (s_._type == proto::surface<point3_container>::e_annulus and
               not s_._measures.empty()) {

        if (s_._measures.size() != 7u) {
            throw std::invalid_argument(
                "surface_sheet_xy(...) - incorrect length of <measures> for "
                "annulus shape.");
        }

        // Special annulus bounds code
        scalar min_r = s_x * s_._measures[0];
        scalar max_r = s_x * s_._measures[1];
        scalar min_phi_rel = s_._measures[2];
        scalar max_phi_rel = s_._measures[3];
        // scalar average_phi = s_._measures[4];
        scalar origin_x = s_x * s_._measures[5];
        scalar origin_y = s_x * s_._measures[6];

        point2 cart_origin = {origin_x, origin_y};

        /// Find inner outer radius at edges in STRIP PC
        ///
        /// @note have a look at Acts/Surfaces/AnnulusBounds.hpp
        /// for more information
        ///
        auto circIx = [](scalar O_x, scalar O_y, scalar r,
                         scalar phi) -> point2 {
            //                      _____________________________________________
            //                     /      2  2                    2    2  2    2
            //     O_x + O_y*m - \/  - O_x *m  + 2*O_x*O_y*m - O_y  + m *r  + r
            // x =
            // --------------------------------------------------------------
            //                                  2
            //                                 m  + 1
            //
            // y = m*x
            //
            scalar m = std::tan(phi);
            point2 dir = {std::cos(phi), std::sin(phi)};
            scalar x1 =
                (O_x + O_y * m -
                 std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) +
                           2 * O_x * O_y * m - std::pow(O_y, 2) +
                           std::pow(m, 2) * std::pow(r, 2) + std::pow(r, 2))) /
                (std::pow(m, 2) + 1);
            scalar x2 =
                (O_x + O_y * m +
                 std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) +
                           2 * O_x * O_y * m - std::pow(O_y, 2) +
                           std::pow(m, 2) * std::pow(r, 2) + std::pow(r, 2))) /
                (std::pow(m, 2) + 1);

            point2 v1 = {x1, m * x1};
            if (v1[0] * dir[0] + v1[1] * dir[1] > 0) {
                return v1;
            }
            return {x2, m * x2};
        };

        auto out_left_s_xy = circIx(origin_x, origin_y, max_r, max_phi_rel);
        auto in_left_s_xy = circIx(origin_x, origin_y, min_r, max_phi_rel);
        auto out_right_s_xy = circIx(origin_x, origin_y, max_r, min_phi_rel);
        auto in_right_s_xy = circIx(origin_x, origin_y, min_r, min_phi_rel);

        std::vector<point2> corners = {in_right_s_xy, in_left_s_xy,
                                       out_right_s_xy, out_left_s_xy};

        so.add_object(draw::line(id_ + "_phi_line_l", {0, 0}, in_left_s_xy,
                                 __m_stroke_guide));
        so.add_object(draw::line(id_ + "_phi_line_r", {0, 0}, in_right_s_xy,
                                 __m_stroke_guide));

        std::vector<std::string> phi_labels = {"phi_min_rel", "phi_max_rel"};
        // Measure the phi values
        for (unsigned int ic = 0; ic < 2; ++ic) {

            const auto& corner = corners[ic];
            scalar in_phi = atan2(corner[1], corner[0]);

            scalar in_r = (0.2 + ic * 0.2) * utils::perp(corners[ic]);
            std::array<scalar, 2> in_start = {in_r, 0.};
            std::array<scalar, 2> in_end = {
                static_cast<scalar>(in_r * std::cos(in_phi)),
                static_cast<scalar>(in_r * std::sin(in_phi))};

            std::array<scalar, 2> in_mend = {
                static_cast<scalar>(__m_font._size +
                                    in_r * std::cos(0.5 * in_phi)),
                static_cast<scalar>(__m_font._size +
                                    in_r * std::sin(0.5 * in_phi))};

            auto in_phi_arc = draw::arc_measure(
                id_ + "_in_phi_" + std::to_string(ic), in_r, in_start, in_end,
                __m_stroke, style::marker(), __m_marker, __m_font,
                phi_labels[ic] + " = " + utils::to_string(in_phi), in_mend);
            so.add_object(in_phi_arc);
        }

        style::marker cart_origin_marker = style::marker({"o"});
        cart_origin_marker._fill = style::fill(style::color{{255, 0, 0}});

        style::font cart_font = __m_font;
        cart_font._fc = style::color{{255, 0, 0}};

        style::stroke cart_stroke = __m_stroke;
        cart_stroke._sc = style::color{{255, 0, 0}};

        style::marker cart_marker = __m_marker;
        cart_marker._fill = style::fill(style::color{{255, 0, 0}});
        cart_marker._stroke = cart_stroke;

        // The origin of the cartesian system
        so.add_object(draw::marker(id_ + "_origin_cart", {origin_x, origin_y},
                                   cart_origin_marker));

        so.add_object(draw::text(
            id_ + "_origin_label",
            {origin_x - 2 * cart_font._size, origin_y - 2 * cart_font._size},
            {"cart_origin = " + utils::to_string(std::array<scalar, 2>{
                                    s_._measures[5], s_._measures[6]})},
            cart_font));

        style::stroke cart_guide = __m_stroke_guide;
        cart_guide._sc = style::color{{255, 0, 0}};

        std::vector<std::string> rs = {"r", "R"};
        // Draw the circles
        for (unsigned int idc = 0; idc < 2; ++idc) {

            const auto& c0 = corners[2 * idc];
            const auto& c1 = corners[2 * idc + 1];

            point2 cart_c0 = point2{c0[0] - origin_x, c0[1] - origin_y};
            point2 cart_c1 = point2{c1[0] - origin_x, c1[1] - origin_y};
            scalar cart_r = utils::perp(cart_c0);

            scalar cart_phi0 = atan2(cart_c0[1], cart_c0[0]) - 0.2;
            scalar cart_phi1 = atan2(cart_c1[1], cart_c1[0]) + 0.25;
            scalar cart_phir = cart_phi1 - (2 * idc + 1) * 0.05;

            point2 start_arc = {
                static_cast<scalar>(origin_x + cart_r * std::cos(cart_phi0)),
                static_cast<scalar>(origin_y + cart_r * std::sin(cart_phi0))};
            point2 end_arc = {
                static_cast<scalar>(origin_x + cart_r * std::cos(cart_phi1)),
                static_cast<scalar>(origin_y + cart_r * std::sin(cart_phi1))};

            point2 end_r = {
                static_cast<scalar>(origin_x + cart_r * std::cos(cart_phir)),
                static_cast<scalar>(origin_y + cart_r * std::sin(cart_phir))};

            auto helper_r =
                draw::arc(id_ + "_arc_helper_" + std::to_string(idc), cart_r,
                          start_arc, end_arc, style::fill(), cart_guide);
            so.add_object(helper_r);

            auto measure_r = draw::measure(
                id_ + "_measure_r_" + std::to_string(idc), cart_origin, end_r,
                cart_stroke, style::marker(), cart_marker, cart_font,
                rs[idc] + " = " +
                    utils::to_string(static_cast<scalar>(cart_r / s_x)),
                {static_cast<scalar>(__m_font._size + end_r[0]),
                 static_cast<scalar>(__m_font._size + end_r[1])});
            so.add_object(measure_r);
        }
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
