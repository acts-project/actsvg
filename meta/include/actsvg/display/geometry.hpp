// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/proto/surface.hpp"

namespace actsvg {

using namespace defaults;

namespace display {

/** Draw the surface with a dedicated view
 *
 * @param id_ the identification of this surface
 * @param s_ the surface type
 * @param v_ the view type
 * @param b_ draw the boolean
 * @param fs_ draw as focus
 * @param sc_ draw at scale
 * @param dt_ draw as template
 *
 * @note template surfaces ignore the view_type::scene range restriction
 */
template <typename surface_type, typename view_type>
svg::object surface(const std::string& id_, const surface_type& s_,
                    const view_type& v_, bool _b = true, bool fs_ = false,
                    bool sc_ = false, bool dt_ = false) {

    svg::object s;

    // If the surface has a template and is it defined
    if (s_._template_object.is_defined()) {

        style::transform draw_transform = s_._transform;
        // No rotation nor shift as template
        if (dt_) {
            draw_transform._tr = {0., 0.};
            draw_transform._rot = {0., 0., 0.};
        }
        // Apply scale or not
        if (not sc_) {
            draw_transform._scale = {1., 1.};
        }

        // Create a surface object from the template
        s = draw::from_template(id_, s_._template_object, s_._fill, s_._stroke,
                                draw_transform);
        return s;
    }

    style::transform draw_transform = fs_ ? style::transform{} : s_._transform;
    draw_transform._scale = s_._transform._scale;

    // Surface directly
    if (s_._type == surface_type::type::e_disc) {

        // x-y view for discs
        if constexpr (std::is_same_v<view_type, views::x_y>) {

            // view activation / deactivation
            if (not utils::inside_range(s_._zparameters[0u],
                                        v_._scene._range[1u])) {
                s._active = false;
                return s;
            }
            // Sector or circle
            if (std::abs((s_._opening[1u] - s_._opening[0u]) - 2 * M_PI) >
                5 * std::numeric_limits<scalar>::epsilon()) {
                auto view_vertices = generators::sector_contour(
                    s_._radii[0u], s_._radii[1u], s_._opening[0u],
                    s_._opening[1u]);
                s = draw::polygon(id_, view_vertices, s_._fill, s_._stroke,
                                  draw_transform);
            } else {

                s = draw::circle(id_, {0., 0.}, s_._radii[1u], s_._fill,
                                 s_._stroke, draw_transform);

                // A ring is present
                if (s_._radii[0u] != 0.) {

                    std::string mask_id = id_ + "_mask";

                    auto s_c_ = s_;
                    s_c_._radii = {0., s_._radii[1u]};

                    svg::object outer_mask =
                        surface(id_ + "_mask_surface_outer", s_c_, v_, false);
                    outer_mask._fill = style::fill{true};
                    outer_mask._stroke = style::stroke{true};
                    outer_mask._attribute_map["fill"] = "white";

                    s_c_._radii = {0., s_._radii[0u]};
                    svg::object inner_mask =
                        surface(id_ + "_mask_surface_inner", s_c_, v_, false);
                    inner_mask._fill = style::fill{true};
                    inner_mask._stroke = style::stroke{true};
                    inner_mask._attribute_map["fill"] = "black";

                    // Create the mask object
                    svg::object mask;
                    mask._fill = style::fill{true};
                    mask._stroke = style::stroke{true};
                    mask._id = mask_id;
                    mask._tag = "mask";
                    mask.add_object(outer_mask);
                    mask.add_object(inner_mask);

                    // Modify the surface
                    s._definitions.push_back(mask);
                    s._attribute_map["mask"] = utils::id_to_url(mask_id);
                }
            }
        }
        // r_z view for discs
        if constexpr (std::is_same_v<view_type, views::z_r>) {
            scalar zpos = s_._zparameters[0u];
            point2 start = {zpos, s_._radii[0u]};
            point2 end = {zpos, s_._radii[1u]};
            s = draw::line(id_, start, end, s_._stroke, draw_transform);
        }

    } else if (s_._type == surface_type::type::e_cylinder) {

        // xy - view
        if constexpr (std::is_same_v<view_type, views::x_y>) {
            // view activation / deactivation
            if (not utils::inside_range(s_._zparameters[0u],
                                        v_._scene._range[1u])) {
                s._active = false;
                return s;
            }

            if (std::abs((s_._opening[1u] - s_._opening[0u]) - 2 * M_PI) >
                5 * std::numeric_limits<scalar>::epsilon()) {
                scalar r = s_._radii[1u];
                point2 start = {r * std::cos(s_._opening[0u]),
                                r * s_._opening[0u]};
                point2 end = {r * std::cos(s_._opening[1u]),
                              r * s_._opening[1u]};
                s = draw::arc(id_, r, start, end, style::fill({}), s_._stroke,
                              draw_transform);
            } else {
                s = draw::circle(id_, {0., 0.}, s_._radii[1u], __w_fill,
                                 s_._stroke, draw_transform);
            }
        }

        // z-r view
        if constexpr (std::is_same_v<view_type, views::z_r>) {
            scalar zpos = s_._zparameters[0u];
            scalar zhalf = s_._zparameters[1u];
            point2 start = {zpos - zhalf, s_._radii[1u]};
            point2 end = {zpos + zhalf, s_._radii[1u]};
            s = draw::line(id_, start, end, s_._stroke, draw_transform);
        }

    } else if (s_._type == surface_type::type::e_straw) {

        // xy view
        if constexpr (std::is_same_v<view_type, views::x_y>) {
            // Skin of the straw
            auto ss = draw::circle(id_, {0., 0.}, s_._radii[1u], s_._fill,
                                   s_._stroke, draw_transform);
            // Wire of the straw
            auto ws = draw::circle(id_, {0., 0.}, s_._radii[0u], __s_fill,
                                   s_._stroke, draw_transform);
            s._tag = "g";
            s.add_object(ss);
            s.add_object(ws);
        }

    } else if (not s_._vertices.empty()) {

        // View activiation, deactivation
        // if only one vertex is within the view range, the surface is shown
        bool view_active = true;
        if constexpr (std::is_same_v<view_type, views::x_y>) {
            view_active = false;
            for (auto v : s_._vertices) {
                if (utils::inside_range(v[2], v_._scene._range[1u])) {
                    view_active = true;
                    break;
                }
            }
        }
        if (not view_active) {
            s._active = view_active;
            return s;
        }
        auto view_vertices = v_.path(s_._vertices);

        if constexpr (std::is_same_v<view_type, views::z_rphi>) {
            // Check if we have to split the surface in phi and
            // draw wiggle
            // - currently supported for rectangular surfaces only
            if (s_._vertices.size() == 4u) {
                scalar min_phi = std::numeric_limits<scalar>::max();
                scalar max_phi = std::numeric_limits<scalar>::min();
                // Pre-emptively split the vertices
                std::vector<typename surface_type::point3_type> n_vertices;
                std::vector<typename surface_type::point3_type> p_vertices;
                // Calculate phi
                for (const auto& v : s_._vertices) {
                    scalar phi = std::atan2(v[1u], v[0u]);
                    min_phi = std::min(min_phi, phi);
                    max_phi = std::max(max_phi, phi);
                    if (phi < 0.) {
                        n_vertices.push_back(v);
                    } else {
                        p_vertices.push_back(v);
                    }
                }
                // Detect phi boundary jump
                if (max_phi - min_phi > 1.5 * M_PI) {
                    // Check first if you have only collected one
                    // negative vertex. It may happen if you have
                    // tilted surfaces. Handle with care and
                    // re-evaluate the intersections with
                    // the grid to show split surfaces correcly
                    if (n_vertices.size() < p_vertices.size()) {
                        std::vector<typename surface_type::point3_type>
                            bottom_vertices = {p_vertices[2u], p_vertices[0u]};
                        std::vector<typename surface_type::point3_type>
                            top_vertices = {n_vertices[0u], n_vertices[0u]};

                        for (size_t index = 0u; index < 2u; index++) {
                            auto p_bottom = bottom_vertices[index];
                            auto p_top = top_vertices[index];
                            float m = (p_bottom[1u] - p_top[1u]) /
                                      (p_bottom[2u] - p_top[2u]);
                            float q = p_bottom[1u] - m * p_bottom[2u];
                            float z = -q / m;

                            typename surface_type::point3_type add = {
                                p_bottom[0u],
                                std::numeric_limits<scalar>::epsilon(), z};
                            p_vertices.push_back(add);
                            add[1u] = -std::numeric_limits<scalar>::epsilon();
                            n_vertices.push_back(add);
                        }

                        auto p_view_vertices = v_.path(p_vertices);
                        auto p_middle = p_view_vertices.back();
                        p_middle[1u] = static_cast<scalar>(
                            1.5 * p_view_vertices.back()[1u] -
                            0.5 * p_view_vertices.front()[1u]);
                        p_middle[0u] = static_cast<scalar>(
                            0.5 * (p_view_vertices.at(p_view_vertices.size() -
                                                      2)[0u] +
                                   p_view_vertices.back()[0u]));
                        p_view_vertices.insert(p_view_vertices.end() - 1,
                                               p_middle);

                        auto s_n = draw::polygon(id_ + std::string("_n_split"),
                                                 v_.path(n_vertices), s_._fill,
                                                 s_._stroke, draw_transform);

                        auto s_p = draw::polygon(id_ + std::string("_p_split"),
                                                 p_view_vertices, s_._fill,
                                                 s_._stroke, draw_transform);

                        // The surface turns into a group of two
                        s._tag = "g";

                        s.add_object(s_n);
                        s.add_object(s_p);
                        return s;
                    }

                    // The surface turns into a group of two
                    s._tag = "g";

                    // Split part on negative side
                    auto n4 = n_vertices[0u];
                    auto n2 = n_vertices[1u];

                    n2[1u] = -std::numeric_limits<scalar>::epsilon();
                    n4[1u] = -std::numeric_limits<scalar>::epsilon();

                    // Create the wiggle in
                    typename surface_type::point3_type n3 = {
                        n2[0u],
                        static_cast<scalar>(0.5 *
                                            (n_vertices[1u][1u] + n2[1u])),
                        static_cast<scalar>(0.5 * (n2[2u] + n4[2u]))};

                    n_vertices.push_back(n2);
                    n_vertices.push_back(n3);
                    n_vertices.push_back(n4);

                    auto s_n = draw::polygon(id_ + std::string("_n_split"),
                                             v_.path(n_vertices), s_._fill,
                                             s_._stroke, draw_transform);
                    // Split art on positive side
                    auto p4 = p_vertices[0u];
                    auto p2 = p_vertices[1u];
                    p2[1u] = std::numeric_limits<scalar>::epsilon();
                    p4[1u] = std::numeric_limits<scalar>::epsilon();
                    typename surface_type::point3_type p3 = {
                        p2[0u], p2[1u],
                        static_cast<scalar>(0.5 * (p2[2u] + p4[2u]))};
                    p_vertices.push_back(p2);
                    p_vertices.push_back(p3);
                    p_vertices.push_back(p4);
                    // Create the wiggle out (on view to not fall into phi trap)
                    auto p_view_vertices = v_.path(p_vertices);
                    const auto& p_v_1 = p_view_vertices[1u];
                    const auto& p_v_2 = p_view_vertices[2u];
                    auto& p_v_3 = p_view_vertices[3u];
                    p_v_3[1u] =
                        static_cast<scalar>(1.5 * p_v_2[1u] - 0.5 * p_v_1[1u]);

                    auto s_p = draw::polygon(id_ + std::string("_p_split"),
                                             p_view_vertices, s_._fill,
                                             s_._stroke, draw_transform);
                    s.add_object(s_n);
                    s.add_object(s_p);
                    return s;
                }
            }
        }

        s = draw::polygon(id_, view_vertices, s_._fill, s_._stroke,
                          draw_transform);
    }

    if (_b) {
        /// Boolean surfaces only supported for x-y view so far
        if constexpr (std::is_same_v<view_type, views::x_y>) {
            if (s_._boolean_surface.size() == 1u and
                s_._boolean_operation == surface_type::boolean::e_subtraction) {
                std::string mask_id = id_ + "_mask";
                // make a new boolean surface
                svg::object outer_mask =
                    surface(id_ + "_mask_surface_outer", s_, v_, false);
                outer_mask._fill = style::fill{true};
                outer_mask._stroke = style::stroke{true};
                outer_mask._attribute_map["fill"] = "white";

                svg::object inner_mask = surface(id_ + "_mask_surface_inner",
                                                 s_._boolean_surface[0], v_);
                inner_mask._fill = style::fill{true};
                inner_mask._stroke = style::stroke{true};
                inner_mask._attribute_map["fill"] = "black";

                // Create the mask object
                svg::object mask;
                mask._fill = style::fill{true};
                mask._stroke = s_._stroke;
                mask._id = mask_id;
                mask._tag = "mask";
                mask.add_object(outer_mask);
                mask.add_object(inner_mask);

                // Modify the surface
                s._definitions.push_back(mask);
                s._attribute_map["mask"] = utils::id_to_url(mask_id);
            }
        }
    }

    return s;
}

/** Draw a portal link
 *
 * @param id_ the indentification of this portal link
 * @param p_ the portal for understanding the span
 * @param link_ the link itself
 * @param v_ the view type
 *
 * @return a single object containing the portal view
 **/
template <typename portal_type, typename view_type>
svg::object portal_link(const std::string& id_,
                        [[maybe_unused]] const portal_type& p_,
                        const typename portal_type::link& link_,
                        const view_type& v_) {
    svg::object l;
    l._tag = "g";
    l._id = id_;

    scalar d_z = static_cast<scalar>(link_._end[2u] - link_._start[2u]);
    // View activation / deactivation
    if constexpr (std::is_same_v<view_type, views::x_y>) {
        if (v_._scene._strict and d_z * v_._scene._view[1] < 0) {
            l._active = false;
            return l;
        }
    }

    scalar d_x = static_cast<scalar>(link_._end[0u] - link_._start[0u]);
    scalar d_y = static_cast<scalar>(link_._end[1u] - link_._start[1u]);
    scalar d_r = std::sqrt(d_x * d_x + d_y * d_y);

    if (std::is_same_v<view_type, views::x_y> and
        d_r <= std::numeric_limits<scalar>::epsilon()) {
        svg::object arr_xy;
        arr_xy._tag = "g";
        arr_xy._id = id_ + "_arrow";
        arr_xy.add_object(draw::circle(id_ + "_arrow_top",
                                       {static_cast<scalar>(link_._start[0u]),
                                        static_cast<scalar>(link_._start[1u])},
                                       link_._end_marker._size,
                                       link_._end_marker._fill));
        // Camera view onto the surface
        if (d_z * v_._scene._view[1] < 0) {
            arr_xy.add_object(
                draw::circle(id_ + "_arrow_top_tip",
                             {static_cast<scalar>(link_._start[0u]),
                              static_cast<scalar>(link_._start[1u])},
                             link_._end_marker._size * 0.1, __w_fill));
        } else {
            scalar d_l_x =
                link_._end_marker._size * 0.9 * std::cos(0.25 * M_PI);
            scalar d_l_y =
                link_._end_marker._size * 0.9 * std::sin(0.25 * M_PI);
            arr_xy.add_object(
                draw::line(id_ + "_arrow_top_cl0",
                           {static_cast<scalar>(link_._start[0u] - d_l_x),
                            static_cast<scalar>(link_._start[1u] - d_l_y)},
                           {static_cast<scalar>(link_._start[0u] + d_l_x),
                            static_cast<scalar>(link_._start[1u] + d_l_y)},
                           __w_stroke));
            arr_xy.add_object(
                draw::line(id_ + "_arrow_top_cl1",
                           {static_cast<scalar>(link_._start[0u] + d_l_x),
                            static_cast<scalar>(link_._start[1u] - d_l_y)},
                           {static_cast<scalar>(link_._start[0u] - d_l_x),
                            static_cast<scalar>(link_._start[1u] + d_l_y)},
                           __w_stroke));
        }
        l.add_object(arr_xy);
        // draw plot
    } else {
        typename portal_type::container_type start_end_3d = {link_._start,
                                                             link_._end};
        auto start_end = v_.path(start_end_3d);

        l.add_object(draw::arrow(id_ + "_arrow", start_end[0u], start_end[1u],
                                 link_._stroke, link_._start_marker,
                                 link_._end_marker));
    }
    return l;
}

/** Draw a portal with a dedicated view
 *
 * @param id_ the identification of this surface
 * @param s_ the surface type
 * @param v_ the view type
 *
 * @return a single object containing the portal view;
 **/
template <typename portal_type, typename view_type>
svg::object portal(const std::string& id_, const portal_type& p_,
                   const view_type& v_) {
    svg::object p;
    p._tag = "g";
    p._id = id_;
    p._fill._sterile = true;
    p._stroke._sterile = true;

    p.add_object(surface(id_ + "_surface", p_._surface, v_));
    for (auto [il, vl] : utils::enumerate(p_._volume_links)) {
        p.add_object(portal_link(id_ + "_volume_link_" + std::to_string(il), p_,
                                 vl, v_));
    }

    return p;
}

/** Draw a volume
 *
 * @param id_ the indentification of this portal link
 * @param dv_ the detector volume
 * @param v_ the view type
 * @param p_ draw the portals
 * @param s_ draw the surfaces
 *
 * @return a single object containing the volume view
 **/
template <typename volume_type, typename view_type>
svg::object volume(const std::string& id_, const volume_type& dv_,
                   const view_type& v_, bool p_ = true, bool s_ = true) {
    svg::object v;
    v._tag = "g";
    v._id = id_;
    v._fill._sterile = true;
    v._stroke._sterile = true;

    // The volume shape - only rz view currently supported
    if (not dv_._vertices.empty() and std::is_same_v<view_type, views::z_r>) {
        auto view_vertices = v_.path(dv_._vertices);
        auto pv = draw::polygon(id_ + "_volume", view_vertices, dv_._fill,
                                dv_._stroke, dv_._transform);
        v.add_object(pv);
    } else {
        if (dv_._type == volume_type::type::e_cylinder and
            dv_._bound_values.size() >= 6u) {
            scalar ri = dv_._bound_values[0u];
            scalar ro = dv_._bound_values[1u];
            scalar zp = dv_._bound_values[2u];
            scalar zh = dv_._bound_values[3u];
            scalar ps = dv_._bound_values[4u];
            scalar ap = dv_._bound_values[5u];
            if constexpr (std::is_same_v<view_type, views::x_y>) {
                // Make a dummy surface and draw it
                typename volume_type::surface_type s;
                s._name = id_ + "_volume";
                s._radii = {ri, ro};
                s._opening = {ap - ps, ap + ps};
                s._zparameters = {zp, zh};
                s._fill = dv_._fill;
                s._stroke = dv_._stroke;
                s._type = decltype(s)::type::e_disc;
                v.add_object(surface(s._name, s, v_));
            }
            if constexpr (std::is_same_v<view_type, views::z_r>) {
                std::vector<point2> view_vertices = {
                    {zp - zh, ri}, {zp + zh, ri}, {zp + zh, ro}, {zp - zh, ro}};
                auto pv = draw::polygon(id_ + "_volume", view_vertices,
                                        dv_._fill, dv_._stroke, dv_._transform);
                v.add_object(pv);
            }
        }
    }

    // Draw the portals
    if (p_) {
        for (auto [ip, p] : utils::enumerate(dv_._portals)) {
            v.add_object(portal(id_ + "_portal_" + std::to_string(ip), p, v_));
        }
    }
    // Draw the surfaces
    if (s_) {
        for (auto [is, s] : utils::enumerate(dv_._v_surfaces)) {
            v.add_object(display::surface(
                id_ + "_surface_" + std::to_string(is), s, v_));
        }
    }
    return v;
}

/** Draw a detector
 *
 * @param id_ the indentification of this portal link
 * @param d_ the detector
 * @param v_ the view type
 *
 * @return a single object containing the volume view
 **/
template <typename detector_type, typename view_type>
svg::object detector(const std::string& id_, const detector_type& d_,
                     const view_type& v_) {
    svg::object d;
    d._tag = "g";
    d._id = id_;
    d._fill._sterile = true;
    d._stroke._sterile = true;

    // Sort the volumes after their depth level, local copy first
    auto volumes = d_._volumes;
    std::sort(volumes.begin(), volumes.end(),
              [](const auto& a_, const auto& b_) {
                  return (a_._depth_level < b_._depth_level);
              });

    // Draw the volume areas first
    for (const auto& v : volumes) {
        d.add_object(volume(v._name, v, v_, false));
    }

    // Collect all the portals (in a named map) - to avoid double drawing of
    // shared ones
    std::map<std::string, typename detector_type::volume_type::portal_type>
        portals;

    for (const auto& v : volumes) {
        for (const auto& p : v._portals) {
            portals[p._name] = p;
        }
    }
    // Now draw the portals
    for (const auto& [n, p] : portals) {
        d.add_object(portal(n, p, v_));
    }

    return d;
}

/** Draw eta lines in a zr view
 *
 * @param id_ the identiier
 * @param zr_ the z range of the detector
 * @param rr_ the r range of the detector
 * @param els_ the stroked eta lines + boolean whether to label
 * @param tr_ a potential transform
 *
 * @return a single object containing the frame
 */
static inline svg::object eta_lines(
    const std::string& id_, scalar zr_, scalar rr_,
    const std::vector<std::tuple<std::vector<scalar>, style::stroke, bool,
                                 style::font>>& els_,
    const style::transform& tr_ = style::transform()) {

    svg::object e;
    e._tag = "g";
    e._id = id_;
    e._transform = tr_;

    auto theta_from_eta = [](scalar eta) -> scalar {
        return static_cast<scalar>(2. * std::atan(std::exp(-eta)));
    };

    scalar theta_cut = std::atan2(rr_, zr_);

    for (auto [iet, elt] : utils::enumerate(els_)) {
        auto stroke = std::get<style::stroke>(elt);
        for (auto [ie, eta] :
             utils::enumerate(std::get<std::vector<scalar>>(elt))) {
            scalar theta = theta_from_eta(eta);
            std::array<scalar, 2> start = {0., 0.};
            std::array<scalar, 2> end;
            if (theta < theta_cut) {
                end = {zr_, static_cast<scalar>(zr_ * std::tan(theta))};
            } else {
                end = {static_cast<scalar>(rr_ * 1 / std::tan(theta)), rr_};
            }
            // Draw the line
            std::string uid = std::to_string(iet) + "_" + std::to_string(ie);
            auto e_line =
                draw::line(id_ + "eta_line_" + uid, start, end, stroke);
            e.add_object(e_line);
            // Label it if told to do so
            if (std::get<bool>(elt)) {
                auto font = std::get<style::font>(elt);
                end[0] +=
                    static_cast<scalar>(std::cos(theta) * 0.5 * font._size);
                end[1] +=
                    static_cast<scalar>(std::sin(theta) * 0.5 * font._size);
                if (eta == 0.) {
                    end[0] -= static_cast<scalar>(0.5 * font._size);
                }
                const std::string eta_code = "\u03B7";
                auto e_text = eta_code + " = " + utils::to_string(eta);
                auto e_label =
                    draw::text(id_ + "eta_label_" + uid, end, {e_text}, font);
                e.add_object(e_label);
            }
        }
    }
    return e;
}

}  // namespace display

}  // namespace actsvg
