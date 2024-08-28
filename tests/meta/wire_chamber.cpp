// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/proto/cluster.hpp"
#include "actsvg/proto/surface.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

/// Helper method to create a wire chamber
///
/// @param chamber_name the name of the chamber
/// @param chamber_offset_y where the straws are posistioned
/// @param r_tube the tube radius
/// @param n_tubes the number of tubes
/// @param n_layers the number of wire layers
///
/// @return a vector of surfaces
std::vector<proto::surface<point3_container>> wire_chamber(
    const std::string chamber_name, scalar chamber_offset_y = 0.,
    scalar r_tube = 20., unsigned int n_tubes = 40u,
    unsigned int n_layers = 4u) {

    // Create the surfaces
    std::vector<proto::surface<point3_container>> wires;

    // packing in r and offset in y & length in x
    scalar r_pack = r_tube * static_cast<scalar>(std::sqrt(3.));
    scalar c_y = (n_layers % 2) ? -(n_layers - 1) / 2 * r_pack
                                : -(n_layers / 2 - 0.5_scalar) * r_pack;
    scalar m_lx = (2 * n_tubes + 1) * r_tube;
    for (unsigned int il = 0u; il < n_layers; ++il) {
        // positioningin y
        scalar t_y = c_y + il * r_pack;
        scalar offset_x = (il % 2) * r_tube;
        for (unsigned it = 0u; it < n_tubes; ++it) {
            // The x position
            scalar t_x = -0.5_scalar * m_lx + (2 * it + 1) * r_tube + offset_x;
            // Creat a wire
            proto::surface<point3_container> wire;
            wire._radii = {1., r_tube};
            wire._zparameters = {0., 200};
            wire._name = chamber_name + "_wire_" + std::to_string(il);
            wire._type = proto::surface<point3_container>::type::e_straw;
            wire._fill = style::fill({{200, 200, 200}, 0.5});
            wire._stroke = style::stroke({{0, 0, 0}}, 1.);
            wire._stroke._hl_width = 2.;
            wire._transform._tr = {t_x, t_y + chamber_offset_y};
            wires.push_back(wire);
        }
    }
    return wires;
}

/// Helper method to find the closest approach of the line
///
///
std::optional<proto::cluster<1u>> drift_cluster(
    const proto::surface<point3_container>& s_, const point2& start_,
    scalar alpha_) {

    auto t = s_._transform._tr;

    scalar d_r = std::abs(std::cos(alpha_) * (start_[1u] - t[1u]) -
                          std::sin(alpha_) * (start_[0u] - t[0u]));

    if (d_r <= s_._radii[1u]) {
        proto::cluster<1u> drift_cluster{};
        drift_cluster._type = proto::cluster<1u>::type::e_drift;
        drift_cluster._coords = {proto::cluster<1u>::coordinate::e_r};
        drift_cluster._measurement = {d_r};
        return drift_cluster;
    }
    return std::nullopt;
}

TEST(proto, wire_chamber) {

    // Files: geometry (+track)
    svg::file gfile;

    // The surfaces
    std::vector<proto::surface<point3_container>> all_wires;
    std::vector<point2> all_wires_positions;
    std::vector<scalar> s_offsets = {-200, 200};
    for (auto [io, s_o] : utils::enumerate(s_offsets)) {
        auto multi_wires =
            wire_chamber("multilayer_" + std::to_string(io), s_o);
        all_wires.insert(all_wires.end(), multi_wires.begin(),
                         multi_wires.end());
    }

    for (const auto& s : all_wires) {
        svg::object wire = display::surface(s._name, s, views::x_y{});
        // Collect the wire positions
        auto tr = s._transform._tr;
        all_wires_positions.push_back({tr[0], tr[1]});
        gfile.add_object(wire);
    }

    std::ofstream rstream;
    rstream.open("test_meta_wire_chamber.svg");
    rstream << gfile;
    rstream.close();
}

TEST(proto, wire_chamber_with_track) {

    // Files: geometry with track
    svg::file tfile;

    std::vector<proto::surface<point3_container>> all_wires;
    std::vector<point2> all_wires_positions;
    std::vector<scalar> s_offsets = {-200, 200};
    for (auto [io, s_o] : utils::enumerate(s_offsets)) {
        auto multi_wires =
            wire_chamber("multilayer_" + std::to_string(io), s_o, 20u, 20u);
        all_wires.insert(all_wires.end(), multi_wires.begin(),
                         multi_wires.end());
    }

    // Track start
    point2 start = {-200, -300};
    scalar alpha = 1.3_scalar;

    for (const auto& s : all_wires) {
        // Collect the wire positions
        const auto& tr = s._transform._tr;
        all_wires_positions.push_back({tr[0u], tr[1u]});
        // Intersect the wire
        auto dc = drift_cluster(s, start, alpha);
        if (dc.has_value()) {
            // Make a copy of the surface
            proto::surface<point3_container> sc = s;
            sc._fill = style::fill({{0, 100, 200}, 0.5});
            // Draw the copy
            tfile.add_object(display::surface(s._name, sc, views::x_y{}));
            // Add the drift circle
            const auto& dcv = dc.value();
            auto dco = draw::circle("dc_wire_" + s._name, {tr[0], tr[1]},
                                    dcv._measurement[0u]);
            tfile.add_object(dco);
        } else {
            tfile.add_object(display::surface(s._name, s, views::x_y{}));
        }
    }

    // Create the track
    style::stroke track_stroke = style::stroke({});
    track_stroke._width = 2;

    style::marker track_head = style::marker({"<<"});

    scalar length = 700.;
    point2 direction = {static_cast<scalar>(std::cos(alpha)),
                        static_cast<scalar>(std::sin(alpha))};
    point2 end = {start[0] + direction[0] * length,
                  start[1] + direction[1] * length};
    auto track = draw::measure("track", start, end, track_stroke,
                               style::marker({""}), track_head);
    tfile.add_object(track);

    std::ofstream rstream;
    rstream.open("test_meta_wire_chamber_with_track.svg");
    rstream << tfile;
    rstream.close();
}
