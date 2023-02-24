// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <array>
#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/proto/surface.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

std::vector<proto::surface<point3_container>> wire_chamber(
    const std::string chamber_name, scalar chamber_offset_y = 0.,
    scalar r_tube = 20., unsigned int n_tubes = 40u,
    unsigned int n_layers = 4u) {

    // Create the surfaces
    std::vector<proto::surface<point3_container>> wires;

    // packing in r and offset in y & length in x
    scalar r_pack = r_tube * std::sqrt(3.);
    scalar c_y = (n_layers % 2) ? -(n_layers - 1) / 2 * r_pack
                                : -(n_layers / 2 - 0.5) * r_pack;
    scalar m_lx = (2 * n_tubes + 1) * r_tube;
    for (unsigned int il = 0u; il < n_layers; ++il) {
        // positioningin y
        scalar t_y = c_y + il * r_pack;
        scalar offset_x = (il % 2) * r_tube;
        for (unsigned it = 0u; it < n_tubes; ++it) {
            // The x position
            scalar t_x = -0.5 * m_lx + (2 * it + 1) * r_tube + offset_x;
            // Creat a wire
            proto::surface<point3_container> wire;
            wire._radii = {1., r_tube};
            wire._zparameters = {0., 200};
            wire._name = chamber_name + "_wire_" + std::to_string(il);
            wire._type = proto::surface<point3_container>::type::e_straw;
            wire._fill = style::fill({{0, 100, 0}, 0.5});
            wire._stroke = style::stroke({{0, 0, 0}}, 1.);
            wire._stroke._hl_width = 2.;
            wire._transform._tr = {t_x, t_y + chamber_offset_y};
            wires.push_back(wire);
        }
    }
    return wires;
}

TEST(proto, wire_chamber) {

    // Create the multiwire chamber
    svg::file rfile;

    std::vector<proto::surface<point3_container>> all_wires;
    std::vector<scalar> s_offsets = {-200, 200};
    for (auto [io, s_o] : utils::enumerate(s_offsets)) {
        auto multi_wires =
            wire_chamber("multilayer_" + std::to_string(io), s_o);
        all_wires.insert(all_wires.end(), multi_wires.begin(),
                         multi_wires.end());
    }

    for (const auto& s : all_wires) {
        svg::object wire = display::surface(s._name, s, views::x_y{});
        rfile.add_object(wire);
    }

    std::ofstream rstream;
    rstream.open("test_meta_wire_chamber.svg");
    rstream << rfile;
    rstream.close();
}
