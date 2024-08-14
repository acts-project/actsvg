// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <array>
#include <fstream>

#include "../common/helix.hpp"
#include "actsvg/display/datamodel.hpp"
#include "actsvg/proto/track.hpp"
#include "actsvg/styles/defaults.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3u>;

TEST(proto, seeds_single_xy) {

    views::x_y xy_view;

    std::vector<scalar> radii = {30., 50., 70.};
    std::vector<svg::object> cylinders = {};
    for (auto [ir, r] : utils::enumerate(radii)) {
        cylinders.push_back(
            draw::circle("l" + std::to_string(ir), {0., 0.}, r, __nn_fill));
    }

    point3 origin = {0., 0., 0.};
    scalar oot = 1._scalar / std::sqrt(3._scalar);
    point3 direction = {oot, oot, oot};
    scalar radius = 100.;
    scalar o = 1.;
    std::array<scalar, 2u> phi_range = {-0.1_scalar, 0.35_scalar * pi};
    proto::trajectory<point3> trj = test::generate_helix<point3>(
        origin, direction, radius, o, phi_range, 100);
    trj._origin_size = 0.;
    trj._origin_fill = defaults::__g_fill;
    trj._origin_stroke = style::stroke{style::color{{0, 0, 0}}, 0.5};
    trj._path_stroke = style::stroke{style::color{{0, 255, 0}}, 0.5};
    trj._path_arrow = style::marker{"<<", 3., style::color{{0, 255, 0}}};

    proto::seed<point3> seed;
    scalar s0x = 30._scalar * std::cos(0.3_scalar * pi);
    scalar s0y = 30._scalar * std::sin(0.3_scalar * pi);
    scalar s1x = 50._scalar * std::cos(0.33_scalar * pi);
    scalar s1y = 50._scalar * std::sin(0.33_scalar * pi);
    scalar s2x = 70._scalar * std::cos(0.365_scalar * pi);
    scalar s2y = 70._scalar * std::sin(0.365_scalar * pi);

    seed._space_points = {point3{s0x, s0y, 0.}, point3{s1x, s1y},
                          point3{s2x, s2y, 0.}};
    seed._trajectory = trj;
    seed._space_point_size = 1.;
    seed._space_point_fill = style::fill(style::color{{0, 255, 0}});

    auto seed_traj_xy = display::seed("seed_xy", seed, xy_view);

    std::ofstream rstream;

    svg::file seed_xy;
    seed_xy.add_objects(cylinders);
    seed_xy.add_object(seed_traj_xy);

    rstream.open("test_meta_single_seed_xy.svg");
    rstream << seed_xy;
    rstream.close();
}
