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
#include "../common/playground.hpp"
#include "actsvg/display/datamodel.hpp"
#include "actsvg/proto/track.hpp"
#include "actsvg/styles/defaults.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3u>;

TEST(proto, helix_trajectory) {

    views::x_y xy_view;

    point3 origin = {0., 0., 0.};
    scalar oot = 1._scalar / std::sqrt(3._scalar);
    point3 direction = {oot, oot, oot};
    scalar radius = 100.;
    scalar o = 1.;
    std::array<scalar, 2u> phi_range = {-0._scalar, 0.75_scalar * pi};
    proto::trajectory<point3> trj_p = test::generate_helix<point3>(
        origin, direction, radius, o, phi_range, 100);
    trj_p._origin_size = 5.;
    trj_p._origin_fill = defaults::__g_fill;
    trj_p._origin_stroke = style::stroke{style::color{{0, 0, 0}}, 1.5};
    trj_p._path_stroke = style::stroke{style::color{{0, 255, 0}}, 1.5};

    auto trj_xy_polyline_p =
        display::trajectory("helix_xy_polyline_p", trj_p, xy_view);

    proto::trajectory<point3> trj_n = test::generate_helix<point3>(
        origin, direction, 1._scalar * radius, -o, phi_range, 100);
    trj_n._origin_size = 5.;
    trj_n._origin_fill = defaults::__g_fill;
    trj_n._origin_stroke = style::stroke{style::color{{0, 0, 0}}, 1.5};
    trj_n._path_arrow =
        style::marker{">", 5.5, style::fill{style::color{{255, 0, 0}}}};
    trj_n._path_stroke = style::stroke{style::color{{255, 0, 0}}, 1.5};

    auto trj_xy_polyline_n =
        display::trajectory("helix_xy_polyline_n", trj_n, xy_view);

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file hfile;
    hfile.add_object(pg);
    hfile.add_object(trj_xy_polyline_p);
    hfile.add_object(trj_xy_polyline_n);

    std::ofstream rstream;
    rstream.open("test_meta_helix_polyline_xy.svg");
    rstream << hfile;
    rstream.close();
}
