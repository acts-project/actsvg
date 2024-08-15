// This file is part of the actsvg packge.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>

#include "../common/playground.hpp"
#include "actsvg/core/draw.hpp"

using namespace actsvg;

TEST(core, bezier_single_segment) {

    scalar ref = 75.;

    std::vector<std::array<point2, 2u>> points_dirs;
    points_dirs.push_back({point2{ref, 0}, point2{0, 1.}});
    points_dirs.push_back({point2{static_cast<scalar>(ref * std::sqrt(1. / 2)),
                                  static_cast<scalar>(ref * std::sqrt(1. / 2))},
                           point2{-static_cast<scalar>(std::sqrt(1. / 2)),
                                  static_cast<scalar>(std::sqrt(1. / 2))}});

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_bezier_single.svg");

    auto c = draw::bezier("bezier", points_dirs);

    auto descr = draw::text("descr", {150, 100}, {"Bezier segment"});

    svg::file of;
    of.add_object(pg);
    of.add_object(c);
    of.add_object(descr);
    fo << of;
    fo.close();
}

TEST(core, bezier_spiral) {

    unsigned int nSegs = 200.;
    scalar r_start = 125.;
    scalar r_end = 75.;
    scalar phi_start = -0.25;
    scalar phi_rot = 4.5_scalar * pi;

    scalar r_step = (r_end - r_start) / nSegs;
    scalar phi_step = phi_rot / nSegs;

    std::vector<std::array<point2, 2u>> points_dirs;
    for (unsigned int is = 0; is < nSegs; ++is) {
        scalar r = r_start + is * r_step;
        scalar phi = phi_start + is * phi_step;
        point2 p = {r * std::cos(phi), r * std::sin(phi)};
        point2 d = {-std::sin(phi), std::cos(phi)};
        points_dirs.push_back({p, d});
    }

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_bezier_spiral.svg");

    auto c = draw::bezier("bezier", points_dirs);
    auto descr =
        draw::text("descr", {150, 100}, {"Sprial drawn with Bezier segments"});

    svg::file of;
    of.add_object(pg);
    of.add_object(c);
    of.add_object(descr);
    fo << of;
    fo.close();
}
