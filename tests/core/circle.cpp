// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>

#include "../common/playground.hpp"
#include "actsvg/core/draw.hpp"
#include "actsvg/core/style.hpp"

using namespace actsvg;

TEST(core, circle_plain) {
    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_circle.svg");

    // Draw the circles
    auto cc = draw::circle("cc", {0., 0.}, 10,
                           style::fill(style::color{{255, 255, 255}}));
    auto c = draw::circle("c", {40., -20.}, 25.,
                          style::fill{style::color{{0, 125, 0}}});
    // Label
    auto descr = draw::text("descr", {100, 100},
                            {"one circle in white with r=10 at (0,0)",
                             "one circle in green with r = 25 at (40,-20)"});

    svg::file of;
    of.add_object(pg);
    of.add_object(cc);
    of.add_object(c);
    of.add_object(descr);
    fo << of;
    fo.close();
}

TEST(core, circle_shifted) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    style::transform t{{100, 100}};

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_circle_shifted.svg");

    auto c = draw::circle("c", {40., -20.}, 25.,
                          style::fill{style::color{{0, 125, 0}}},
                          style::stroke{style::color{{0, 0, 0}}}, t);
    // Label
    auto descr = draw::text("descr", {20, 140},
                            {"shifted circle to (140, 80) with r=25"});

    svg::file of;
    of.add_object(pg);
    of.add_object(c);
    of.add_object(descr);
    fo << of;
    fo.close();
}

TEST(core, circle_scaled) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    style::transform t{{0, 0}};
    t._scale = {25, 25};

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_circle_scaled.svg");

    auto c = draw::circle("c", {-10., -10.}, 4,
                          style::fill{style::color{{0, 125, 0}}},
                          style::stroke{style::color{{0, 0, 0}}}, t);

    // Label
    auto descr = draw::text("descr", {20, 140},
                            {"scaled circle at (-250, -250) with r=50"});

    svg::file of;
    of.add_object(pg);
    of.add_object(c);
    of.add_object(descr);
    fo << of;
    fo.close();
}
