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

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_circle.svg");

    fo << ftemplate._html_head;
    fo << ftemplate._svg_head;
    fo << " width=\"900\" height=\"900\" viewBox=\"-450 -450 900 900\"";
    fo << ftemplate._svg_def_end;
    // Add the playground
    fo << pg;
    // Draw the circles
    fo << draw::circle("cc", {0., 0.}, 10,
                       style::fill(style::color{{255, 255, 255}}));
    fo << draw::circle("c", {40., -20.}, 25.,
                       style::fill{style::color{{0, 125, 0}}});
    // Close the file
    fo << ftemplate._svg_tail;
    fo << ftemplate._html_tail;
    fo.close();
}

TEST(core, circle_shifted) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    style::transform t{{100, 100}};

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_circle_shifted.svg");

    fo << ftemplate._html_head;
    fo << ftemplate._svg_head;
    fo << " width=\"900\" height=\"900\" viewBox=\"-450 -450 900 900\"";
    fo << ftemplate._svg_def_end;
    // Add the playground
    fo << pg;
    // Add the circle
    fo << draw::circle("c", {40., -20.}, 25.,
                       style::fill{style::color{{0, 125, 0}}},
                       style::stroke{style::color{{0, 0, 0}}}, t);
    // Close the file
    fo << ftemplate._svg_tail;
    fo << ftemplate._html_tail;
    fo.close();
}

TEST(core, circle_scaled) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    style::transform t{{0, 0}};
    t._scale = {10, 10};

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_circle_scaled.svg");

    fo << ftemplate._html_head;
    fo << ftemplate._svg_head;
    fo << " width=\"900\" height=\"900\" viewBox=\"-450 -450 900 900\"";
    fo << ftemplate._svg_def_end;
    // Add the playground
    fo << pg;
    // Add the circle
    fo << draw::circle("c", {4., -2.}, 2.5,
                       style::fill{style::color{{0, 125, 0}}},
                       style::stroke{style::color{{0, 0, 0}}}, t);
    // Close the file
    fo << ftemplate._svg_tail;
    fo << ftemplate._html_tail;
    fo.close();
}
