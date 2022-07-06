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

using namespace actsvg;

TEST(core, line) {

    auto line = draw::line("l", {40., -20.}, {80., 100.},
                           style::stroke{style::color{{255, 0, 0}}, 2});

    std::array<scalar, 2> ref_x_range = {40., 80.};
    std::array<scalar, 2> ref_y_range = {-100., 20.};
    ASSERT_TRUE(line._x_range == ref_x_range);
    ASSERT_TRUE(line._y_range == ref_y_range);
}

TEST(core, line_plain) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_line.svg");

    fo << ftemplate._html_head;
    fo << ftemplate._svg_head;
    fo << " width=\"900\" height=\"900\" viewBox=\"-450 -450 900 900\"";
    fo << ftemplate._svg_def_end;
    // Add the playground
    fo << pg;
    // Add the line
    fo << draw::line("l", {40, -20.}, {80., 100.},
                     style::stroke{style::color{{255, 0, 0}}, 2});
    // Close the file
    fo << ftemplate._svg_tail;
    fo << ftemplate._html_tail;
    fo.close();
}

TEST(core, line_shifted) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_line_shifted.svg");

    style::transform t{{100, 100}};

    fo << ftemplate._html_head;
    fo << ftemplate._svg_head;
    fo << " width=\"900\" height=\"900\" viewBox=\"-450 -450 900 900\"";
    fo << ftemplate._svg_def_end;
    // Add the playground
    fo << pg;
    // Add the line
    fo << draw::line("l", {40, -20.}, {80., 100.},
                     style::stroke{style::color{{0, 255, 0}}, 2}, t);
    // Close the file
    fo << ftemplate._svg_tail;
    fo << ftemplate._html_tail;
    fo.close();
}

TEST(core, line_scaled) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_line_scaled.svg");

    style::transform t;
    t._scale = {10, 10};

    fo << ftemplate._html_head;
    fo << ftemplate._svg_head;
    fo << " width=\"900\" height=\"900\" viewBox=\"-450 -450 900 900\"";
    fo << ftemplate._svg_def_end;
    // Add the playground
    fo << pg;
    // Add the line
    fo << draw::line("l", {4, -2.}, {8., 10.},
                     style::stroke{style::color{{0, 0, 255}}, 2}, t);
    // Close the file
    fo << ftemplate._svg_tail;
    fo << ftemplate._html_tail;
    fo.close();
}
