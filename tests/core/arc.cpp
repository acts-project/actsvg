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

TEST(core, arc_plain) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_arc.svg");

    fo << ftemplate._html_head;
    fo << ftemplate._svg_head;
    fo << " width=\"900\" height=\"900\" viewBox=\"-450 -450 900 900\"";
    fo << ftemplate._svg_def_end;
    // Add the playground
    fo << pg;

    scalar phi_min = -0.25;
    scalar phi_max = 0.75;
    scalar r = 125.;

    point2 start = { r * std::cos(phi_min), r * std::sin(phi_min)};
    point2 end = { r * std::cos(phi_max), r * std::sin(phi_max)};

    // Add the line
    fo << draw::arc("a", r, start, end, style::fill(),
                     style::stroke{{{255, 0, 0}}, 2});
    // Close the file
    fo << ftemplate._svg_tail;
    fo << ftemplate._html_tail;
    fo.close();
}
