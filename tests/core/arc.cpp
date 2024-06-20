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
    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_arc.svg");

    scalar phi_min = -0.25;
    scalar phi_max = 0.75;
    scalar r = 125.;

    point2 start = {r * std::cos(phi_min), r * std::sin(phi_min)};
    point2 end = {r * std::cos(phi_max), r * std::sin(phi_max)};

    // Add the line
    auto a = draw::arc("a", r, start, end, style::fill(),
                       style::stroke{{{255, 0, 0}}, 2});

    auto descr = draw::text("descr", {150, 100},
                            {"arc at r = 125", "with phi in [-0.25,0.75]"});

    svg::file of;
    of.add_object(pg);
    of.add_object(a);
    of.add_object(descr);
    fo << of;
    fo.close();
}
