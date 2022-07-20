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

TEST(core, copy_triangle) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_copy_triangle.svg");

    std::vector<std::array<scalar, 2u>> triangle = {
        {-100, -100}, {-100, 100}, {100, 100}};

    auto original = draw::polygon("original", triangle);
    original._sterile = true;

    // A red filled copy
    style::fill red_fill(style::color({{255, 0, 0}}));
    auto red_copy = draw::from_template("red_copy", original, red_fill);

    // A blue filled copy, shifted
    style::transform shift_x100;
    shift_x100._tr = {100, 0};
    style::fill blue_fill(style::color({{0, 0, 255}, 0.5}));
    auto blue_copy = draw::from_template("blue_copy", original, blue_fill,
                                         style::stroke(), shift_x100);

    // A green filled copy, shifted and rotated
    style::transform shift_x200_r45;
    shift_x200_r45._tr = {200, 0};
    shift_x200_r45._rot[0] = -45.;
    style::fill green_fill(style::color({{0, 255, 0}, 0.5}));
    auto green_copy = draw::from_template("green_copy", original, green_fill,
                                          style::stroke(), shift_x200_r45);

    // Label
    auto descr = draw::text(
        "descr", {-350, 300},
        {"[-] original (sterile) triganle: not shown", "[x] red copy at (0,0)",
         "[x] blue copy shifted by (100, 0)",
         "[x] green copy shifted by (200,0), rotated for 45 degrees"});

    svg::file of;
    of.add_object(pg);
    of.add_object(red_copy);
    of.add_object(blue_copy);
    of.add_object(green_copy);
    of.add_object(descr);
    fo << of;
    fo.close();
}
