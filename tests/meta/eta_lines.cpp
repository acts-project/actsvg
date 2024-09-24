// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "../common/playground.hpp"
#include "actsvg/core.hpp"
#include "actsvg/display/geometry.hpp"

using namespace actsvg;

TEST(meta, eta_lines_positive_side) {

    auto pg = test::playground({-300, -10}, {300, 120});

    svg::file etafile;
    etafile.add_object(pg);

    scalar r = 120;
    scalar z = 300;

    std::vector<scalar> eta_main = {0, 1, 2, 3};
    style::stroke stroke_main;
    bool label_main = true;
    style::font font_main;

    std::vector<scalar> eta_half = {0.5, 1.5, 2.5};
    style::stroke stroke_half;
    stroke_half._dasharray = {2, 2};
    bool label_half = false;
    style::font font_half;

    auto eta_lines = display::eta_lines(
        "eta_lines_", z, r,
        {std::tie(eta_main, stroke_main, label_main, font_main),
         std::tie(eta_half, stroke_half, label_half, font_half)});

    etafile.add_object(eta_lines);
    std::ofstream tstream;
    tstream.open("test_display_eta_lines_ps.svg");
    tstream << etafile;
    tstream.close();
}

TEST(meta, eta_lines_full_detector) {

    auto pg = test::playground({-300, -10}, {300, 120});

    svg::file etafile;
    etafile.add_object(pg);

    scalar r = 120;
    scalar z = 300;

    std::vector<scalar> eta_main = {-3, -2, -1, 0, 1, 2, 3};
    style::stroke stroke_main;
    bool label_main = true;
    style::font font_main;

    std::vector<scalar> eta_half = {0.5, 1.5, 2.5};
    style::stroke stroke_half;
    stroke_half._dasharray = {2, 2};
    bool label_half = false;
    style::font font_half;

    auto eta_lines = display::eta_lines(
        "eta_lines_", z, r,
        {std::tie(eta_main, stroke_main, label_main, font_main),
         std::tie(eta_half, stroke_half, label_half, font_half)});

    etafile.add_object(eta_lines);
    std::ofstream tstream;
    tstream.open("test_display_eta_lines_fd.svg");
    tstream << etafile;
    tstream.close();
}
