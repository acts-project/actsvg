// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>

#include "../common/playground.hpp"
#include "actsvg/core.hpp"

using namespace actsvg;

TEST(draw, rectangle) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file tfile0;

    style::fill red(style::color{{255, 0, 0}});
    red._fc._highlight = {"mouseover", "mouseout"};
    red._fc._hl_rgb = {{0, 255, 0}};

    style::fill green(style::color{{0, 255, 0}});
    green._fc._highlight = {"mouseover", "mouseout"};
    green._fc._hl_rgb = {{0, 0, 255}};

    auto r0 = draw::rectangle("r0", {0, 0}, 20, 10, red, style::stroke());
    auto r1 = draw::rectangle("r1", {100, 100}, 15, 40, green, style::stroke());

    tfile0._objects.push_back(pg);
    tfile0._objects.push_back(r0);
    tfile0._objects.push_back(r1);

    std::ofstream tstream;
    tstream.open("test_core_rectangle.svg");
    tstream << tfile0;
    tstream.close();
}