// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>

#include "actsvg/core.hpp"

using namespace actsvg;

TEST(draw, triangle) {

    svg::file tfile0;
    std::vector<std::array<scalar, 2u>> triangle = {
        {100, 400}, {200, 600}, {600, 200}};

    auto tsvg0 = draw::polygon("t0", triangle);
    tfile0._objects.push_back(tsvg0);

    std::ofstream tstream;
    tstream.open("triangle.svg");
    tstream << tfile0;
    tstream.close();

    style::fill red_fill_hl_green({255, 0, 0});
    red_fill_hl_green._fc._highlight = {"mouseover", "mouseout"};
    red_fill_hl_green._fc._hl_rgb = {0, 255, 0};

    svg::file tfile1;
    auto tsvg1 = draw::polygon("t1", triangle, red_fill_hl_green);
    tfile1._objects.push_back(tsvg1);

    tstream.open("triangle_highlight.svg");
    tstream << tfile1;
    tstream.close();
}
