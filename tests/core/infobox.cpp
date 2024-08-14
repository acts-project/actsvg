// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

#include "../common/playground.hpp"
#include "actsvg/core.hpp"

using namespace actsvg;

TEST(core, info_box) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    std::string title = "info box";
    style::fill title_blue;
    title_blue._fc._rgb = {0, 0, 150};
    title_blue._fc._opacity = 0.8_scalar;
    style::font title_font;
    title_font._size = 24;
    title_font._fc = style::color{{255, 255, 255}};

    std::vector<std::string> info = {"this is an info box", "a = 101",
                                     "this line should define the length"};

    style::fill info_blue;
    info_blue._fc._rgb = {0, 0, 150};
    info_blue._fc._opacity = 0.4_scalar;
    style::font info_font;
    info_font._size = 18;

    style::stroke stroke;

    auto p = draw::circle("test_circle", {-100, 100}, 15., title_blue, stroke);

    auto t =
        draw::text("test_instruction", {-100, 120},
                   {"Move the mouse cursor over the the circle"}, info_font);

    auto info_box = draw::connected_info_box("test_box_", {100, 100}, title,
                                             title_blue, title_font, info,
                                             info_blue, info_font, stroke, p);

    svg::file ifile;
    ifile.add_object(pg);
    ifile.add_object(p);
    ifile.add_object(t);
    ifile.add_object(info_box);

    std::ofstream tstream;
    tstream.open("test_core_infobox.svg");
    tstream << ifile;
    tstream.close();
}
