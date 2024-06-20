// This file is part of the actsvg packge.
//
// Copyright (C) 2024 CERN for the benefit of the ACTS project
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

TEST(core, gradient_box_linear_x) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    style::gradient gradient;
    gradient._id = "g_def";
    gradient._stops = {
        style::gradient::stop{0., style::color{style::rgb{0, 0, 255}}},
        style::gradient::stop{0.25, style::color{style::rgb{0, 255, 0}}},
        style::gradient::stop{0.75, style::color{style::rgb{255, 255, 0}}},
        style::gradient::stop{1., style::color{style::rgb{255, 0, 0}}}};

    auto g = draw::gradient_as_object(gradient);

    // Create the box
    svg::object gbox;
    gbox._tag = "rect";
    gbox._id = "gbox";
    gbox._sterile = true;
    gbox._attribute_map["x"] = std::to_string(0.);
    gbox._attribute_map["y"] = std::to_string(0.);
    gbox._attribute_map["width"] = std::to_string(200);
    gbox._attribute_map["height"] = std::to_string(200);
    gbox._attribute_map["fill"] = "url(#" + gradient._id + ")";

    svg::file mfile;
    mfile.add_object(pg);
    mfile.add_object(g);
    mfile.add_object(gbox);

    std::ofstream tstream;
    tstream.open("test_core_gradient_box.svg");
    tstream << mfile;
    tstream.close();
}

TEST(core, gradient_box_vertical) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    auto g = draw::gradient_box(
        "g_box", {100, 20}, {20, 200},
        {{style::gradient::stop{0., style::color{style::rgb{0, 0, 255}}}, 0.},
         {style::gradient::stop{0.45, style::color{style::rgb{0, 255, 0}}},
          2.0},
         {style::gradient::stop{0.75, style::color{style::rgb{255, 255, 0}}},
          2.6},
         {style::gradient::stop{1., style::color{style::rgb{255, 0, 0}}},
          5.2}});

    svg::file mfile;
    mfile.add_object(pg);
    mfile.add_object(g);

    std::ofstream tstream;
    tstream.open("test_core_gradient_box_vertical.svg");
    tstream << mfile;
    tstream.close();
}

TEST(core, gradient_box_vertical_label) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    auto g = draw::gradient_box(
        "g_box", {100, 20}, {20, 200},
        {{style::gradient::stop{0., style::color{style::rgb{0, 0, 255}}}, 0.},
         {style::gradient::stop{0.45, style::color{style::rgb{0, 255, 0}}},
          2.0},
         {style::gradient::stop{0.75, style::color{style::rgb{255, 255, 0}}},
          2.6},
         {style::gradient::stop{1., style::color{style::rgb{255, 0, 0}}},
          5.2}});

    svg::file mfile;
    mfile.add_object(pg);
    mfile.add_object(g);

    std::ofstream tstream;
    tstream.open("test_core_gradient_box_vertical_label.svg");
    tstream << mfile;
    tstream.close();
}

TEST(core, gradient_box_horizontal) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    auto g = draw::gradient_box(
        "g_box", {50, 50}, {200, 20},
        {{style::gradient::stop{0., style::color{style::rgb{0, 0, 255}}}, 0.},
         {style::gradient::stop{0.25, style::color{style::rgb{0, 255, 0}}},
          1.3},
         {style::gradient::stop{0.75, style::color{style::rgb{255, 255, 0}}},
          2.6},
         {style::gradient::stop{1., style::color{style::rgb{255, 0, 0}}},
          5.2}});

    svg::file mfile;
    mfile.add_object(pg);
    mfile.add_object(g);

    std::ofstream tstream;
    tstream.open("test_core_gradient_box_horizontal.svg");
    tstream << mfile;
    tstream.close();
}

TEST(core, gradient_box_horizontal_label) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    auto g = draw::gradient_box(
        "g_box", {50, 50}, {200, 20},
        {{style::gradient::stop{0., style::color{style::rgb{0, 0, 255}}}, 0.},
         {style::gradient::stop{0.25, style::color{style::rgb{0, 255, 0}}},
          1.3},
         {style::gradient::stop{0.75, style::color{style::rgb{255, 255, 0}}},
          2.6},
         {style::gradient::stop{1., style::color{style::rgb{255, 0, 0}}}, 5.2}},
        style::label{"arb. unit", style::label::horizontal::right,
                     style::label::vertical::top});

    svg::file mfile;
    mfile.add_object(pg);
    mfile.add_object(g);

    std::ofstream tstream;
    tstream.open("test_core_gradient_box_horizontal_label.svg");
    tstream << mfile;
    tstream.close();
}