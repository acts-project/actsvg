// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include "actsvg/core/draw.hpp"
#include "actsvg/core/style.hpp"
#include "../common/playground.hpp"

using namespace actsvg;

TEST(text, unconnected_text) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_text.svg");

    fo << ftemplate._html_head;
    fo << ftemplate._svg_head;
    fo << " width=\"900\" height=\"900\" viewBox=\"-450 -450 900 900\"";
    fo << ftemplate._svg_def_end;
    // Add the playground
    fo << pg;
 
    style::color red{{255, 0, 0}};
    style::font fs;
    fs._family = "Arial";
    fs._fc = red;

    // Add the text
    fo << draw::text("t0", {10, 10}, {"Arial test text at (10,10)"}, fs);
    // Close the file
    fo << ftemplate._svg_tail;
    fo << ftemplate._html_tail;
    fo.close();

}


TEST(text, multiline_text) {

    svg::file ftemplate;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    style::color red{{255, 0, 0}};
    style::font fs;
    fs._family = "Arial";
    fs._fc = red;

    // Write out the file
    std::ofstream fo;
    fo.open("test_core_text_multiline.svg");

    fo << ftemplate._html_head;
    fo << ftemplate._svg_head;
    fo << " width=\"900\" height=\"900\" viewBox=\"-450 -450 900 900\"";
    fo << ftemplate._svg_def_end;
    // Add the playground
    fo << pg;

    // Add the text
    fo << draw::text("t0", {100, 100}, {"line 0", "line 1"}, fs);
    // Close the file
    fo << ftemplate._svg_tail;
    fo << ftemplate._html_tail;
    fo.close();
}