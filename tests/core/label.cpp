// This file is part of the actsvg packge.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
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

namespace {
// Helper function
void test_label(style::label::horizontal h, style::label::vertical v,
                const std::string& name) {
    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Write out the file
    std::ofstream fo;
    fo.open(std::string("test_core_label_" + name + ".svg").c_str());

    auto box = draw::rectangle("box", {200, 100}, 150, 20,
                               style::fill{style::color{{255, 0, 0}}});
    auto [llc, urc] = box.generate_bounding_box();
    auto l = style::label{name, h, v};
    l.place(llc, urc);

    svg::file of;
    of.add_object(pg);
    of.add_object(box);
    of.add_object(draw::label("label", l));
    fo << of;
    fo.close();
}
}  // namespace

TEST(core, label_left_bottom) {
    test_label(style::label::horizontal::left, style::label::vertical::bottom,
               "left_bottom");
}

TEST(core, label_right_bottom) {
    test_label(style::label::horizontal::right, style::label::vertical::bottom,
               "right_bottom");
}

TEST(core, label_left_top) {
    test_label(style::label::horizontal::left, style::label::vertical::top,
               "left_top");
}

TEST(core, label_right_top) {
    test_label(style::label::horizontal::right, style::label::vertical::top,
               "right_top");
}

TEST(core, label_center_center) {
    test_label(style::label::horizontal::center, style::label::vertical::center,
               "center_center");
}

TEST(core, label_left_center) {
    test_label(style::label::horizontal::left, style::label::vertical::center,
               "left_center");
}

TEST(core, label_right_center) {
    test_label(style::label::horizontal::right, style::label::vertical::center,
               "right_center");
}