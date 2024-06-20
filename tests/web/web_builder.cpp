// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/web/web_builder.hpp"

#include <gtest/gtest.h>

#include <filesystem>

#include "actsvg/core/draw.hpp"

using namespace actsvg;

TEST(web, web_builder) {
    // Draw the circles
    auto cc = draw::circle("circle 1", {0., 0.}, 10,
                           style::fill(style::color{{255, 255, 255}}));
    auto c = draw::circle("circle 2", {40., -20.}, 25.,
                          style::fill{style::color{{0, 125, 0}}});

    web::web_builder wb;
    wb.build(std::filesystem::current_path() / "test_web_web_builder",
             std::vector{c, cc});
}
