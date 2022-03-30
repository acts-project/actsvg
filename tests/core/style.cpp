// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fstream>
#include <gtest/gtest.h>

#include "actsvg/core/svg.hpp"
#include "actsvg/core/style.hpp"

#include <iostream>
#include <sstream>

using namespace actsvg;

TEST(svg, fill_style)
{

    svg::object colored{"red_object"};

    style::color red{{255, 0, 0}};
    style::fill red_fill{red};
    red_fill.attach_attributes(colored);

    std::cout << colored << std::endl;

    svg::object highlight{"red_highlight_object"};

    style::color red_hl{{255, 0, 0}};
    red_hl._hl_rgb = {0, 255, 0};
    red_hl._highlight = {"mouseover", "mouseout"};
    style::fill red_hl_fill{red_hl};
    red_hl_fill.attach_attributes(highlight);

    std::cout << highlight << std::endl;
}

TEST(svg, stroke_style)
{

    svg::object stroked{"stroked_object"};

    style::color black{{0, 0, 0}};
    style::stroke black_stroke{black, 1, {3, 1, 3}};

    black_stroke.attach_attributes(stroked);

    std::cout << stroked << std::endl;
}

TEST(svg, transform)
{
    svg::object translated{"translated"};
    style::transform t0{{1., 2., 0.}};
    t0.attach_attributes(translated);
    std::cout << translated << std::endl;

    svg::object rotated{"rotated"};
    style::transform t1{{0., 0., 1.}};
    t1.attach_attributes(rotated);
    std::cout << rotated << std::endl;

    svg::object translated_rotated{"translated_rotated"};
    style::transform t2{{3., 2., 1.}};
    t2.attach_attributes(translated_rotated);
    std::cout << translated_rotated << std::endl;

    svg::object translated_rotated_scaled{"translated_rotated_scaled"};
    style::transform t3{{3., 2., 1.}};
    t3._scale = {100., 1.};
    t3.attach_attributes(translated_rotated_scaled);
    std::cout << translated_rotated_scaled << std::endl;

}
