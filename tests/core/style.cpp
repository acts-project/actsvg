// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core/style.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include "actsvg/core/svg.hpp"

using namespace actsvg;

TEST(svg, fill_style) {

    svg::object colored{"red_object"};

    style::color red{{255, 0, 0}};
    style::fill red_fill{red};
    colored._fill = red_fill;

    std::cout << colored << std::endl;

    svg::object highlight{"red_highlight_object"};

    style::color red_hl{{255, 0, 0}};
    red_hl._hl_rgb = {0, 255, 0};
    red_hl._highlight = {"mouseover", "mouseout"};
    style::fill red_hl_fill{red_hl};
    highlight._fill = red_hl_fill;

    std::cout << highlight << std::endl;
}

TEST(svg, stroke_style) {

    svg::object stroked{"stroked_object"};

    style::color black{{0, 0, 0}};
    style::stroke black_stroke{black, 1, {3, 1, 3}};

    stroked._stroke = black_stroke;
    std::cout << stroked << std::endl;
}

TEST(svg, transform) {
    svg::object translated{"translated"};
    style::transform t0{{1., 2., 0.}};
    translated._transform = t0;
    std::cout << translated << std::endl;

    svg::object rotated{"rotated"};
    style::transform t1{{0., 0., 1.}};
    rotated._transform = t1;
    std::cout << rotated << std::endl;

    svg::object translated_rotated{"translated_rotated"};
    style::transform t2{{3., 2., 1.}};
    translated_rotated._transform = t2;
    std::cout << translated_rotated << std::endl;

    svg::object translated_rotated_scaled{"translated_rotated_scaled"};
    style::transform t3{{3., 2., 1.}};
    t3._scale = {100., 1.};
    translated_rotated_scaled._transform = t3;
    std::cout << translated_rotated_scaled << std::endl;
}
