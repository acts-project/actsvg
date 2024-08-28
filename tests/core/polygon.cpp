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

TEST(draw, triangle) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file tfile0;
    std::vector<std::array<scalar, 2u>> triangle = {
        {-100, -100}, {100, -100}, {0, 100}};

    auto tsvg0 = draw::polygon("t0", triangle);

    tfile0._objects.push_back(pg);
    tfile0._objects.push_back(tsvg0);

    std::ofstream tstream;
    tstream.open("test_core_triangle.svg");
    tstream << tfile0;
    tstream.close();
}

TEST(draw, triangle_shifted) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});
    svg::file tfile0;
    std::vector<std::array<scalar, 2u>> triangle = {
        {-100, -100}, {100, -100}, {0, 100}};

    style::transform shift;
    shift._tr = {50, 50};

    auto tsvg0 =
        draw::polygon("t0", triangle, style::fill(), style::stroke(), shift);
    tfile0._objects.push_back(pg);
    tfile0._objects.push_back(tsvg0);

    std::ofstream tstream;
    tstream.open("test_core_triangle_shifted.svg");
    tstream << tfile0;
    tstream.close();
}

TEST(draw, triangle_rotated) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});
    svg::file tfile0;
    std::vector<std::array<scalar, 2u>> triangle = {
        {-100, -100}, {100, -100}, {0, 100}};

    style::transform rotation;
    rotation._rot[0] = {45};

    auto tsvg0 =
        draw::polygon("t0", triangle, style::fill(), style::stroke(), rotation);
    tfile0._objects.push_back(pg);
    tfile0._objects.push_back(tsvg0);

    std::ofstream tstream;
    tstream.open("test_core_triangle_rotated.svg");
    tstream << tfile0;
    tstream.close();
}

TEST(draw, triangle_rotated_shifted) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});
    svg::file tfile0;
    std::vector<std::array<scalar, 2u>> triangle = {
        {-100, -100}, {100, -100}, {0, 100}};

    style::transform rotation_shifted;
    rotation_shifted._rot[0] = {45};
    rotation_shifted._tr = {-100, 100};

    auto tsvg0 = draw::polygon("t0", triangle, style::fill(), style::stroke(),
                               rotation_shifted);
    tfile0._objects.push_back(pg);
    tfile0._objects.push_back(tsvg0);

    std::ofstream tstream;
    tstream.open("test_core_triangle_rotated_shifted.svg");
    tstream << tfile0;
    tstream.close();
}

TEST(draw, triangle_scaled) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});
    svg::file tfile0;
    std::vector<std::array<scalar, 2u>> triangle = {
        {-100, -100}, {100, -100}, {0, 100}};

    style::transform shift;
    shift._scale = {2, 2};

    auto tsvg0 =
        draw::polygon("t0", triangle, style::fill(), style::stroke(), shift);
    tfile0._objects.push_back(pg);
    tfile0._objects.push_back(tsvg0);

    std::ofstream tstream;
    tstream.open("test_core_triangle_scaled.svg");
    tstream << tfile0;
    tstream.close();
}

TEST(draw, triangle_highligh) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    std::vector<std::array<scalar, 2u>> triangle = {
        {-10, -30}, {80, 70}, {10, 200}};

    style::fill red_fill_hl_green(style::color{{255, 0, 0}});
    red_fill_hl_green._fc._highlight = {"mouseover", "mouseout"};
    red_fill_hl_green._fc._hl_rgb = {{0, 255, 0}};

    svg::file tfile1;
    auto tsvg1 = draw::polygon("t1", triangle, red_fill_hl_green);
    tfile1.add_object(pg);
    tfile1.add_object(tsvg1);

    std::ofstream tstream;
    tstream.open("test_core_triangle_highlight.svg");
    tstream << tfile1;
    tstream.close();
}

TEST(draw, disc_sector) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    style::fill blue_fill(style::color{{0, 0, 200}});
    blue_fill._fc._opacity = 0.5;

    auto blue_sector_vertices =
        generators::sector_contour(120., 230., -0.15_scalar, 0.15_scalar);
    auto blue_sector =
        draw::polygon("blue_sector", blue_sector_vertices, blue_fill);

    style::fill red_fill(style::color{{200, 0, 0}});
    red_fill._fc._opacity = 0.5;
    auto red_sector_vertices = generators::sector_contour(60., 75., -2.0, 2.0);
    auto red_sector =
        draw::polygon("red_sector", red_sector_vertices, red_fill);

    style::fill green_fill(style::color{{0, 200, 0}});
    green_fill._fc._opacity = 0.5;

    auto green_sector_vertices =
        generators::sector_contour(30., 70., 2.0, -2.0);
    auto green_sector =
        draw::polygon("red_sector", green_sector_vertices, green_fill);

    svg::file tfile1;
    tfile1.add_object(pg);
    tfile1.add_object(blue_sector);
    tfile1.add_object(red_sector);
    tfile1.add_object(green_sector);

    std::ofstream tstream;
    tstream.open("test_core_sector.svg");
    tstream << tfile1;
    tstream.close();
}
