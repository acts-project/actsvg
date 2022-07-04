// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "../common/playground.hpp"
#include "actsvg/core.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/proto/surface.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_collection = std::vector<point3>;

TEST(proto, rectanglular_surface) {

    point3_collection rectangle_vertices = {{-30., -70., 123.},
                                            {30., -70., 123.},
                                            {30., 70., 123.},
                                            {-30., 70., 123.}};

    proto::surface<point3_collection> rectangle;
    rectangle._vertices = rectangle_vertices;
    rectangle._name = "rectangle surface";
    rectangle._type = proto::surface<point3_collection>::e_rectangle;
    rectangle._fill = style::fill({{0, 0, 100}, 0.5});

    svg::object rec = display::surface("rectangle", rectangle, views::x_y{});

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file rfile;
    rfile.add_object(pg);
    // standard markers
    rfile.add_object(rec);

    std::ofstream rstream;
    rstream.open("test_meta_rectangle.svg");
    rstream << rfile;
    rstream.close();
}

TEST(proto, rectanglular_subtracted_surface) {

    point3_collection subtracted_vertices = {{-30., -70., 123.},
                                             {30., -70., 123.},
                                             {30., 70., 123.},
                                             {-30., 70., 123.}};

    point3_collection rectangle_vertices = {{-130., -170., 123.},
                                            {130., -170., 123.},
                                            {130., 170., 123.},
                                            {-130., 170., 123.}};

    proto::surface<point3_collection> rectangle;
    rectangle._vertices = rectangle_vertices;
    rectangle._name = "rectangle surface";
    rectangle._type = proto::surface<point3_collection>::e_rectangle;
    rectangle._fill = style::fill({{0, 0, 100}, 0.5});

    proto::surface<point3_collection> subtracted_rectangle;
    subtracted_rectangle._vertices = subtracted_vertices;

    rectangle._boolean_surface = { subtracted_rectangle };
    rectangle._boolean_operation =
        proto::surface<point3_collection>::e_subtraction;

    svg::object rec = display::surface("rectangle", rectangle, views::x_y{});

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file rfile;
    rfile.add_object(pg);
    rfile.add_object(rec);

    std::ofstream rstream;
    rstream.open("test_meta_subtracted_rectangle.svg");
    rstream << rfile;
    rstream.close();
}


TEST(proto, full_disc) {

    proto::surface<point3_collection> disc;
    disc._radii = {0., 150.};
    disc._name = "disc surface";
    disc._type = proto::surface<point3_collection>::e_disc;
    disc._fill = style::fill({{0, 100, 0}, 0.5});

    svg::object fd = display::surface("disc", disc, views::x_y{});

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file rfile;
    rfile.add_object(pg);
    rfile.add_object(fd);

    std::ofstream rstream;
    rstream.open("test_meta_full_disc.svg");
    rstream << rfile;
    rstream.close();
}

TEST(proto, full_ring) {

    proto::surface<point3_collection> ring;
    ring._radii = {50., 150.};
    ring._name = "ring surface";
    ring._type = proto::surface<point3_collection>::e_disc;
    ring._fill = style::fill({{0, 100, 0}, 0.5});

    svg::object fr = display::surface("ring", ring, views::x_y{});

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file rfile;
    rfile.add_object(pg);
    rfile.add_object(fr);

    std::ofstream rstream;
    rstream.open("test_meta_full_ring.svg");
    rstream << rfile;
    rstream.close();
}

TEST(proto, wedge) {

    proto::surface<point3_collection> wedge;
    wedge._radii = {0., 150.};
    wedge._opening = {-0.25, 0.25};
    wedge._name = "wedge surface";
    wedge._type = proto::surface<point3_collection>::e_disc;
    wedge._fill = style::fill({{0, 100, 0}, 0.5});

    svg::object w = display::surface("wedge", wedge, views::x_y{});

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file rfile;
    rfile.add_object(pg);
    rfile.add_object(w);

    std::ofstream rstream;
    rstream.open("test_meta_wedge.svg");
    rstream << rfile;
    rstream.close();
}


TEST(proto, sector) {

    proto::surface<point3_collection> sector;
    sector._radii = {50., 150.};
    sector._opening = {-0.25, 0.25};
    sector._name = "sector surface";
    sector._type = proto::surface<point3_collection>::e_disc;
    sector._fill = style::fill({{0, 100, 0}, 0.5});

    svg::object s = display::surface("sector", sector, views::x_y{});

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file rfile;
    rfile.add_object(pg);
    rfile.add_object(s);

    std::ofstream rstream;
    rstream.open("test_meta_sector.svg");
    rstream << rfile;
    rstream.close();
}

