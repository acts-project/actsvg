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
#include <tuple>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

TEST(display, sheet_rectangle) {

    proto::surface<point3_container> rectangle;
    rectangle._type = proto::surface<point3_container>::type::e_rectangle;
    rectangle._vertices = {
        {-8.5, -34, 0.}, {8.5, -34, 0.}, {8.5, 34., 0.}, {-8.5, 34., 0.}};
    rectangle._measures = {8.5, 34.};

    svg::object surface_sheet =
        display::surface_sheet_xy("sheet_rectangle", rectangle);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("sheet_rectangle.svg");
    sout << surface_file;
    sout.close();
}

TEST(display, sheet_trapezoid) {

    proto::surface<point3_container> trapezoid;
    trapezoid._type = proto::surface<point3_container>::type::e_trapez;
    trapezoid._vertices = {
        {-8.5, -34, 0.}, {8.5, -34, 0.}, {14.5, 34., 0.}, {-14.5, 34., 0.}};
    trapezoid._measures = {8.5, 14.5, 34.};

    svg::object surface_sheet =
        display::surface_sheet_xy("sheet_trapezoid", trapezoid);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("sheet_trapezoid.svg");
    sout << surface_file;
    sout.close();
}

TEST(display, sheet_diamond) {

    proto::surface<point3_container> diamond;
    diamond._type = proto::surface<point3_container>::type::e_diamond;
    diamond._vertices = {{-8.5, -34, 0.}, {8.5, -34, 0.},   {20., 0., 0.},
                         {14.5, 20., 0.}, {-14.5, 20., 0.}, {-20., 0., 0.}};
    diamond._measures = {8.5, 20., 14.5, 34., 20.};

    svg::object surface_sheet =
        display::surface_sheet_xy("sheet_diamond", diamond);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("sheet_diamond.svg");
    sout << surface_file;
    sout.close();
}

TEST(discplay, sheet_polygon) {

    proto::surface<point3_container> polygon;
    polygon._type = proto::surface<point3_container>::type::e_polygon;
    polygon._vertices = {{-40., -40., 0.},
                         {60., -40., 0.},
                         {100., 80., 0.},
                         {0., 150, 0.},
                         {-100., 70., 0.}};

    for (const auto& v : polygon._vertices) {
        polygon._measures.push_back(v[0]);
        polygon._measures.push_back(v[1]);
    }

    svg::object surface_sheet =
        display::surface_sheet_xy("sheet_polygon", polygon);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("sheet_polygon.svg");
    sout << surface_file;
    sout.close();
}

TEST(display, sheet_full_disc) {

    proto::surface<point3_container> full_disc;
    full_disc._type = proto::surface<point3_container>::type::e_disc;
    full_disc._radii = {0., 30.};
    full_disc._measures = {30.};

    svg::object surface_sheet =
        display::surface_sheet_xy("sheet_full_disc", full_disc);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("sheet_full_disc.svg");
    sout << surface_file;
    sout.close();
}

TEST(display, sheet_full_ring) {

    proto::surface<point3_container> full_ring;
    full_ring._type = proto::surface<point3_container>::type::e_disc;
    full_ring._radii = {10., 30.};
    full_ring._measures = {10., 30.};

    svg::object surface_sheet =
        display::surface_sheet_xy("sheet_full_ring", full_ring);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("sheet_full_ring.svg");
    sout << surface_file;
    sout.close();
}

TEST(display, sheet_wedge) {

    proto::surface<point3_container> wedge;
    wedge._type = proto::surface<point3_container>::type::e_disc;
    wedge._radii = {0., 30.};
    wedge._opening = {-0.1_scalar, 0.7_scalar};
    wedge._measures = {0., 30., 0.8};

    svg::object surface_sheet =
        display::surface_sheet_xy("sheet_wedge", wedge, {600, 600}, false);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("sheet_wedge.svg");
    sout << surface_file;
    sout.close();
}

TEST(display, sheet_sector) {

    proto::surface<point3_container> sector;
    sector._type = proto::surface<point3_container>::type::e_disc;
    sector._radii = {10., 30.};
    sector._opening = {1.0, 1.7_scalar};
    sector._measures = {10., 30., 0.7};

    svg::object surface_sheet =
        display::surface_sheet_xy("sheet_sector", sector, {600, 600}, false);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("sheet_sector.svg");
    sout << surface_file;
    sout.close();
}

TEST(discplay, sheet_annulus) {

    proto::surface<point3_container> annulus;
    annulus._type = proto::surface<point3_container>::type::e_annulus;

    auto convert =
        [](const std::array<double, 3>& in) -> std::array<scalar, 3> {
        return {static_cast<scalar>(in[0]), static_cast<scalar>(in[1]),
                static_cast<scalar>(in[2])};
    };

    std::vector<std::array<double, 3>> raw = {
        {4.00364, 3.67003, 0.}, {3.82937, 4.28028, 0.}, {3.60204, 4.87282, 0.},
        {3.32341, 5.44303, 0.}, {2.99565, 5.9865, 0.},  {2.62131, 6.49899, 0.},
        {2.20329, 6.97652, 0.}, {1.74484, 7.41539, 0.}, {2.9259, 12.4347, 0.},
        {3.95183, 11.7812, 0.}, {4.90633, 11.0272, 0.}, {5.77959, 10.1804, 0.},
        {6.56264, 9.24955, 0.}, {7.24743, 8.24421, 0.}, {7.82692, 7.17472, 0.}};

    std::transform(raw.begin(), raw.end(),
                   std::back_inserter(annulus._vertices), convert);

    annulus._measures = {7.2, 12.0, 0.74195, 1.33970, 0., -3., 2.};

    // Focus on the sheet
    svg::object surface_sheet =
        display::surface_sheet_xy("sheet_annulus", annulus);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("sheet_annulus.svg");
    sout << surface_file;
    sout.close();
}
