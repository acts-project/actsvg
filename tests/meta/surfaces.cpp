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
using point3_container = std::vector<point3>;

TEST(proto, rectanglular_surface) {

    point3_container rectangle_vertices = {{-30., -70., 123.},
                                           {30., -70., 123.},
                                           {30., 70., 123.},
                                           {-30., 70., 123.}};

    proto::surface<point3_container> rectangle;
    rectangle._vertices = rectangle_vertices;
    rectangle._name = "rectangle surface";
    rectangle._type = proto::surface<point3_container>::type::e_rectangle;
    rectangle._fill = style::fill({{0, 0, 100}, 0.5});

    svg::object rec = display::surface("rectangle", rectangle, views::x_y{});

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file rfile;
    rfile.add_object(pg);
    rfile.add_object(rec);

    std::ofstream rstream;
    rstream.open("test_meta_rectangle.svg");
    rstream << rfile;
    rstream.close();
}

TEST(proto, rectanglular_subtracted_surface) {

    point3_container subtracted_vertices = {{-30., -70., 123.},
                                            {30., -70., 123.},
                                            {30., 70., 123.},
                                            {-30., 70., 123.}};

    point3_container rectangle_vertices = {{-130., -170., 123.},
                                           {130., -170., 123.},
                                           {130., 170., 123.},
                                           {-130., 170., 123.}};

    proto::surface<point3_container> rectangle;
    rectangle._vertices = rectangle_vertices;
    rectangle._name = "rectangle surface";
    rectangle._type = proto::surface<point3_container>::type::e_rectangle;
    rectangle._fill = style::fill({{0, 0, 100}, 0.5});

    proto::surface<point3_container> subtracted_rectangle;
    subtracted_rectangle._vertices = subtracted_vertices;

    rectangle._boolean_surface = {subtracted_rectangle};
    rectangle._boolean_operation =
        proto::surface<point3_container>::boolean::e_subtraction;

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

TEST(proto, rectanglular_surfaces_at_phi_poles) {

    scalar d_phi = 0.23_scalar;
    point3_container around_0_vertices = {
        {70., static_cast<scalar>(70. * std::sin(d_phi)), -100.},
        {70., static_cast<scalar>(70. * std::sin(-d_phi)), -100.},
        {70., static_cast<scalar>(70. * std::sin(-d_phi)), 100.},
        {70., static_cast<scalar>(70. * std::sin(d_phi)), 100.}};

    proto::surface<point3_container> rectangle_0;
    rectangle_0._vertices = around_0_vertices;
    rectangle_0._name = "rectangle surface around phi = 0";
    rectangle_0._type = proto::surface<point3_container>::type::e_rectangle;
    rectangle_0._fill = style::fill({{0, 0, 100}, 0.5});
    rectangle_0._radii[0u] = 70.;

    point3_container around_pi_vertices = {
        {-70., static_cast<scalar>(70. * std::sin(-M_PI + d_phi)), -200.},
        {-70., static_cast<scalar>(70. * std::sin(M_PI - d_phi)), -200.},
        {-70., static_cast<scalar>(70. * std::sin(M_PI - d_phi)), -100.},
        {-70., static_cast<scalar>(70. * std::sin(-M_PI + d_phi)), -100.}};

    proto::surface<point3_container> rectangle_pi;
    rectangle_pi._vertices = around_pi_vertices;
    rectangle_pi._name = "rectangle surface around phi = M_PI";
    rectangle_pi._type = proto::surface<point3_container>::type::e_rectangle;
    rectangle_pi._fill = style::fill({{0, 100, 0}, 0.5});
    rectangle_pi._radii[0u] = 70.;

    svg::object rec_0_xy =
        display::surface("rectangle", rectangle_0, views::x_y{});
    svg::object rec_pi_xy =
        display::surface("rectangle", rectangle_pi, views::x_y{});

    svg::object rec_0_zphi = display::surface(
        "rectangle", rectangle_0, views::z_rphi{views::scene(), 70.});

    svg::object rec_pi_zphi = display::surface(
        "rectangle", rectangle_pi, views::z_rphi{views::scene(), 70.});

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::file rfile_xy;
    rfile_xy.add_object(pg);
    rfile_xy.add_object(rec_0_xy);
    rfile_xy.add_object(rec_pi_xy);

    svg::file rfile_zphi;
    rfile_zphi.add_object(pg);
    rfile_zphi.add_object(rec_0_zphi);
    rfile_zphi.add_object(rec_pi_zphi);

    std::ofstream rstream;
    rstream.open("test_meta_rectangle_phi_xy.svg");
    rstream << rfile_xy;
    rstream.close();

    rstream.open("test_meta_rectangle_phi_zphi.svg");
    rstream << rfile_zphi;
    rstream.close();
}

TEST(proto, full_disc) {

    proto::surface<point3_container> disc;
    disc._radii = {0., 150.};
    disc._zparameters = {50., 0.};
    disc._name = "disc surface";
    disc._type = proto::surface<point3_container>::type::e_disc;
    disc._fill = style::fill({{0, 100, 0}, 0.5});
    disc._stroke = style::stroke({{0, 0, 0}}, 2.);
    disc._stroke._hl_width = 2.;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Test the disc in x-y view
    svg::object fd_xy = display::surface("disc", disc, views::x_y{});

    svg::file rfile_xy;
    rfile_xy.add_object(pg);
    rfile_xy.add_object(fd_xy);

    std::ofstream rstream;
    rstream.open("test_meta_full_disc_xy.svg");
    rstream << rfile_xy;
    rstream.close();

    // Test the disc in z-r view
    svg::object fd_zr = display::surface("disc", disc, views::z_r{});

    svg::file rfile_zr;
    rfile_zr.add_object(pg);
    rfile_zr.add_object(fd_zr);

    rstream.open("test_meta_full_disc_zr.svg");
    rstream << rfile_zr;
    rstream.close();
}

TEST(proto, full_discs_view_range) {

    proto::surface<point3_container> disc;
    disc._radii = {0., 150.};
    disc._zparameters = {50., 0.};
    disc._name = "disc surface";
    disc._type = proto::surface<point3_container>::type::e_disc;
    disc._fill = style::fill({{0, 100, 0}, 0.5});
    disc._stroke = style::stroke({{0, 0, 0}}, 2.);
    disc._stroke._hl_width = 2.;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    views::x_y xy_range;
    xy_range._scene._range[1] = {0., 10.};

    // Test the disc in x-y view
    svg::object fd_xy = display::surface("disc_view_hidden", disc, xy_range);

    svg::file rfile_xy;
    rfile_xy.add_object(pg);
    rfile_xy.add_object(fd_xy);

    std::ofstream rstream;
    rstream.open("test_meta_full_view_hidden_disc_xy.svg");
    rstream << rfile_xy;
    rstream.close();
}

TEST(proto, full_ring) {

    proto::surface<point3_container> ring;
    ring._radii = {50., 150.};
    ring._zparameters = {50., 0.};
    ring._name = "ring surface";
    ring._type = proto::surface<point3_container>::type::e_disc;
    ring._fill = style::fill({{0, 100, 0}, 0.5});
    ring._stroke = style::stroke({{0, 0, 0}}, 2.);
    ring._stroke._hl_width = 2.;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    svg::object fr_xy = display::surface("ring", ring, views::x_y{});

    svg::file rfile_xy;
    rfile_xy.add_object(pg);
    rfile_xy.add_object(fr_xy);

    std::ofstream rstream;
    rstream.open("test_meta_full_ring_xy.svg");
    rstream << rfile_xy;
    rstream.close();

    svg::object fr_zr = display::surface("ring", ring, views::z_r{});

    svg::file rfile_zr;
    rfile_zr.add_object(pg);
    rfile_zr.add_object(fr_zr);

    rstream.open("test_meta_full_ring_zr.svg");
    rstream << rfile_zr;
    rstream.close();
}

TEST(proto, wedge) {

    proto::surface<point3_container> wedge;
    wedge._radii = {0., 150.};
    wedge._opening = {-0.25, 0.25};
    wedge._zparameters = {50., 0.};
    wedge._name = "wedge surface";
    wedge._type = proto::surface<point3_container>::type::e_disc;
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

    proto::surface<point3_container> sector;
    sector._radii = {50., 150.};
    sector._opening = {-0.25, 0.25};
    sector._zparameters = {50., 0.};
    sector._name = "sector surface";
    sector._type = proto::surface<point3_container>::type::e_disc;
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

TEST(proto, full_cylinder) {

    proto::surface<point3_container> cylinder;
    cylinder._radii = {0., 150.};
    cylinder._name = "cylinder surface";
    cylinder._type = proto::surface<point3_container>::type::e_cylinder;
    cylinder._fill = style::fill({{0, 100, 0}, 0.5});
    cylinder._stroke = style::stroke({{0, 0, 0}}, 2.);
    cylinder._stroke._hl_width = 2.;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Full cylinder in x-y
    svg::object fc_xy = display::surface("cylinder", cylinder, views::x_y{});

    svg::file rfile_xy;
    rfile_xy.add_object(pg);
    rfile_xy.add_object(fc_xy);

    std::ofstream rstream;
    rstream.open("test_meta_full_cylinder_xy.svg");
    rstream << rfile_xy;
    rstream.close();

    // Full cylinder in z-r
    svg::object fc_zr = display::surface("cylinder", cylinder, views::z_r{});

    svg::file rfile_zr;
    rfile_zr.add_object(pg);
    rfile_zr.add_object(fc_zr);

    rstream.open("test_meta_full_cylinder_zr.svg");
    rstream << rfile_zr;
    rstream.close();
}

TEST(proto, sector_cylinder) {

    proto::surface<point3_container> sector_cylinder;
    sector_cylinder._radii = {0., 150.};
    sector_cylinder._opening = {-0.2_scalar, 0.5};
    sector_cylinder._zparameters = {10., 200};
    sector_cylinder._name = "sectoral cylinder surface";
    sector_cylinder._type = proto::surface<point3_container>::type::e_cylinder;
    sector_cylinder._fill = style::fill({{0, 100, 0}, 0.5});
    sector_cylinder._stroke = style::stroke({{0, 0, 0}}, 2.);
    sector_cylinder._stroke._hl_width = 2.;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Full cylinder in x-y
    svg::object sc_xy =
        display::surface("sectoral_cylinder", sector_cylinder, views::x_y{});

    svg::file rfile_xy;
    rfile_xy.add_object(pg);
    rfile_xy.add_object(sc_xy);

    std::ofstream rstream;
    rstream.open("test_meta_sectoral_cylinder_xy.svg");
    rstream << rfile_xy;
    rstream.close();

    // Full cylinder in z-r
    svg::object sc_zr =
        display::surface("sectoral_cylinder", sector_cylinder, views::z_r{});

    svg::file rfile_zr;
    rfile_zr.add_object(pg);
    rfile_zr.add_object(sc_zr);

    rstream.open("test_meta_sectoral_cylinder_zr.svg");
    rstream << rfile_zr;
    rstream.close();
}

TEST(proto, straw_surface) {

    proto::surface<point3_container> straw;
    straw._radii = {2., 20.};
    straw._zparameters = {0., 200};
    straw._name = "straw surface";
    straw._type = proto::surface<point3_container>::type::e_straw;
    straw._fill = style::fill({{0, 100, 0}, 0.5});
    straw._stroke = style::stroke({{0, 0, 0}}, 1.);
    straw._stroke._hl_width = 2.;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Straw surface in x-y
    svg::object sc_xy = display::surface("straw", straw, views::x_y{});

    svg::file rfile_xy;
    rfile_xy.add_object(pg);
    rfile_xy.add_object(sc_xy);

    std::ofstream rstream;
    rstream.open("test_meta_straw_xy.svg");
    rstream << rfile_xy;
    rstream.close();
}

TEST(proto, annulus_surface) {

    scalar min_radius = 120.;
    scalar max_radius = 190.0;

    scalar min_phi = 0.74195_scalar;
    scalar max_phi = 1.33970_scalar;
    scalar avg_phi = 0.5_scalar * (min_phi + max_phi);

    scalar center_x = -20.;
    scalar center_y = 20.;

    // construct the annulus
    proto::surface<point3_container> annulus;
    annulus._measures = {min_radius, max_radius, min_phi, max_phi,
                         avg_phi,    center_x,   center_y};
    annulus._type = proto::surface<point3_container>::type::e_annulus;
    annulus._fill = style::fill({{0, 100, 0}, 0.5});
    annulus._stroke = style::stroke({{0, 0, 0}}, 2.);
    annulus._stroke._hl_width = 2.;

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Annulus in x-y
    svg::object an_xy = display::surface("annulus", annulus, views::x_y{});

    svg::file rfile_xy;
    rfile_xy.add_object(pg);
    rfile_xy.add_object(an_xy);

    std::ofstream rstream;
    rstream.open("test_meta_annulus_xy.svg");
    rstream << rfile_xy;
    rstream.close();
}
