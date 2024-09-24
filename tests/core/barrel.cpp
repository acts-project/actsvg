// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <array>
#include <fstream>
#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/data/odd_pixel_barrel.hpp"

using namespace actsvg;

auto barrel_modules = data::generate_barrel_modules();

TEST(barrel, x_y_view) {

    svg::file barrel_file;

    // Draw the surfaces
    style::fill module_color{style::color{{28, 156, 168}}};
    module_color._fc._opacity = 0.5;
    module_color._fc._hl_rgb = {{10, 200, 10}};

    style::stroke stroke_color{style::color{{8, 76, 87}}};
    stroke_color._width = 0.5;

    views::x_y x_y_view;

    std::vector<svg::object> modules;
    for (auto [m, bm] : utils::enumerate(barrel_modules)) {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = x_y_view.path(bm);
        modules.push_back(
            draw::polygon(m_id, module_contour, module_color, stroke_color));
    }

    // Add the surfaces
    barrel_file.add_objects(modules);

    // File output
    std::ofstream barrel_stream;
    barrel_stream.open("test_core_barrel_xy.svg");
    barrel_stream << barrel_file;
    barrel_stream.close();
}

TEST(barrel, z_phi_view) {

    svg::file barrel_file;

    // Draw the surfaces
    style::fill module_color{style::color{{28, 156, 168}}};
    module_color._fc._opacity = 0.5;
    module_color._fc._highlight = {"mouseover", "mouseout"};
    module_color._fc._hl_rgb = {{10, 200, 10}};

    style::stroke stroke_color{style::color{{8, 76, 87}}};
    stroke_color._width = 0.5;

    style::transform scale;
    scale._scale = {1, 150};

    views::z_phi z_phi_view;

    std::vector<svg::object> modules;
    for (auto [m, bm] : utils::enumerate(barrel_modules)) {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = z_phi_view.path(bm);
        modules.push_back(draw::polygon(m_id, module_contour, module_color,
                                        stroke_color, scale));
    }

    // Add the surfaces
    barrel_file.add_objects(modules);

    // File output
    std::ofstream barrel_stream;
    barrel_stream.open("test_core_barrel_zphi.svg");
    barrel_stream << barrel_file;
    barrel_stream.close();
}

TEST(barrel, z_phi_view_grid) {

    svg::file barrel_file;

    // Draw the surfaces
    style::fill module_color{style::color{{28, 156, 168}}};
    module_color._fc._opacity = 0.5;
    module_color._fc._highlight = {"mouseover", "mouseout"};
    module_color._fc._hl_rgb = {{10, 200, 10}};

    style::stroke stroke_color{style::color{{8, 76, 87}}};
    stroke_color._width = 0.5;

    style::transform scale;
    scale._scale = {1, 150};

    // Let's generate a grid & draw it
    style::fill grid_color{style::color{{200, 200, 200}}};
    grid_color._fc._opacity = 0.25;
    grid_color._fc._highlight = {"mouseover", "mouseout"};
    grid_color._fc._hl_rgb = {{255, 0, 0}};

    style::stroke grid_stroke{style::color{{255, 0, 0}}};
    grid_stroke._width = 0.5;
    grid_stroke._dasharray = {1, 1};

    views::z_phi z_phi_view;

    // Create the module objects
    std::vector<svg::object> modules;
    for (auto [m, bm] : utils::enumerate(barrel_modules)) {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = z_phi_view.path(bm);
        modules.push_back(draw::polygon(m_id, module_contour, module_color,
                                        stroke_color, scale));
    }

    // Find out the min/max of the z values
    scalar z_min = std::numeric_limits<scalar>::max();
    scalar z_max = std::numeric_limits<scalar>::min();
    for (auto& m : modules) {
        z_min = std::min(z_min, m._x_range[0]);
        z_max = std::max(z_max, m._x_range[1]);
    }

    // Grid construction: z values
    std::vector<scalar> z_values;
    unsigned int z_tiles = 14;
    z_values.reserve(z_tiles);
    scalar z_step = (z_max - z_min) / z_tiles;
    for (unsigned int iz = 0; iz <= z_tiles; ++iz) {
        scalar z_value = z_min + iz * z_step;
        z_values.push_back(z_value);
    }

    // Grid construction: phi values
    std::vector<scalar> phi_values;
    unsigned int n_sectors = 48;
    phi_values.reserve(n_sectors);
    scalar phi_step = 2_scalar * pi / n_sectors;
    for (unsigned int is = 0; is <= n_sectors; ++is) {
        scalar c_phi = -pi + is * phi_step;
        phi_values.push_back(c_phi);
    }

    // Construct the grid
    auto grid = draw::tiled_cartesian_grid("z_phi_", z_values, phi_values,
                                           grid_color, grid_stroke, scale);

    // Add the surfaces
    barrel_file.add_objects(modules);
    // Add the grid tiles
    barrel_file.add_objects(grid._sub_objects);

    // File output
    std::ofstream barrel_stream;
    barrel_stream.open("test_core_barrel_grid_zphi.svg");
    barrel_stream << barrel_file;
    barrel_stream.close();
}
