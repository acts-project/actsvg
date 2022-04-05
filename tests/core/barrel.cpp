// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fstream>
#include <gtest/gtest.h>

#include "actsvg/actsvg.hpp"
#include "actsvg/data/odd_pixel_barrel.hpp"

#include <array>
#include <vector>
#include <string>
#include <fstream>

using namespace actsvg;

using rectangle = std::array<std::array<scalar, 3u>, 4u>;

std::vector<rectangle> generate_barrel_modules()
{
    std::vector<rectangle> modules;
    size_t number_of_modules = data::odd_pixel_barrel.size() / 4u;
    modules.reserve(number_of_modules);
    for (size_t im = 0; im < number_of_modules; ++im)
    {
        modules.push_back({data::odd_pixel_barrel[4 * im], data::odd_pixel_barrel[4 * im + 1],
                           data::odd_pixel_barrel[4 * im + 2], data::odd_pixel_barrel[4 * im + 3]});
    }
    return modules;
}

auto barrel_modules = generate_barrel_modules();

TEST(barrel, x_y_view)
{

    svg::file barrel_file;

    // Draw the surfaces
    style::fill module_color{{28, 156, 168}};
    module_color._fc._opacity = 0.5;
    module_color._fc._hl_rgb = {{10, 200, 10}};

    style::stroke stroke_color{{8, 76, 87}};
    stroke_color._width = 0.5;

    views::x_y x_y_view;

    std::vector<svg::object> modules;
    for (auto [m, bm] : utils::enumerate(barrel_modules))
    {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = x_y_view(bm);
        modules.push_back(draw::polygon(module_contour, m_id, module_color, stroke_color));
    }
   
    // Add the surfaces
    barrel_file._objects.insert(barrel_file._objects.end(), modules.begin(), modules.end());

    // File output
    std::ofstream barrel_stream;
    barrel_stream.open("odd_pixel_barrel_xy.svg");
    barrel_stream << barrel_file;
    barrel_stream.close();
}

TEST(barrel, z_phi_view)
{

    svg::file barrel_file;

    // Draw the surfaces
    style::fill module_color{{28, 156, 168}};
    module_color._fc._opacity = 0.5;
    module_color._fc._highlight = {"mouseover", "mouseout"};
    module_color._fc._hl_rgb = {{10, 200, 10}};

    style::stroke stroke_color{{8, 76, 87}};
    stroke_color._width = 0.5;

    style::transform scale;
    scale._scale = { 1, 150 };

    views::z_phi z_phi_view;

    std::vector<svg::object> modules;
    for (auto [m, bm] : utils::enumerate(barrel_modules))
    {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = z_phi_view(bm);
        modules.push_back(draw::polygon(module_contour, m_id, module_color, stroke_color, scale));
    }
   
    // Add the surfaces
    barrel_file._objects.insert(barrel_file._objects.end(), modules.begin(), modules.end());

    // File output
    std::ofstream barrel_stream;
    barrel_stream.open("odd_pixel_barrel_zphi.svg");
    barrel_stream << barrel_file;
    barrel_stream.close();
}

TEST(barrel, z_phi_view_grid)
{

    svg::file barrel_file;

    // Draw the surfaces
    style::fill module_color{{28, 156, 168}};
    module_color._fc._opacity = 0.5;
    module_color._fc._highlight = {"mouseover", "mouseout"};
    module_color._fc._hl_rgb = {{10, 200, 10}};

    style::stroke stroke_color{{8, 76, 87}};
    stroke_color._width = 0.5;

    style::transform scale;
    scale._scale = { 1, 150 };

    // Let's generate a grid & draw it
    style::fill grid_color{{200, 200, 200}};
    grid_color._fc._opacity = 0.25;
    grid_color._fc._highlight = {"mouseover", "mouseout"};
    grid_color._fc._hl_rgb = {{255, 0, 0}};

    style::stroke grid_stroke{{255, 0, 0}};
    grid_stroke._width = 0.5;
    grid_stroke._dasharray = {1, 1};


    views::z_phi z_phi_view;

    // Create the module objects 
    std::vector<svg::object> modules;
    for (auto [m, bm] : utils::enumerate(barrel_modules))
    {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = z_phi_view(bm);
        modules.push_back(draw::polygon(module_contour, m_id, module_color, stroke_color, scale));
    }

    // Find out the min/max of the z values 
    scalar z_min = std::numeric_limits<scalar>::max();
    scalar z_max = std::numeric_limits<scalar>::min();
    for (auto& m : modules){
        z_min = std::min(z_min, m._x_range[0]);
        z_max = std::max(z_max, m._x_range[1]);
    }

    // Grid construction: z values
    std::vector<scalar> z_values;
    unsigned int z_tiles = 14;
    z_values.reserve(z_tiles);
    scalar z_step = (z_max - z_min)/z_tiles;
    for (unsigned int iz = 0; iz <= z_tiles; ++iz)
    {
        scalar z_value = z_min + iz * z_step;
        z_values.push_back(z_value);
    }

    // Grid construction: phi values
    std::vector<scalar> phi_values;
    unsigned int n_sectors = 48;
    phi_values.reserve(n_sectors);
    scalar phi_step = 2 * M_PI / n_sectors;
    for (unsigned int is = 0; is <= n_sectors; ++is)
    {
        scalar c_phi = -M_PI + is * phi_step;
        phi_values.push_back(c_phi);
    }

    // Construct the grid
    auto grid_tiles = draw::z_phi_grid(z_values, phi_values, grid_color, grid_stroke, scale);

    // Make the connections

    // Add the surfaces
    barrel_file._objects.insert(barrel_file._objects.end(), modules.begin(), modules.end());
    // Add the grid tiles
    barrel_file._objects.insert(barrel_file._objects.end(), grid_tiles.begin(), grid_tiles.end());

    // File output
    std::ofstream barrel_stream;
    barrel_stream.open("odd_pixel_barrel_grid_zphi.svg");
    barrel_stream << barrel_file;
    barrel_stream.close();
}

