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
#include "actsvg/data/odd_pixel_endcap.hpp"

using namespace actsvg;

auto endcap_modules = data::generate_endcap_modules();

TEST(endcap, z_r_view) {

    svg::file ec_file;

    // Draw the surfaces
    style::fill module_color{style::color{{28, 156, 168}}};
    module_color._fc._opacity = 0.5;
    module_color._fc._highlight = {"mouseover", "mouseout"};
    module_color._fc._hl_rgb = {{255, 0, 0}};

    style::stroke stroke_color{style::color{{8, 76, 87}}};
    stroke_color._width = 0.5;

    views::z_r z_r_view;

    std::vector<svg::object> modules;
    for (auto [m, ecm] : utils::enumerate(endcap_modules)) {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = z_r_view.path(ecm);
        modules.push_back(
            draw::polygon(m_id, module_contour, module_color, stroke_color));
    }

    // Add the surfaces
    ec_file.add_objects(modules);

    // File output
    std::ofstream ec_stream;
    ec_stream.open("test_core_endcap_zr.svg");
    ec_stream << ec_file;
    ec_stream.close();
}

TEST(endcap, x_y_view) {

    svg::file ec_file;
    ec_file._height = 800;
    ec_file._width = 800;

    // Draw the surfaces
    style::fill module_color{style::color{{28, 156, 168}}};
    module_color._fc._opacity = 0.5;
    module_color._fc._hl_rgb = {{10, 200, 10}};
    module_color._fc._highlight = {"mouseover", "mouseout"};

    style::stroke stroke_color{style::color{{8, 76, 87}}};
    stroke_color._width = 0.5;

    views::x_y x_y_view;

    std::vector<svg::object> modules;
    std::vector<svg::object> labels;
    for (auto [m, ecm] : utils::enumerate(endcap_modules)) {
        std::string m_id = std::string("m") + std::to_string(m);
        std::string t_id = std::string("t") + std::to_string(m);
        auto module_contour = x_y_view.path(ecm);
        auto module =
            draw::polygon(m_id, module_contour, module_color, stroke_color);
        modules.push_back(module);
        std::string module_txt = "Module " + std::to_string(m);
        std::string center_txt =
            "Center (" + utils::to_string(module._barycenter[0]);
        center_txt += __c + utils::to_string(module._barycenter[1]);
        center_txt += ")";
        auto ctext = draw::connected_text(
            t_id, {module._barycenter[0], module._barycenter[1]},
            {module_txt, center_txt}, style::font(), style::transform(),
            module);
        labels.push_back(ctext);
    }

    // Add the surfaces
    ec_file.add_objects(modules);
    ec_file.add_objects(labels);

    // File output
    std::ofstream ec_stream;
    ec_stream.open("test_core_endcap_xy.svg");
    ec_stream << ec_file;
    ec_stream.close();
}

TEST(endcap, x_y_view_grid) {

    svg::file ec_file;
    ec_file._height = 800;
    ec_file._width = 800;

    // Draw the surfaces
    style::fill module_color{style::color{{28, 156, 168}}};
    module_color._fc._opacity = 0.5;
    module_color._fc._hl_rgb = {{10, 200, 10}};

    style::stroke stroke_color{style::color{{8, 76, 87}}};
    stroke_color._width = 0.5;

    style::font font_style;
    font_style._size = 8;

    views::x_y x_y_view;

    std::vector<svg::object> modules;
    for (auto [m, ecm] : utils::enumerate(endcap_modules)) {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = x_y_view.path(ecm);
        modules.push_back(
            draw::polygon(m_id, module_contour, module_color, stroke_color));
    }

    // Let's generate a grid & draw it
    style::fill grid_color{style::color{{200, 200, 200}}};
    grid_color._fc._opacity = 0.25;
    grid_color._fc._highlight = {"mouseover", "mouseout"};
    grid_color._fc._hl_rgb = {{255, 0, 0}};

    style::stroke grid_stroke{style::color{{255, 0, 0}}};
    grid_stroke._width = 0.5;
    grid_stroke._dasharray = {1, 1};

    std::vector<scalar> r_values = {42., 108., 174.};
    std::vector<scalar> phi_values;
    unsigned int n_sectors = 48;
    phi_values.reserve(n_sectors);
    scalar phi_step = 2_scalar * pi / n_sectors;
    for (unsigned int is = 0; is <= n_sectors; ++is) {
        scalar c_phi = -pi + is * phi_step;
        phi_values.push_back(c_phi);
    }
    auto grid_sectors = draw::tiled_polar_grid("r_phi_", r_values, phi_values,
                                               grid_color, grid_stroke);

    // Create the connection map here' simply with some tolerances
    scalar close_by_phi = 0.25;

    std::vector<std::vector<size_t>> associations;
    for (auto [ig, g] : utils::enumerate(grid_sectors._sub_objects)) {
        std::vector<size_t> sector_associations;
        for (auto [is, s] : utils::enumerate(modules)) {
            // phi matching only
            scalar g_phi = std::atan2(g._barycenter[1], g._barycenter[0]);
            scalar s_phi = std::atan2(s._barycenter[1], s._barycenter[0]);
            if (std::abs(g_phi - s_phi) < close_by_phi or
                std::abs(g_phi - s_phi) > (2 * M_PI - close_by_phi)) {
                sector_associations.push_back(is);
            }
        };
        associations.push_back(sector_associations);
        std::cout << std::endl;
    }

    // Build the connectors
    auto c_objects = connectors::connect_action(grid_sectors._sub_objects,
                                                modules, associations);
    ec_file.add_objects(c_objects);

    // Add the surfaces
    ec_file.add_objects(modules);
    // Add the grid sectors
    ec_file.add_objects(grid_sectors._sub_objects);

    // File output
    std::ofstream ec_stream;
    ec_stream.open("test_core_endcap_grid_xy.svg");
    ec_stream << ec_file;
    ec_stream.close();
}
