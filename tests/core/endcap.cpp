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
#include "actsvg/data/odd_pixel_ec.hpp"

using namespace actsvg;

using trapezoid = std::array<std::array<scalar, 3u>, 4u>;

std::vector<trapezoid> generate_endcap_modules() {
    std::vector<trapezoid> modules;
    size_t number_of_modules = data::odd_pixel_ec.size() / 4u;
    modules.reserve(number_of_modules);
    for (size_t im = 0; im < number_of_modules; ++im) {
        modules.push_back(
            {data::odd_pixel_ec[4 * im], data::odd_pixel_ec[4 * im + 1],
             data::odd_pixel_ec[4 * im + 2], data::odd_pixel_ec[4 * im + 3]});
    }
    return modules;
}

auto endcap_modules = generate_endcap_modules();

TEST(endcap, z_r_view) {

    svg::file ec_file;

    // Draw the surfaces
    style::fill module_color{{28, 156, 168}};
    module_color._fc._opacity = 0.5;
    module_color._fc._highlight = {"mouseover", "mouseout"};
    module_color._fc._hl_rgb = {{255, 0, 0}};

    style::stroke stroke_color{{8, 76, 87}};
    stroke_color._width = 0.5;

    views::z_r z_r_view;

    std::vector<svg::object> modules;
    for (auto [m, ecm] : utils::enumerate(endcap_modules)) {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = z_r_view(ecm);
        modules.push_back(
            draw::polygon(m_id, module_contour, module_color, stroke_color));
    }

    // Add the surfaces
    ec_file._objects.insert(ec_file._objects.end(), modules.begin(),
                            modules.end());

    // File output
    std::ofstream ec_stream;
    ec_stream.open("odd_pixel_ec_zr.svg");
    ec_stream << ec_file;
    ec_stream.close();
}

TEST(endcap, x_y_view) {

    svg::file ec_file;
    ec_file._height = 800;
    ec_file._width = 800;

    // Draw the surfaces
    style::fill module_color{{28, 156, 168}};
    module_color._fc._opacity = 0.5;
    module_color._fc._hl_rgb = {{10, 200, 10}};
    module_color._fc._highlight = {"mouseover", "mouseout"};

    style::stroke stroke_color{{8, 76, 87}};
    stroke_color._width = 0.5;

    views::x_y x_y_view;

    std::vector<svg::object> modules;
    std::vector<svg::object> labels;
    for (auto [m, ecm] : utils::enumerate(endcap_modules)) {
        std::string m_id = std::string("m") + std::to_string(m);
        std::string t_id = std::string("t") + std::to_string(m);
        auto module_contour = x_y_view(ecm);
        auto module =
            draw::polygon(m_id, module_contour, module_color, stroke_color);
        modules.push_back(module);
        std::string module_txt = "Module " + std::to_string(m);
        std::string center_txt =
            "Center (" + std::to_string(module._real_barycenter[0]);
        center_txt += __c + std::to_string(module._real_barycenter[1]);
        center_txt += ")";
        std::vector<std::string> text = {"Module " + std::to_string(m),
                                         "Center "};

        auto ctext = draw::connected_text(
            t_id, module._real_barycenter, {module_txt, center_txt},
            style::font(), style::transform(), module);
        labels.push_back(ctext);
    }

    // Add the surfaces
    ec_file._objects.insert(ec_file._objects.end(), modules.begin(),
                            modules.end());
    ec_file._objects.insert(ec_file._objects.end(), labels.begin(),
                            labels.end());

    // File output
    std::ofstream ec_stream;
    ec_stream.open("odd_pixel_ec_xy.svg");
    ec_stream << ec_file;
    ec_stream.close();
}

TEST(endcap, x_y_view_grid) {

    svg::file ec_file;
    ec_file._height = 800;
    ec_file._width = 800;

    // Draw the surfaces
    style::fill module_color{{28, 156, 168}};
    module_color._fc._opacity = 0.5;
    module_color._fc._hl_rgb = {{10, 200, 10}};

    style::stroke stroke_color{{8, 76, 87}};
    stroke_color._width = 0.5;

    style::font font_style;
    font_style._size = 8;

    views::x_y x_y_view;

    std::vector<svg::object> modules;
    for (auto [m, ecm] : utils::enumerate(endcap_modules)) {
        std::string m_id = std::string("m") + std::to_string(m);
        auto module_contour = x_y_view(ecm);
        modules.push_back(
            draw::polygon(m_id, module_contour, module_color, stroke_color));
    }

    // Let's generate a grid & draw it
    style::fill grid_color{{200, 200, 200}};
    grid_color._fc._opacity = 0.25;
    grid_color._fc._highlight = {"mouseover", "mouseout"};
    grid_color._fc._hl_rgb = {{255, 0, 0}};

    style::stroke grid_stroke{{255, 0, 0}};
    grid_stroke._width = 0.5;
    grid_stroke._dasharray = {1, 1};

    std::vector<scalar> r_values = {42., 108., 174.};
    std::vector<scalar> phi_values;
    unsigned int n_sectors = 48;
    phi_values.reserve(n_sectors);
    scalar phi_step = 2 * M_PI / n_sectors;
    for (unsigned int is = 0; is <= n_sectors; ++is) {
        scalar c_phi = -M_PI + is * phi_step;
        phi_values.push_back(c_phi);
    }
    auto grid_sectors = draw::tiled_polar_grid("r_phi_",  r_values, phi_values,
                                         grid_color, grid_stroke);

    // Create the connection map here' simply with some tolerances
    scalar close_by_r = 75.;
    scalar close_by_phi = 0.1;

    std::vector<std::vector<size_t>> associations;
    for (auto [ig, g] : utils::enumerate(grid_sectors)) {
        std::vector<size_t> sector_associations;
        for (auto [is, s] : utils::enumerate(modules)) {
            // phi matching only
            scalar g_phi =
                std::atan2(g._real_barycenter[1], g._real_barycenter[0]);
            scalar s_phi =
                std::atan2(s._real_barycenter[1], s._real_barycenter[0]);
            if (std::abs(g_phi - s_phi) < 0.25 or
                std::abs(g_phi - s_phi) > (2 * M_PI - 0.25)) {
                sector_associations.push_back(is);
                std::cout << is << ",";
            }
        };
        associations.push_back(sector_associations);
        std::cout << std::endl;
    }

    // Build the connectors
    connectors::connect_objects(grid_sectors, modules, associations);

    // Add the surfaces
    ec_file._objects.insert(ec_file._objects.end(), modules.begin(),
                            modules.end());
    // Add the grid sectors
    ec_file._objects.insert(ec_file._objects.end(), grid_sectors.begin(),
                            grid_sectors.end());

    // File output
    std::ofstream ec_stream;
    ec_stream.open("odd_pixel_ec_grid_xy.svg");
    ec_stream << ec_file;
    ec_stream.close();
}
