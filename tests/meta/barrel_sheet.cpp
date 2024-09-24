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

#include "../common/playground.hpp"
#include "actsvg/core.hpp"
#include "actsvg/data/odd_pixel_barrel.hpp"
#include "actsvg/meta.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

namespace {

/** Helper to generate a barrel
 *
 * @param t_ use template surfaces
 * @param b_s_ create also backside surfaces
 *
 **/
proto::volume<point3_container> generate_barrel(
    bool t_ = false, bool b_s_ = false) noexcept(false) {

    // Create the endcap volume
    proto::volume<point3_container> barrel;
    barrel._name = "Barrel";
    barrel._surface_grid._type = proto::grid::e_z_phi;

    scalar hx = 8.4_scalar;
    scalar hy = 36.;
    scalar barrel_r = 32.;
    scalar tilt_phi = 0.14_scalar;
    int staves = 16;
    int stave_modules = 14;

    // The template
    proto::surface<point3_container> barrel_module_template;
    barrel_module_template._vertices = {
        {-hx, -hy, 0.}, {hx, -hy, 0.}, {hx, hy, 0.}, {-hx, hy, 0.}};
    barrel_module_template._measures = {hx, hy};
    barrel_module_template._type =
        proto::surface<point3_container>::type::e_rectangle;
    barrel_module_template._sf_type =
        proto::surface<point3_container>::sf_type::e_sensitive;

    scalar offz =
        (0.5_scalar * stave_modules - 0.5_scalar) * (1.98_scalar * hy);
    scalar nexz = offz / (0.5_scalar * stave_modules - 0.5_scalar);

    // A view for the template
    views::z_phi z_phi_view;
    views::x_y x_y_view;

    scalar delta_phi = 2 * pi / staves;

    /// The surface containers
    std::vector<proto::surface<point3_container>> regular;
    std::vector<proto::surface<point3_container>> backside;

    /// The bin association vector
    std::vector<std::vector<size_t>> bin_associations;

    size_t cmodule = 0;
    for (int iz = 0; iz < stave_modules; ++iz) {
        scalar tr_z = -offz + iz * nexz;
        if (iz == 0) {
            barrel._surface_grid._edges_0.push_back(tr_z - hy);
        }
        barrel._surface_grid._edges_0.push_back(tr_z + hy);

        for (int is = 0; is < staves; ++is, cmodule++) {
            scalar cphi = -pi + is * delta_phi;
            if (iz == 0) {
                if (is == 0) {
                    barrel._surface_grid._edges_1.push_back(
                        cphi + 0.5_scalar * delta_phi);
                }
                barrel._surface_grid._edges_1.push_back(cphi +
                                                        1.5_scalar * delta_phi);
            }

            // The position of the module
            scalar tr_x = barrel_r * std::cos(cphi);
            scalar tr_y = barrel_r * std::sin(cphi);
            std::array<scalar, 3> tr = {tr_x, tr_y, tr_z};

            std::vector<size_t> barrel_assoc = {cmodule, cmodule + 1,
                                                cmodule + 2};
            bin_associations.push_back(barrel_assoc);

            // The columns of the rotation matrix
            scalar mphi = 0.5_scalar * pi + cphi + tilt_phi;
            std::array<scalar, 3> col_x = {std::cos(mphi), std::sin(mphi), 0.};
            std::array<scalar, 3> col_y = {0., 0., 1.};
            if (not t_) {

                point3 ll =
                    utils::add(tr, utils::add(utils::scale(col_y, -hy),
                                              utils::scale(col_x, -hx)));
                point3 lr = utils::add(tr, utils::add(utils::scale(col_y, -hy),
                                                      utils::scale(col_x, hx)));
                point3 hr = utils::add(tr, utils::add(utils::scale(col_y, hy),
                                                      utils::scale(col_x, hx)));
                point3 hl =
                    utils::add(tr, utils::add(utils::scale(col_y, hy),
                                              utils::scale(col_x, -hx)));

                proto::surface<point3_container> barrel_module;
                barrel_module._name =
                    "module_" + std::to_string(is) + "_" + std::to_string(iz);
                barrel_module._vertices = {ll, lr, hr, hl};
                barrel_module._measures = barrel_module_template._measures;
                barrel_module._type = barrel_module_template._type;
                barrel_module._sf_type = barrel_module_template._sf_type;
                barrel_module._aux_info["module_info"] = {
                    std::string("Module " + std::to_string(is) + " " +
                                std::to_string(iz)),
                    std::string("center z/phi = " + utils::to_string(tr_z) +
                                "/" + std::to_string(cphi))};
                barrel_module._aux_info["grid_info"] = {
                    std::string("* module " + std::to_string(iz) + " " +
                                std::to_string(is))};

                barrel_module._template_object = display::surface(
                    barrel_module._name + "_template", barrel_module_template,
                    x_y_view, {true, true, true, true});

                regular.push_back(barrel_module);

                // Create a backside module
                if (b_s_) {
                    // The position of the module
                    tr_x = (barrel_r + 2._scalar) * std::cos(cphi);
                    tr_y = (barrel_r + 2._scalar) * std::sin(cphi);
                    tr = {tr_x, tr_y, tr_z};

                    ll = utils::add(tr, utils::add(utils::scale(col_y, -hy),
                                                   utils::scale(col_x, -hx)));
                    lr = utils::add(tr, utils::add(utils::scale(col_y, -hy),
                                                   utils::scale(col_x, hx)));
                    hr = utils::add(tr, utils::add(utils::scale(col_y, hy),
                                                   utils::scale(col_x, hx)));
                    hl = utils::add(tr, utils::add(utils::scale(col_y, hy),
                                                   utils::scale(col_x, -hx)));

                    proto::surface<point3_container> barrel_module_backside;
                    barrel_module_backside._name = "module_backside_" +
                                                   std::to_string(is) + "_" +
                                                   std::to_string(iz);

                    barrel_module_backside._vertices = {ll, lr, hr, hl};
                    barrel_module_backside._measures =
                        barrel_module_template._measures;
                    barrel_module_backside._type = barrel_module_template._type;
                    barrel_module_backside._sf_type =
                        barrel_module_template._sf_type;
                    barrel_module_backside._sf_type =
                        barrel_module_template._sf_type;
                    barrel_module_backside._aux_info["module_info"] = {
                        std::string("Module (backside) " + std::to_string(is) +
                                    " " + std::to_string(iz)),
                        std::string("center z/phi = " + utils::to_string(tr_z) +
                                    "/" + std::to_string(cphi))};
                    barrel_module_backside._aux_info["grid_info"] = {
                        std::string("* module " + std::to_string(iz) + " " +
                                    std::to_string(is))};

                    barrel_module_backside._template_object =
                        display::surface(barrel_module._name + "_template",
                                         barrel_module_template, x_y_view,
                                         {true, true, true, true});

                    backside.push_back(barrel_module_backside);
                }
            }
        }
    }

    // Store the surfaces
    barrel._surfaces.push_back(regular);
    barrel._grid_associations.push_back(bin_associations);
    if (not backside.empty()) {
        barrel._surfaces.push_back(backside);
        barrel._grid_associations.push_back(bin_associations);
    }
    return barrel;
}
}  // namespace

TEST(display, barrel_x_y_view) {

    auto barrel = generate_barrel();

    // Set a playground
    auto pg = test::playground({-100, -100}, {100, 100});

    svg::file barrel_file;
    barrel_file._width = 1000;

    // Draw the surfaces
    style::fill module_color{style::color{{28, 156, 168}}};
    module_color._fc._opacity = 0.5;
    module_color._fc._hl_rgb = {{10, 200, 10}};

    style::stroke stroke_color{style::color{{8, 76, 87}}};
    stroke_color._width = 0.5;

    views::x_y x_y_view;

    std::vector<svg::object> modules;
    for (auto [s, surfaces] : utils::enumerate(barrel._surfaces)) {
        for (auto [m, bm] : utils::enumerate(surfaces)) {
            std::string m_id =
                std::string("m_") + std::to_string(s) + "_" + std::to_string(m);
            auto module_contour = x_y_view.path(bm._vertices);
            modules.push_back(draw::polygon(m_id, module_contour, module_color,
                                            stroke_color));
        }
    }
    // Add the surfaces
    barrel_file.add_objects(modules);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_barrel_xy_view.svg");
    eout << barrel_file;
    eout.close();
}

TEST(display, barrel_sheet_module_info) {

    auto barrel = generate_barrel();

    // Create the sheet
    svg::object barrel_sheet = display::barrel_sheet(
        "barrel_sheet", barrel, {600, 600}, display::e_module_info);

    svg::file barrel_file_fs;
    barrel_file_fs._width = 1000;
    barrel_file_fs.add_object(barrel_sheet);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_barrel_module_info.svg");
    eout << barrel_file_fs;
    eout.close();
}

TEST(display, barrel_sheet_grid_info) {

    auto barrel = generate_barrel();

    // Create the sheet
    svg::object barrel_sheet = display::barrel_sheet(
        "barrel_sheet", barrel, {600, 600}, display::e_grid_info);

    svg::file barrel_file_fs;
    barrel_file_fs._width = 1000;
    barrel_file_fs.add_object(barrel_sheet);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_barrel_grid_info.svg");
    eout << barrel_file_fs;
    eout.close();
}

TEST(display, barrel_sheet_module_info_backside) {

    auto barrel = generate_barrel(true);

    // Create the sheet
    svg::object barrel_sheet = display::barrel_sheet(
        "barrel_sheet_backside", barrel, {600, 600}, display::e_module_info);

    svg::file barrel_file_fs;
    barrel_file_fs._width = 1000;
    barrel_file_fs.add_object(barrel_sheet);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_barrel_module_info_backside.svg");
    eout << barrel_file_fs;
    eout.close();
}
