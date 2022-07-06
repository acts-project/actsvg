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
#include "actsvg/data/odd_pixel_barrel.hpp"
#include "actsvg/meta.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

// Helper method to generate the barrel volume description
template <typename point3_container_type = point3_container>
proto::volume<point3_container_type> generate_barrel_volume() {

    proto::volume<point3_container> barrel;

    scalar z_min = std::numeric_limits<scalar>::max();
    scalar z_max = std::numeric_limits<scalar>::min();

    // ODD module templates
    proto::surface<point3_container> barrel_module_template;
    barrel_module_template._vertices = {
        {-8.4, -36, 0.}, {8.4, -36, 0.}, {8.4, 36., 0.}, {-8.4, 36., 0.}};
    barrel_module_template._measures = {8.4, 36.};
    barrel_module_template._type =
        proto::surface<point3_container>::e_rectangle;
    barrel._template_surfaces.push_back(barrel_module_template);

    size_t number_of_modules = data::odd_pixel_barrel.size() / 4u;
    barrel._surfaces.reserve(number_of_modules);
    for (size_t im = 0; im < number_of_modules; ++im) {
        // Register the templates - always first
        barrel._templates.push_back(0u);
        // Create the module surface
        proto::surface<point3_container> barrel_module;
        barrel_module._name = "module_" + std::to_string(im);
        barrel_module._vertices = {data::odd_pixel_barrel[4 * im],
                                   data::odd_pixel_barrel[4 * im + 1],
                                   data::odd_pixel_barrel[4 * im + 2],
                                   data::odd_pixel_barrel[4 * im + 3]};

        // Loop again for the z_min/z_max estimation
        for (size_t io = 0; io < 4; ++io) {
            scalar z = data::odd_pixel_barrel[4 * im + io][2];
            z_min = std::min(z_min, z);
            z_max = std::max(z_max, z);
        }

        // Add some descriptive text
        barrel_module._info = {"module #" + std::to_string(im),
                               display::center_string(barrel_module._vertices)};
        barrel._surfaces.push_back(barrel_module);
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
    scalar phi_step = 2 * M_PI / n_sectors;
    for (unsigned int is = 0; is <= n_sectors; ++is) {
        scalar c_phi = -M_PI + is * phi_step;
        phi_values.push_back(c_phi);
    }

    barrel._surface_grid._edges_0 = z_values;
    barrel._surface_grid._edges_1 = phi_values;

    // Create the associations by simple matching
    if (barrel._surface_grid._associations.empty()) {
        views::z_phi z_phi_view;

        for (unsigned int iz = 0; iz < z_tiles; ++iz) {
            scalar z_value = z_min + iz * z_step;
            for (unsigned int iphi = 0; iphi < n_sectors; ++iphi) {
                scalar phi_value = -M_PI + iphi * phi_step;

                std::cout << "bin " << iz << " x " << iphi << std::endl;

                std::map<unsigned long, unsigned long> module_associations;
                for (auto [is, s] : utils::enumerate(barrel._surfaces)) {
                    auto vertices = z_phi_view(s._vertices);
                    // Any touching vertex counts + central value
                    point2 center = {0., 0.};
                    for (auto v : vertices) {
                        center[0] += v[0];
                        center[1] += v[1];
                    }
                    center[0] /= vertices.size();
                    center[1] /= vertices.size();
                    vertices.push_back(center);

                    for (auto v : vertices) {
                        if (std::abs(z_value - v[0]) < 1.0 * z_step) {
                            scalar phi = v[1];
                            if (std::abs(phi - phi_value) < 1.0 * phi_step or
                                std::abs(phi - phi_value) >
                                    (2 * M_PI - phi_step)) {
                                std::cout
                                    << "- diff z =  " << z_value - v[0]
                                    << " / diff phi =  " << phi - phi_value
                                    << std::endl;
                                module_associations[is] = is;
                            }
                        }
                    }
                }
                std::vector<unsigned long> module_associtations_sl;
                for (auto [key, value] : module_associations) {
                    module_associtations_sl.push_back(key);
                }
                barrel._surface_grid._associations.push_back(
                    module_associtations_sl);
            }
        }
    }

    return barrel;
}

auto barrel = generate_barrel_volume<>();

TEST(display, barrel_sheet_module_info) {

    barrel._name = "ODD Pixel Barrel (sample)";

    // Create the sheet
    svg::object barrel_sheet = display::barrel_sheet(
        "sheet_odd_barrel", barrel, {600, 600}, display::e_module_info);

    svg::file barrel_file;
    barrel_file._width = 1000;
    barrel_file.add_object(barrel_sheet);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_odd_barrel_module_info.svg");
    eout << barrel_file;
    eout.close();
}

TEST(display, barrel_sheet_grid_info) {

    barrel._name = "ODD Pixel Barrel (sample)";

    // Create the sheet
    svg::object barrel_sheet = display::barrel_sheet(
        "sheet_odd_barrel", barrel, {600, 600}, display::e_grid_info);

    svg::file barrel_file;
    barrel_file._width = 1000;
    barrel_file.add_object(barrel_sheet);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_odd_barrel_grid_info.svg");
    eout << barrel_file;
    eout.close();
}
