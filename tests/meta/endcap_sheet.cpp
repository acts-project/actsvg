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
#include "actsvg/data/odd_pixel_ec.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

// Helper method to generate the endcap volume description
template <typename point3_container_type = point3_container>
proto::volume<point3_container_type> generate_endcap_volume() {

    proto::volume<point3_container> endcap;

    // ODD module templates
    proto::surface<point3_container> endcap_inner_template;
    endcap_inner_template._vertices = {
        {-8.5, -34, 0.}, {8.5, -34, 0.}, {14.5, 34., 0.}, {-14.5, 34., 0.}};
    endcap_inner_template._measures = {8.5, 14.5, 34.};
    endcap._template_surfaces.push_back(endcap_inner_template);

    proto::surface<point3_container> endcap_outer_template;
    endcap_outer_template._vertices = {
        {-10.5, -34., 0.}, {10.5, -34., 0.}, {16.5, 34., 0.}, {-16.5, 34., 0}};
    endcap_outer_template._measures = {10.5, 16.5, 34.};
    endcap._template_surfaces.push_back(endcap_outer_template);

    size_t number_of_modules = data::odd_pixel_ec.size() / 4u;
    endcap._surfaces.reserve(number_of_modules);
    for (size_t im = 0; im < number_of_modules; ++im) {
        // Register the templates - first 24 are of template 0
        if (im < 24) {
            endcap._templates.push_back(0u);
        } else {
            endcap._templates.push_back(1u);
        }
        // Create the module surface
        proto::surface<point3_container> endcap_module;
        endcap_module._name = "module_" + std::to_string(im);
        endcap_module._vertices = {
            data::odd_pixel_ec[4 * im], data::odd_pixel_ec[4 * im + 1],
            data::odd_pixel_ec[4 * im + 2], data::odd_pixel_ec[4 * im + 3]};

        // Add some descriptive text
        endcap_module._info = {"module #" + std::to_string(im),
                               display::center_string(endcap_module._vertices)};
        endcap._surfaces.push_back(endcap_module);
    }

    // Let's create the grid
    std::vector<scalar> r_values = {42., 108., 174.};
    std::vector<scalar> phi_values;
    unsigned int n_sectors = 48;
    phi_values.reserve(n_sectors);
    scalar phi_step = 2 * M_PI / n_sectors;
    for (unsigned int is = 0; is <= n_sectors; ++is) {
        scalar c_phi = -M_PI + is * phi_step;
        phi_values.push_back(c_phi);
    }
    endcap._surface_grid._edges_0 = r_values;
    endcap._surface_grid._edges_1 = phi_values;

    endcap._surface_grid._associations = data::odd_pixel_ec_assoc;

    return endcap;
}

auto endcap = generate_endcap_volume<>();

TEST(display, endcap_sheet_module_info) {

    endcap._name = "ODD Pixel Endcap (sample)";

    // Create the sheet
    svg::object endcap_sheet =
        display::endcap_sheet("odd_endcap_sheet", endcap, {600, 600}, display::e_module_info);

    svg::file endcap_file;
    endcap_file._width = 1000;
    endcap_file.add_object(endcap_sheet);

    // Write out the file
    std::ofstream eout;
    eout.open("odd_endcap_sheet_module_info.svg");
    eout << endcap_file;
    eout.close();
}

TEST(display, endcap_sheet_grid_info) {

    endcap._name = "ODD Pixel Endcap (sample)";

    // Create the sheet
    svg::object endcap_sheet =
        display::endcap_sheet("odd_endcap_sheet", endcap, {600, 600}, display::e_grid_info);

    svg::file endcap_file;
    endcap_file._width = 1000;
    endcap_file.add_object(endcap_sheet);

    // Write out the file
    std::ofstream eout;
    eout.open("odd_endcap_sheet_grid_info.svg");
    eout << endcap_file;
    eout.close();
}
