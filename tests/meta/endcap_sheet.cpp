// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <exception>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/data/odd_pixel_ec.hpp"
#include "actsvg/meta.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

namespace {

proto::volume<point3_container> generate_endcap(
    const proto::surface<point3_container>::type& m_t_) noexcept(false) {

    unsigned int sectors = 16;
    scalar half_opening = static_cast<scalar>(M_PI / sectors);

    proto::surface<point3_container> inner_template_surface;
    inner_template_surface._name = "inner_surface_";

    proto::surface<point3_container> middle_template_surface;
    middle_template_surface._name = "middle_surface_";

    proto::surface<point3_container> outer_template_surface;
    outer_template_surface._name = "outer_surface_";

    std::vector<scalar> rad_pos = {};

    // Disc case
    if (m_t_ == proto::surface<point3_container>::e_disc) {
        inner_template_surface._type = m_t_;
        inner_template_surface._radii = {100, 200};
        inner_template_surface._opening = {
            static_cast<scalar>(-half_opening - 0.05),
            static_cast<scalar>(half_opening + 0.05)};
        inner_template_surface._measures = {
            100, 200, static_cast<scalar>(half_opening + 0.05), 0.};

        middle_template_surface._type =
            proto::surface<point3_container>::e_disc;
        middle_template_surface._radii = {190, 320};
        middle_template_surface._opening = {
            static_cast<scalar>(-half_opening - 0.05),
            static_cast<scalar>(half_opening + 0.05)};
        middle_template_surface._measures = {
            190, 320, static_cast<scalar>(half_opening + 0.05), 0.};

        outer_template_surface._type = proto::surface<point3_container>::e_disc;
        outer_template_surface._radii = {310, 450};
        outer_template_surface._opening = {
            static_cast<scalar>(-half_opening - 0.05),
            static_cast<scalar>(half_opening + 0.05)};
        outer_template_surface._measures = {
            310, 450, static_cast<scalar>(half_opening + 0.05), 0.};

    } else if (m_t_ == proto::surface<point3_container>::e_trapez) {
        inner_template_surface._type = m_t_;
        inner_template_surface._measures = {26., 40., 55.};
        inner_template_surface._vertices = {
            {-inner_template_surface._measures[0],
             -inner_template_surface._measures[2], 0},
            {inner_template_surface._measures[0],
             -inner_template_surface._measures[2], 0},
            {inner_template_surface._measures[1],
             inner_template_surface._measures[2], 0},
            {-inner_template_surface._measures[1],
             inner_template_surface._measures[2], 0}};

        middle_template_surface._type = m_t_;
        middle_template_surface._measures = {40., 68., 75.};
        middle_template_surface._vertices = {
            {-middle_template_surface._measures[0],
             -middle_template_surface._measures[2], 0},
            {middle_template_surface._measures[0],
             -middle_template_surface._measures[2], 0},
            {middle_template_surface._measures[1],
             middle_template_surface._measures[2], 0},
            {-middle_template_surface._measures[1],
             middle_template_surface._measures[2], 0}};

        outer_template_surface._type = m_t_;
        outer_template_surface._measures = {68., 100., 75.};
        outer_template_surface._vertices = {
            {-outer_template_surface._measures[0],
             -outer_template_surface._measures[2], 0},
            {outer_template_surface._measures[0],
             -outer_template_surface._measures[2], 0},
            {outer_template_surface._measures[1],
             outer_template_surface._measures[2], 0},
            {-outer_template_surface._measures[1],
             outer_template_surface._measures[2], 0}};

        rad_pos = {150., 265., 400.};
    }
    // Template surfaces
    std::vector<proto::surface<point3_container>> template_surfaces;
    template_surfaces.push_back(inner_template_surface);
    template_surfaces.push_back(middle_template_surface);
    template_surfaces.push_back(outer_template_surface);

    // Actual surfaces
    std::vector<proto::surface<point3_container>> surfaces;
    views::x_y x_y_view;

    // Create the endcap volume
    proto::volume<point3_container> endcap;
    endcap._name = "Endcap with templates";

    // Loop over the templates and place them
    for (const auto [is, ts] : utils::enumerate(template_surfaces)) {

        // Create the template object
        auto template_object =
            display::surface(ts._name + "_template", ts, x_y_view);
        template_object._sterile = true;
        // Get the radiii
        scalar ri = ts._radii[0];
        scalar ro = ts._radii[1];

        for (unsigned int isc = 0; isc < sectors; ++isc) {

            // Create the surface from a template
            proto::surface<point3_container> s =
                proto::surface<point3_container>::from_template(
                    ts, template_object, ts._name + std::to_string(isc));

            // Phi
            scalar phi = isc * 2 * half_opening;
            phi -= phi > M_PI ? 2 * M_PI : 0.;

            /// Add some descriptive text
            s._info = {
                "module_" + std::to_string(is) + "_" + std::to_string(isc),
                "phi =" + utils::to_string(phi)};

            // The transform
            style::transform t;
            if (m_t_ == proto::surface<point3_container>::e_trapez) {
                t._tr[0] = static_cast<scalar>(rad_pos[is] * std::cos(phi));
                t._tr[1] = static_cast<scalar>(rad_pos[is] * std::sin(phi));
                s._info.push_back(std::string("center = ") +
                                  utils::to_string(t._tr));
                t._rot[0] = static_cast<scalar>(phi * 180 / M_PI - 90.);
            } else {
                t._rot[0] = static_cast<scalar>(phi * 180 / M_PI);
            }

            s._transform = t;

            if (m_t_ == proto::surface<point3_container>::e_disc) {
                // Let's fill some vertices
                scalar phi_low = phi - half_opening;
                scalar phi_high = phi + half_opening;
                scalar cos_phi_low = std::cos(phi_low);
                scalar sin_phi_low = std::sin(phi_low);
                scalar cos_phi_high = std::cos(phi_high);
                scalar sin_phi_high = std::cos(phi_high);
                point3 A = {ri * cos_phi_low, ri * sin_phi_low, 0.};
                point3 B = {ri * std::cos(phi), ri * std::sin(phi), 0.};
                point3 C = {ri * cos_phi_high, ri * sin_phi_high, 0.};
                point3 D = {ro * cos_phi_high, ro * sin_phi_high, 0.};
                point3 E = {ro * std::cos(phi), ro * std::sin(phi), 0.};
                point3 F = {ro * cos_phi_low, ro * sin_phi_low, 0.};
                s._vertices = {A, B, C, D, E, F};
            }

            // Set the template
            s._template_object = template_object;
            surfaces.push_back(s);
        }
    }

    // Add the surfaces
    endcap._surfaces = surfaces;
    return endcap;
}

}  // namespace

// Helper method to generate the endcap volume description
template <typename point3_container_type = point3_container>
proto::volume<point3_container_type> generate_endcap_volume() {

    proto::volume<point3_container> endcap;

    // ODD module templates
    proto::surface<point3_container> endcap_inner_template;
    endcap_inner_template._vertices = {
        {-8.5, -34, 0.}, {8.5, -34, 0.}, {14.5, 34., 0.}, {-14.5, 34., 0.}};
    endcap_inner_template._measures = {8.5, 14.5, 34.};

    proto::surface<point3_container> endcap_outer_template;
    endcap_outer_template._vertices = {
        {-10.5, -34., 0.}, {10.5, -34., 0.}, {16.5, 34., 0.}, {-16.5, 34., 0}};
    endcap_outer_template._measures = {10.5, 16.5, 34.};

    size_t number_of_modules = data::odd_pixel_ec.size() / 4u;
    endcap._surfaces.reserve(number_of_modules);
    for (size_t im = 0; im < number_of_modules; ++im) {
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

TEST(display, endcap_sheet_grid_info) {

    endcap._name = "ODD Pixel Endcap (sample)";

    // Create the sheet
    svg::object endcap_sheet = display::endcap_sheet(
        "sheet_odd_endcap", endcap, {600, 600}, display::e_grid_info);

    svg::file endcap_file;
    endcap_file._width = 1000;
    endcap_file.add_object(endcap_sheet);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_odd_endcap_grid_info.svg");
    eout << endcap_file;
    eout.close();
}

TEST(display, endcap_sheet_sec_module_info_ref) {

    auto sector_endcap =
        generate_endcap(proto::surface<point3_container>::e_disc);

    // Create the sheet
    svg::object endcap_sheet_fs =
        display::endcap_sheet("sector_endcap_sheet", sector_endcap, {600, 600},
                              display::e_module_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_sector_module_info.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_trap_module_info_ref) {

    auto trapez_endcap =
        generate_endcap(proto::surface<point3_container>::e_trapez);

    // Create the sheet
    svg::object endcap_sheet_fs =
        display::endcap_sheet("trapez_endcap_sheet", trapez_endcap, {600, 600},
                              display::e_module_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_trapez_module_info.svg");
    eout << endcap_file_fs;
    eout.close();
}
