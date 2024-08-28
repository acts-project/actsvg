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
#include "actsvg/data/odd_pixel_endcap.hpp"
#include "actsvg/meta.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

namespace {

/** Generate an endcap with either disc or trapezoid volumes
 *
 * @param m_t_ is the surface type
 * @param d_s_ is a boolean steering if a double-sided endcap is created
 * @param p_s_ is a boolean steering if a passive surface is built
 * @param f_t_ is a boolean indicating if a template module is used (or not)
 *
 * @return a proto volume for the endcap
 *
 **/
proto::volume<point3_container> generate_endcap(
    const proto::surface<point3_container>::type& m_t_, bool d_s_ = false,
    bool p_s_ = false, bool f_t_ = true) noexcept(false) {

    unsigned int sectors = 16;
    scalar half_opening = static_cast<scalar>(M_PI / sectors);

    // Create the endcap volume
    proto::volume<point3_container> endcap;
    endcap._name = f_t_ ? "Endcap with templates" : "Endcap with individuals";

    // The proto surfaces for the templates
    proto::surface<point3_container> inner_template_surface;
    inner_template_surface._name = "inner_surface_";

    proto::surface<point3_container> middle_template_surface;
    middle_template_surface._name = "middle_surface_";

    proto::surface<point3_container> outer_template_surface;
    outer_template_surface._name = "outer_surface_";

    std::vector<scalar> rad_pos = {};

    endcap._surface_grid._type = proto::grid::e_r_phi;
    endcap._surface_grid._edges_0 = {};
    endcap._surface_grid._edges_1 = {};

    // Surface batches
    std::vector<proto::surface<point3_container>> regular;
    std::vector<proto::surface<point3_container>> support;
    std::vector<proto::surface<point3_container>> backside;

    // Disc case
    if (m_t_ == proto::surface<point3_container>::type::e_disc) {

        inner_template_surface._type = m_t_;
        inner_template_surface._sf_type =
            proto::surface<point3_container>::sf_type::e_sensitive;
        inner_template_surface._radii = {100, 200};
        inner_template_surface._opening = {
            static_cast<scalar>(-half_opening - 0.05),
            static_cast<scalar>(half_opening + 0.05)};
        inner_template_surface._measures = {
            100, 200, static_cast<scalar>(half_opening + 0.05), 0.};

        middle_template_surface._type =
            proto::surface<point3_container>::type::e_disc;
        middle_template_surface._sf_type =
            proto::surface<point3_container>::sf_type::e_sensitive;
        middle_template_surface._radii = {190, 320};
        middle_template_surface._opening = {
            static_cast<scalar>(-half_opening - 0.05),
            static_cast<scalar>(half_opening + 0.05)};
        middle_template_surface._measures = {
            190, 320, static_cast<scalar>(half_opening + 0.05), 0.};

        outer_template_surface._type =
            proto::surface<point3_container>::type::e_disc;
        outer_template_surface._sf_type =
            proto::surface<point3_container>::sf_type::e_sensitive;
        outer_template_surface._radii = {310, 450};
        outer_template_surface._opening = {
            static_cast<scalar>(-half_opening - 0.05),
            static_cast<scalar>(half_opening + 0.05)};
        outer_template_surface._measures = {
            310, 450, static_cast<scalar>(half_opening + 0.05), 0.};

        endcap._surface_grid._edges_0 = {100., 195., 305., 450.};

    } else if (m_t_ == proto::surface<point3_container>::type::e_trapez) {
        inner_template_surface._type = m_t_;
        inner_template_surface._sf_type =
            proto::surface<point3_container>::sf_type::e_sensitive;
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
        middle_template_surface._sf_type =
            proto::surface<point3_container>::sf_type::e_sensitive;
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
        outer_template_surface._sf_type =
            proto::surface<point3_container>::sf_type::e_sensitive;
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
        endcap._surface_grid._edges_0 = {95., 200., 310., 480.};
    }
    // Template surfaces
    std::vector<proto::surface<point3_container>> template_surfaces;
    template_surfaces.push_back(inner_template_surface);
    template_surfaces.push_back(middle_template_surface);
    template_surfaces.push_back(outer_template_surface);

    // A view for the template
    views::x_y x_y_view;

    std::vector<std::vector<int>> assoc_rows = {{0, 1}, {-1, 0, 1}, {-1, 0}};

    std::vector<std::vector<size_t>> bin_associations;
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

            // Current phi-value
            scalar phi = -pi + isc * 2_scalar * half_opening;

            // Create the surface
            proto::surface<point3_container> s;
            // (a) from template
            if (f_t_) {
                s = proto::surface<point3_container>::from_template(
                    ts, template_object, ts._name + "_" + std::to_string(isc));
            } else {
                // (b) as individual objects
                s._name = ts._name + std::to_string(isc);
                s._type = ts._type;
                s._sf_type = ts._sf_type;
                s._aux_info = ts._aux_info;
                s._vertices = ts._vertices;
                s._measures = ts._measures;
                s._radii = ts._radii;
                s._opening = ts._opening;
                s._fill = ts._fill;
                s._stroke = ts._stroke;
            }

            // Fill the grid phi-values
            if (is == 0) {
                if (isc == 0) {
                    endcap._surface_grid._edges_1.push_back(phi - half_opening);
                }
                endcap._surface_grid._edges_1.push_back(phi + half_opening);
            }

            // Fill the grid associations
            std::vector<size_t> bin_assoc;
            size_t current_bin = isc + is * sectors;
            size_t next_lower =
                current_bin > 0 ? current_bin - 1u : sectors - 1u;
            size_t next_higher =
                current_bin + 1u == sectors ? 0 : current_bin + 1u;
            std::vector<size_t> neighbors = {next_lower, current_bin,
                                             next_higher};

            for (auto ioff : assoc_rows[is]) {
                for (auto n : neighbors) {
                    int assoc = ioff * int(sectors) + n;
                    bin_assoc.push_back(assoc);
                }
            }
            bin_associations.push_back(bin_assoc);

            /// Add some descriptive text
            s._aux_info["module_info"] = {
                "module " + std::to_string(is) + " " + std::to_string(isc),
                "phi =" + utils::to_string(phi)};

            s._aux_info["grid_info"] = {"* module " + std::to_string(is) + " " +
                                        std::to_string(isc)};

            // The transform
            style::transform t;
            if (m_t_ == proto::surface<point3_container>::type::e_trapez) {
                t._tr[0] = static_cast<scalar>(rad_pos[is] * std::cos(phi));
                t._tr[1] = static_cast<scalar>(rad_pos[is] * std::sin(phi));
                s._aux_info["module_info"] = {std::string("center = ") +
                                              utils::to_string(t._tr)};
                t._rot[0] = static_cast<scalar>(phi * 180 / M_PI - 90.);
            } else {
                t._rot[0] = static_cast<scalar>(phi * 180 / M_PI);
            }

            s._transform = t;

            if (m_t_ == proto::surface<point3_container>::type::e_disc) {
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

            if (d_s_) {
                s._aux_info["side"] = {"+"};
                // make a copy and call flag it as negative side
                auto so = s;
                so._name = s._name + "_backside";
                // Aux info
                so._aux_info["side"] = {"-"};
                so._aux_info["grid_info"] = {"* bs module " +
                                             std::to_string(is) + " " +
                                             std::to_string(isc)};
                // Measures
                so._measures = s._measures;
                so._radii = s._radii;
                so._opening = s._opening;
                // Change fill
                so._fill = __bs_fill;
                // Set the template
                backside.push_back(so);
            }
            // Set the template
            regular.push_back(s);
        }
    }

    // Store the surface vectors
    endcap._surfaces.push_back(regular);
    endcap._grid_associations.push_back(bin_associations);
    // The passive surface if existing
    if (p_s_) {
        // The proto surfaces for the templates
        proto::surface<point3_container> support_disk;
        support_disk._name = "support_disk";
        support_disk._type = proto::surface<point3_container>::type::e_disc;
        support_disk._sf_type =
            proto::surface<point3_container>::sf_type::e_passive;
        scalar r_min = 100.;
        scalar r_max = 450.;
        point3 A = {r_max, 0., 0.};
        point3 B = {0., r_max, 0.};
        point3 C = {-r_max, 0., 0.};
        point3 D = {0., -r_max, 0.};
        support_disk._vertices = {A, B, C, D};
        support_disk._radii = {r_min, r_max};
        support_disk._measures = {r_min, r_max, static_cast<scalar>(M_PI), 0.};
        // Style
        support_disk._fill = __ss_fill;
        // Auxiliary info
        support_disk._aux_info["module_info"] = {"support disk"};
        support_disk._aux_info["grid_info"] = {"* suport disk"};
        // Store them
        support.push_back(support_disk);
        endcap._surfaces.push_back(support);
        // Create single bin association for the support disk
        std::vector<std::vector<size_t>> single_association;
        for (unsigned int its = 0; its < template_surfaces.size(); ++its) {
            for (unsigned int isc = 0; isc < sectors; ++isc) {
                single_association.push_back({0u});
            }
        }
        endcap._grid_associations.push_back(single_association);
    }
    // The backside surfaces if existing
    if (not backside.empty()) {
        endcap._surfaces.push_back(backside);
        endcap._grid_associations.push_back(bin_associations);
    }

    return endcap;
}

}  // namespace

TEST(display, endcap_sheet_sec_module_info_tmp) {

    auto sector_endcap =
        generate_endcap(proto::surface<point3_container>::type::e_disc);

    // Create the sheet
    svg::object endcap_sheet_fs =
        display::endcap_sheet("sector_endcap_sheet", sector_endcap, {600, 600},
                              display::e_module_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_t_sector_module_info.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_sec_module_info_ind) {

    auto sector_endcap = generate_endcap(
        proto::surface<point3_container>::type::e_disc, false, false);

    // Create the sheet
    svg::object endcap_sheet_fs =
        display::endcap_sheet("sector_endcap_sheet", sector_endcap, {600, 600},
                              display::e_module_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_i_sector_module_info.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_sec_module_grid_tmp) {

    auto sector_endcap =
        generate_endcap(proto::surface<point3_container>::type::e_disc);

    // Create the sheet
    svg::object endcap_sheet_fs = display::endcap_sheet(
        "sector_endcap_sheet", sector_endcap, {600, 600}, display::e_grid_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_t_sector_grid_info.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_sec_module_grid_ind) {

    auto sector_endcap = generate_endcap(
        proto::surface<point3_container>::type::e_disc, false, false, false);

    // Create the sheet
    svg::object endcap_sheet_fs = display::endcap_sheet(
        "sector_endcap_sheet", sector_endcap, {600, 600}, display::e_grid_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_i_sector_grid_info.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_trap_module_info_tmp) {

    auto trapez_endcap =
        generate_endcap(proto::surface<point3_container>::type::e_trapez);

    // Create the sheet
    svg::object endcap_sheet_fs =
        display::endcap_sheet("trapez_endcap_sheet", trapez_endcap, {600, 600},
                              display::e_module_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_t_trapez_module_info.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_trap_module_info_ind) {

    auto trapez_endcap = generate_endcap(
        proto::surface<point3_container>::type::e_trapez, false, false, false);

    // Create the sheet
    svg::object endcap_sheet_fs =
        display::endcap_sheet("trapez_endcap_sheet", trapez_endcap, {600, 600},
                              display::e_module_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_i_trapez_module_info.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_trap_grid_info_tmp) {

    auto trapez_endcap =
        generate_endcap(proto::surface<point3_container>::type::e_trapez);

    // Create the sheet
    svg::object endcap_sheet_fs = display::endcap_sheet(
        "trapez_endcap_sheet", trapez_endcap, {600, 600}, display::e_grid_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_t_trapez_grid_info.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_trap_grid_info_ind) {

    auto trapez_endcap = generate_endcap(
        proto::surface<point3_container>::type::e_trapez, false, false, false);

    // Create the sheet
    svg::object endcap_sheet_fs = display::endcap_sheet(
        "trapez_endcap_sheet", trapez_endcap, {600, 600}, display::e_grid_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_i_trapez_grid_info.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_trap_grid_info_tmp_backside) {

    auto trapez_endcap =
        generate_endcap(proto::surface<point3_container>::type::e_trapez, true);

    trapez_endcap._name = "Endcap with backside modules";

    // Create the sheet
    svg::object endcap_sheet_fs =
        display::endcap_sheet("trapez_endcap_sheet_w_backside", trapez_endcap,
                              {600, 600}, display::e_grid_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_t_trapez_grid_info_backside.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_trap_grid_info_tmp_backside_passive) {

    auto trapez_endcap = generate_endcap(
        proto::surface<point3_container>::type::e_trapez, true, true);

    trapez_endcap._name = "Endcap with backside modules and passive plate";

    // Create the sheet
    svg::object endcap_sheet_fs =
        display::endcap_sheet("trapez_endcap_sheet_w_backside_passive",
                              trapez_endcap, {600, 600}, display::e_grid_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_t_trapez_grid_info_backside_passive.svg");
    eout << endcap_file_fs;
    eout.close();
}

TEST(display, endcap_sheet_trap_module_info_tmp_backside_passive) {

    auto trapez_endcap = generate_endcap(
        proto::surface<point3_container>::type::e_trapez, true, true);

    trapez_endcap._name = "Endcap with backside modules and passive plate";

    // Create the sheet
    svg::object endcap_sheet_fs = display::endcap_sheet(
        "trapez_endcap_sheet_w_backside_passive", trapez_endcap, {600, 600},
        display::e_module_info);

    svg::file endcap_file_fs;
    endcap_file_fs._width = 1000;
    endcap_file_fs.add_object(endcap_sheet_fs);

    // Write out the file
    std::ofstream eout;
    eout.open("sheet_endcap_t_trapez_module_info_backside_passive.svg");
    eout << endcap_file_fs;
    eout.close();
}
