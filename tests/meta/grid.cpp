// This file is part of the actsvg packge.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/proto/grid.hpp"

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "actsvg/display/grids.hpp"
#include "actsvg/styles/defaults.hpp"

using namespace actsvg;

TEST(proto, grid_x_y) {

    proto::grid g;
    g._type = proto::grid::type::e_x_y;
    ASSERT_TRUE(g._type == proto::grid::type::e_x_y);

    g._edges_0 = {-200, -100, 0, 100, 200};
    g._edges_1 = {-500, -100, 0, 200.};

    svg::object d_grid_x_y = display::grid("grid_x_y", g);

    svg::file gfile_x_y;
    gfile_x_y.add_object(d_grid_x_y);

    std::ofstream rstream;
    rstream.open("test_meta_grid_x_y.svg");
    rstream << gfile_x_y;
    rstream.close();
}

TEST(proto, grid_z_phi_closed) {

    proto::grid g;
    g._type = proto::grid::type::e_z_phi;
    g._reference_r = 100.;
    ASSERT_TRUE(g._type == proto::grid::type::e_z_phi);

    g._edges_0 = {-200, -100, 0, 100, 200};
    g._edges_1 = {-0.75 * M_PI, -0.5 * M_PI, -0.25 * M_PI,
                  0.,           0.25 * M_PI, 0.5 * M_PI};
    // Create connection text
    size_t ig = 0u;
    std::vector<svg::object> tile_ids;
    for (size_t i1 = 0; i1 + 1u < g._edges_1.size(); ++i1) {
        for (size_t i0 = 0; i0 + 1u < g._edges_0.size(); ++i0) {
            std::string global_id = "g = [";
            global_id += std::to_string(ig++) + std::string("]");
            std::string local_id = "l = [";
            local_id += std::to_string(i0) + std::string(", ");
            local_id += std::to_string(i1) + std::string("]");
            std::string tile_name = std::string("tile_") + std::to_string(i0) +
                                    std::string("_") + std::to_string(i1);
            point2 text_pos = {static_cast<scalar>(0.9 * g._edges_0[i0] +
                                                   0.1 * g._edges_0[i0 + 1u]),
                               static_cast<scalar>((0.3 * g._edges_1[i1] +
                                                    0.7 * g._edges_1[i1 + 1u]) *
                                                   g._reference_r)};
            tile_ids.push_back(
                draw::text(tile_name, text_pos, {global_id, local_id}));
        }
    }

    g._connections = tile_ids;

    svg::object d_grid_z_phi = display::grid("grid_z_phi", g);

    svg::file gfile_z_phi;
    gfile_z_phi.add_object(d_grid_z_phi);

    std::ofstream rstream;
    rstream.open("test_meta_grid_z_phi.svg");
    rstream << gfile_z_phi;
    rstream.close();
}

TEST(proto, grid_r_phi_open) {

    proto::grid g;
    ASSERT_TRUE(g._type == proto::grid::type::e_r_phi);

    g._edges_0 = {200., 300., 400., 500.};
    g._edges_1 = {-0.75 * M_PI, -0.5 * M_PI, -0.25 * M_PI,
                  0.,           0.25 * M_PI, 0.5 * M_PI};

    svg::object d_grid_r_phi = display::grid("grid_r_phi", g);

    svg::file gfile_r_phi;
    gfile_r_phi.add_object(d_grid_r_phi);

    std::ofstream rstream;
    rstream.open("test_meta_grid_r_phi_open.svg");
    rstream << gfile_r_phi;
    rstream.close();
}

TEST(proto, grid_r_phi_full) {

    proto::grid g;
    ASSERT_TRUE(g._type == proto::grid::type::e_r_phi);

    g._edges_0 = {200., 300., 400., 500.};
    g._edges_1 = {-M_PI,       -0.75 * M_PI, -0.5 * M_PI, -0.25 * M_PI, 0.,
                  0.25 * M_PI, 0.5 * M_PI,   0.75 * M_PI, M_PI};

    // Create connection text
    size_t ig = 0u;
    std::vector<svg::object> tile_ids;
    for (size_t i1 = 0; i1 + 1u < g._edges_1.size(); ++i1) {
        for (size_t i0 = 0; i0 + 1u < g._edges_0.size(); ++i0) {
            scalar r = 0.5_scalar * (g._edges_0[i0] + g._edges_0[i0 + 1u]);
            scalar phi = 0.5_scalar * (g._edges_1[i1] + g._edges_1[i1 + 1u]);
            std::string global_id = "g = [";
            global_id += std::to_string(ig++) + std::string("]");
            std::string local_id = "l = [";
            local_id += std::to_string(i0) + std::string(", ");
            local_id += std::to_string(i1) + std::string("]");
            std::string tile_name = std::string("tile_") + std::to_string(i0) +
                                    std::string("_") + std::to_string(i1);

            point2 text_pos = {r * std::cos(phi), r * std::sin(phi)};
            tile_ids.push_back(
                draw::text(tile_name, text_pos, {global_id, local_id}));
        }
    }

    g._connections = tile_ids;

    svg::object d_grid_r_phi = display::grid("grid_r_phi", g);

    svg::file gfile_r_phi;
    gfile_r_phi.add_object(d_grid_r_phi);

    std::ofstream rstream;
    rstream.open("test_meta_grid_r_phi_full.svg");
    rstream << gfile_r_phi;
    rstream.close();
}

TEST(proto, grid_r_phi_full_connections) {

    proto::grid g;
    ASSERT_TRUE(g._type == proto::grid::type::e_r_phi);

    g._edges_0 = {200., 300., 400., 500.};
    g._edges_1 = {-M_PI,       -0.75 * M_PI, -0.5 * M_PI, -0.25 * M_PI, 0.,
                  0.25 * M_PI, 0.5 * M_PI,   0.75 * M_PI, M_PI};

    // Create connection text
    std::vector<svg::object> rows;
    svg::object row0 = draw::text("row_0", {0., 0.}, {"row 0"});
    svg::object row1 = draw::text("row_1", {0., 0.}, {"row 1"});
    svg::object row2 = draw::text("row_2", {0., 0.}, {"row 2"});
    g._connections = {row0, row1, row2};

    std::vector<std::vector<size_t> > associations;
    for (size_t i0 = 0; i0 + 1u < g._edges_0.size(); ++i0) {
        for (size_t i1 = 0; i1 + 1u < g._edges_1.size(); ++i1) {
            associations.push_back({i0});
        }
    }

    g._connection_types = {connectors::type::e_highlight};
    g._connection_associations = associations;

    svg::object d_grid_r_phi = display::grid("grid_r_phi", g);

    svg::file gfile_r_phi;
    gfile_r_phi.add_object(d_grid_r_phi);

    std::ofstream rstream;
    rstream.open("test_meta_grid_r_phi_connections.svg");
    rstream << gfile_r_phi;
    rstream.close();
}
