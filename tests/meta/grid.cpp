// This file is part of the actsvg packge.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/proto/grid.hpp"

#include <gtest/gtest.h>

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

    svg::object d_grid_r_phi = display::grid("grid_r_phi", g);

    svg::file gfile_r_phi;
    gfile_r_phi.add_object(d_grid_r_phi);

    std::ofstream rstream;
    rstream.open("test_meta_grid_r_phi_full.svg");
    rstream << gfile_r_phi;
    rstream.close();
}
