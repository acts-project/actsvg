// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include "actsvg/core.hpp"
#include "actsvg/display/materials.hpp"
#include "actsvg/proto/grid.hpp"
#include "actsvg/proto/material.hpp"

using namespace actsvg;

TEST(proto, surface_material_x_y) {

    proto::grid grid;
    grid._type = proto::grid::type::e_x_y;
    grid._edges_0 = {-200, -100, 0, 100, 200};
    grid._edges_1 = {-100, 0, 100.};

    // Make a material matrix
    std::vector<std::vector<proto::material_slab> > m_matrix;
    for (std::size_t ib1 = 1u; ib1 <= grid._edges_1.size(); ++ib1) {
        std::vector<proto::material_slab> m_matrix_row;
        m_matrix_row.reserve(grid._edges_0.size());
        for (std::size_t ib0 = 1u; ib0 <= grid._edges_0.size(); ++ib0) {
            proto::material_slab s;
            // Create some material properties
            s._properties = {static_cast<scalar>(ib1 * 10 + ib0 * 1),
                             static_cast<scalar>(ib1 * 3 + ib0 * 1),
                             static_cast<scalar>(ib0 * 2),
                             static_cast<scalar>(ib0),
                             static_cast<scalar>(ib1 * 0.05),
                             1.};
            m_matrix_row.push_back(s);
        }
        m_matrix.push_back(m_matrix_row);
    }

    // Let's create an empty bin]
    m_matrix[1u][2u] = proto::material_slab{};

    // Create the surface material
    proto::surface_material m{m_matrix, grid};

    m._material_ranges = proto::material_ranges(m._material_matrix);

    style::gradient gradient;
    gradient._stops = {
        style::gradient::stop{0., style::color{style::rgb{0, 0, 255}}},
        style::gradient::stop{1., style::color{style::rgb{255, 0, 0}}}};
    m._gradient = gradient;

    auto sm_x_y = display::surface_material("x_y_material", m);

    svg::file file_sm_x_y;
    file_sm_x_y.add_object(sm_x_y);

    std::ofstream rstream;
    rstream.open("test_meta_surface_material_xy.svg");
    rstream << file_sm_x_y;
    rstream.close();
}

TEST(proto, surface_material_r_phi) {

    proto::grid grid;
    grid._type = proto::grid::type::e_r_phi;
    grid._edges_0 = {10, 25, 50, 75, 100, 125, 150};
    grid._edges_1 = {};
    for (unsigned int i = 0; i < 41; ++i) {
        grid._edges_1.push_back(-pi + i * pi / 20._scalar);
    }

    // Make a material matrix - with some sparice rows
    std::vector<std::vector<proto::material_slab> > m_matrix;
    for (std::size_t ib1 = 1u; ib1 <= grid._edges_1.size(); ++ib1) {
        std::vector<proto::material_slab> m_matrix_row;
        m_matrix_row.reserve(grid._edges_0.size());
        for (std::size_t ib0 = 1u; ib0 <= grid._edges_0.size(); ++ib0) {
            proto::material_slab s;
            // Create some material properties
            if ((ib0 == 2 || ib0 == 4) && (ib1 % 2)) {
                s._properties = {0., 0., 0., 0., 0., 0.};
            } else {
                s._properties = {static_cast<scalar>(ib1 * 1 + ib0 * 10),
                                 static_cast<scalar>(ib1 * 1 + ib0 * 3),
                                 static_cast<scalar>(ib0 * 2),
                                 static_cast<scalar>(ib0),
                                 static_cast<scalar>(ib1 * 0.05),
                                 1.};
            }
            m_matrix_row.push_back(s);
        }
        m_matrix.push_back(m_matrix_row);
    }

    // Create the surface material
    proto::surface_material m{m_matrix, grid};
    m._material_ranges = proto::material_ranges(m._material_matrix);

    style::gradient gradient;
    gradient._stops = {
        style::gradient::stop{0., style::color{style::rgb{0, 0, 255}}},
        style::gradient::stop{0.05, style::color{style::rgb{0, 255, 0}}},
        style::gradient::stop{0.1, style::color{style::rgb{255, 255, 0}}},
        style::gradient::stop{1., style::color{style::rgb{255, 0, 0}}}};
    m._gradient = gradient;
    m._info_pos = {-200, 200};

    m._gradient_pos = {200, -150};
    m._gradient_box = {20, 300};
    m._gradient_font._size = 10;

    m._gradient_label = style::label{"t/X0", style::label::horizontal::right,
                                     style::label::vertical::top};

    auto sm_r_phi = display::surface_material("r_phi_material", m);

    svg::file file_sm_r_phi;
    file_sm_r_phi.add_object(sm_r_phi);

    std::ofstream rstream;
    rstream.open("test_meta_surface_material_r_phi.svg");
    rstream << file_sm_r_phi;
    rstream.close();
}
