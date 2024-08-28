// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>

#include "../common/playground.hpp"
#include "actsvg/core.hpp"

using namespace actsvg;

TEST(draw, markers) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Standard oriented makers: black
    style::marker m_default;
    auto m0 = draw::marker("m0", {50, 10}, m_default);

    style::marker m_s20;
    m_s20._size = 20;
    auto m1 = draw::marker("m1", {50, 50}, m_s20);

    style::marker m_arr{"<"};
    m_arr._size = 10;
    auto m2 = draw::marker("m2", {100, 80}, m_arr);

    style::marker m_darr{"<<"};
    m_darr._size = 10;
    auto m3 = draw::marker("m3", {100, 110}, m_darr);

    style::marker m_marr{"|<"};
    m_marr._size = 10;
    auto m4 = draw::marker("m4", {150, 170}, m_marr);

    style::marker m_mdarr{"|<<"};
    m_mdarr._size = 10;
    auto m5 = draw::marker("m5", {150, 210}, m_mdarr);

    // Flipped in red
    style::marker m_arr_r{"<"};
    m_arr_r._size = 10;
    m_arr_r._fill._fc = style::color{{255, 0, 0}};
    auto m2_r = draw::marker("m2_r", {100, -80}, m_arr_r, pi);

    style::marker m_darr_r{"<<"};
    m_darr_r._size = 10;
    m_darr_r._fill._fc = style::color{{255, 0, 0}};
    auto m3_r = draw::marker("m3_r", {100, -110}, m_darr_r, pi);

    style::marker m_marr_r{"|<"};
    m_marr_r._size = 10;
    m_marr_r._fill._fc = style::color{{255, 0, 0}};
    auto m4_r = draw::marker("m4_r", {150, -170}, m_marr_r, pi);

    style::marker m_mdarr_r{"|<<"};
    m_mdarr_r._size = 10;
    m_mdarr_r._fill._fc = style::color{{255, 0, 0}};
    auto m5_r = draw::marker("m5_r", {150, -210}, m_mdarr_r, pi);

    // 45 degree outwards pointing in blue
    style::marker m_arr_b{"<"};
    m_arr_b._size = 10;
    m_arr_b._fill._fc = style::color{{0, 0, 255}};
    auto m2_b = draw::marker("m2_r", {-100, 80}, m_arr_b, 0.75_scalar * pi);

    style::marker m_darr_b{"<<"};
    m_darr_b._size = 10;
    m_darr_b._fill._fc = style::color{{0, 0, 255}};
    auto m3_b = draw::marker("m3_r", {-100, 110}, m_darr_b, 0.75_scalar * pi);

    style::marker m_marr_b{"|<"};
    m_marr_b._size = 10;
    m_marr_b._fill._fc = style::color{{0, 0, 255}};
    auto m4_b = draw::marker("m4_r", {-150, 170}, m_marr_b, 0.75_scalar * pi);

    style::marker m_mdarr_b{"|<<"};
    m_mdarr_b._size = 10;
    m_mdarr_b._fill._fc = style::color{{0, 0, 255}};
    auto m5_b = draw::marker("m5_r", {-150, 210}, m_mdarr_b, 0.75_scalar * pi);

    // back-to-back measure test
    style::marker m_mdarr_g{"|<<"};
    m_mdarr_g._size = 10;
    m_mdarr_g._fill._fc = style::color{{0, 255, 0}};
    auto f_g = draw::marker("f_g", {-100, -100}, m_default);
    auto m5_f_g = draw::marker("m5_f_g", {-100, -100}, m_mdarr_g);
    auto m5_b_g = draw::marker("m5_f_g", {-100, -100}, m_mdarr_g, -pi);

    svg::file mfile;
    mfile.add_object(pg);
    // standard markers
    mfile.add_object(m0);
    mfile.add_object(m1);
    mfile.add_object(m2);
    mfile.add_object(m3);
    mfile.add_object(m4);
    mfile.add_object(m5);
    // flipped markers
    mfile.add_object(m2_r);
    mfile.add_object(m3_r);
    mfile.add_object(m4_r);
    mfile.add_object(m5_r);
    // outwards markers
    mfile.add_object(m2_b);
    mfile.add_object(m3_b);
    mfile.add_object(m4_b);
    mfile.add_object(m5_b);
    // back-to-back test
    mfile.add_object(f_g);
    mfile.add_object(m5_f_g);
    mfile.add_object(m5_b_g);

    std::ofstream tstream;
    tstream.open("test_core_markers.svg");
    tstream << mfile;
    tstream.close();
}
