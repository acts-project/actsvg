// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

TEST(proto, cylindrical_volume) {

    // Create and define a volume
    proto::volume<point3_container> v;
    v._index = 0u;

    ASSERT_TRUE(v._name == "unnamed");
    ASSERT_TRUE(v._surfaces.empty());
    ASSERT_TRUE(v._portals.empty());

    // Negative endcap portal: nec
    proto::portal<point3_container> nec;
    proto::surface<point3_container> s_nec;
    s_nec._type = proto::surface<point3_container>::type::e_disc;
    s_nec._radii = {0., 40.};
    s_nec._zparameters = {-400., 0.};

    // Assign the surface & link to self
    nec._surface = s_nec;
    proto::portal<point3_container>::link link_to_self;
    link_to_self._start = {20., 0., -400.};
    link_to_self._end = {20., 0., -380.};
    nec._volume_links = {link_to_self};

    // Positive endcap portal: pec
    proto::portal<point3_container> pec;
    proto::surface<point3_container> s_pec;
    s_pec._type = proto::surface<point3_container>::type::e_disc;
    s_pec._radii = {0., 40.};
    s_pec._zparameters = {400., 0.};

    // Assign the surface & link to self
    pec._surface = s_pec;
    link_to_self._link_index = 0u;
    link_to_self._start = {20., 0., 400.};
    link_to_self._end = {20., 0., 380.};
    pec._volume_links = {link_to_self};

    // Cover cylinder : c
    proto::portal<point3_container> c;
    proto::surface<point3_container> s_c;
    s_c._type = proto::surface<point3_container>::type::e_cylinder;
    s_c._radii = {0., 40.};
    s_c._zparameters = {0., 400.};

    // Assign the surface & link to self
    c._surface = s_c;
    link_to_self._start = {40., 0., 0.};
    link_to_self._end = {20., 0., 0.};
    c._volume_links = {link_to_self};

    v._portals = {nec, c, pec};

    // Set up the volume parameters
    v._bound_values = {0., 40., 0., 400., M_PI, 0.};
    v._type = decltype(v)::type::e_cylinder;
    v._fill._fc._rgb = {0, 0, 0};
    v._fill._fc._opacity = 0.1;

    // Test the volume in x-y view
    svg::object v_xy = display::volume("cylinder_volume", v, views::x_y{});

    svg::file rfile_xy;
    rfile_xy.add_object(v_xy);

    std::ofstream rstream;
    rstream.open("test_meta_cylinder_volume_xy.svg");
    rstream << rfile_xy;
    rstream.close();

    // Test the disc in z-r view
    svg::object v_zr = display::volume("cylinder_volume", v, views::z_r{});

    svg::file rfile_zr;
    rfile_zr.add_object(v_zr);

    rstream.open("test_meta_cylinder_volume_zr.svg");
    rstream << rfile_zr;
    rstream.close();

    style::color red({{255, 0, 0}});
    red._opacity = 0.1;
    std::vector<style::color> volumeColors = {red};
    v.colorize(volumeColors);

    svg::object v_red_zr = display::volume("cylinder_volume", v, views::z_r{});
    svg::file rfile_red_zr;
    rfile_red_zr.add_object(v_red_zr);

    rstream.open("test_meta_cylinder_volume_red_zr.svg");
    rstream << rfile_red_zr;
    rstream.close();
}
