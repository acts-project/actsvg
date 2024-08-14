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

    ASSERT_TRUE(v._name == "unnamed_volume");
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
    v._bound_values = {0., 40., 0., 400., pi, 0.};
    v._type = decltype(v)::type::e_cylinder;
    v._fill._fc._rgb = {0, 0, 0};
    v._fill._fc._opacity = 0.1_scalar;

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
    red._opacity = 0.1_scalar;
    std::vector<style::color> volumeColors = {red};
    v.colorize(volumeColors);

    svg::object v_red_zr = display::volume("cylinder_volume", v, views::z_r{});
    svg::file rfile_red_zr;
    rfile_red_zr.add_object(v_red_zr);

    rstream.open("test_meta_cylinder_volume_red_zr.svg");
    rstream << rfile_red_zr;
    rstream.close();

    // sectoral cylinder
    proto::volume<point3_container> v_sect;
    v_sect._index = 0u;

    // Negative endcap portal: nec
    proto::portal<point3_container> nec_sect;
    proto::surface<point3_container> s_nec_sect;
    s_nec_sect._type = proto::surface<point3_container>::type::e_disc;
    s_nec_sect._radii = {10., 110.};
    s_nec_sect._opening = {-0.25 * pi, 0.25 * pi};
    s_nec_sect._zparameters = {-400., 0.};
    s_nec_sect._fill = __nn_fill;

    // Assign the surface & link to self
    nec_sect._surface = s_nec_sect;
    link_to_self._start = {60., 0., -400.};
    link_to_self._end = {60., 0., -380.};
    nec_sect._volume_links = {link_to_self};

    // Positive endcap portal: pec
    proto::portal<point3_container> pec_sect;
    proto::surface<point3_container> s_pec_sect;
    s_pec_sect._type = proto::surface<point3_container>::type::e_disc;
    s_pec_sect._radii = {10., 110.};
    s_pec_sect._opening = {-0.25 * pi, 0.25 * pi};
    s_pec_sect._zparameters = {400., 0.};
    s_pec_sect._fill = __nn_fill;

    // Assign the surface & link to self
    pec_sect._surface = s_pec_sect;
    link_to_self._link_index = 0u;
    link_to_self._start = {60., 0., 400.};
    link_to_self._end = {60., 0., 380.};
    pec_sect._volume_links = {link_to_self};

    // Inner cylinder : ci_sect
    proto::portal<point3_container> ci_sect;
    proto::surface<point3_container> s_ci_sect;
    s_ci_sect._type = proto::surface<point3_container>::type::e_cylinder;
    s_ci_sect._radii = {10.};
    s_ci_sect._opening = {-0.25 * pi, 0.25 * pi};
    s_ci_sect._zparameters = {0., 400.};

    // Assign the surface & link to self
    ci_sect._surface = s_ci_sect;
    link_to_self._start = {10., 0., 0.};
    link_to_self._end = {30., 0., 0.};
    ci_sect._volume_links = {link_to_self};

    // Outer cylinder : co_sect
    proto::portal<point3_container> co_sect;
    proto::surface<point3_container> s_co_sect;
    s_co_sect._type = proto::surface<point3_container>::type::e_cylinder;
    s_co_sect._radii = {110.};
    s_co_sect._opening = {-0.25 * pi, 0.25 * pi};
    s_co_sect._zparameters = {0., 400.};

    // Assign the surface & link to self
    co_sect._surface = s_co_sect;
    link_to_self._start = {110., 0., 0.};
    link_to_self._end = {90., 0., 0.};
    co_sect._volume_links = {link_to_self};

    // Negative sector
    proto::portal<point3_container> neg_sect;
    proto::surface<point3_container> s_neg_sect;
    s_neg_sect._type = proto::surface<point3_container>::type::e_rectangle;

    scalar cos_n_phi = std::cos(s_co_sect._opening[0u]);
    scalar sin_n_phi = std::sin(s_co_sect._opening[0u]);
    s_neg_sect._vertices = {{10 * cos_n_phi, 10 * sin_n_phi, -400.},
                            {10 * cos_n_phi, 10 * sin_n_phi, 400},
                            {110 * cos_n_phi, 110 * sin_n_phi, 400},
                            {110 * cos_n_phi, 110 * sin_n_phi, -400}};

    // Assign the surface & link to self
    neg_sect._surface = s_neg_sect;
    link_to_self._start = {cos_n_phi * 60, sin_n_phi * 60, 0.};
    scalar d_neg_x = -20 * sin_n_phi;
    scalar d_neg_y = 20 * cos_n_phi;
    link_to_self._end = {60 * cos_n_phi + d_neg_x, 60 * sin_n_phi + d_neg_y,
                         0.};
    neg_sect._volume_links = {link_to_self};

    // Positive sector
    proto::portal<point3_container> pos_sect;
    proto::surface<point3_container> s_pos_sect;
    s_pos_sect._type = proto::surface<point3_container>::type::e_rectangle;

    scalar cos_p_phi = std::cos(s_co_sect._opening[1u]);
    scalar sin_p_phi = std::sin(s_co_sect._opening[1u]);
    s_pos_sect._vertices = {{10 * cos_p_phi, 10 * sin_p_phi, -400.},
                            {10 * cos_p_phi, 10 * sin_p_phi, 400},
                            {110 * cos_p_phi, 110 * sin_p_phi, 400},
                            {110 * cos_p_phi, 110 * sin_p_phi, -400}};

    // Assign the surface & link to self
    pos_sect._surface = s_pos_sect;
    link_to_self._start = {cos_p_phi * 60, sin_p_phi * 60, 0.};
    scalar d_pos_x = 20 * sin_p_phi;
    scalar d_pos_y = -20 * cos_n_phi;
    link_to_self._end = {60 * cos_p_phi + d_pos_x, 60 * sin_p_phi + d_pos_y,
                         0.};
    pos_sect._volume_links = {link_to_self};

    v_sect._portals = {nec_sect, pec_sect, ci_sect,
                       co_sect,  neg_sect, pos_sect};

    // Set up the volume parameters
    v_sect._bound_values = {10., 110., 0., 400., 0.25 * pi, 0.};
    v_sect._type = decltype(v_sect)::type::e_cylinder;
    v_sect._fill._fc._rgb = {0, 0, 255};
    v_sect._fill._fc._opacity = 0.1_scalar;

    // Test the volume in x-y view
    svg::object v_xy_sect =
        display::volume("cylinder_volume", v_sect, views::x_y{});

    svg::file rfile_xy_sect;
    rfile_xy_sect.add_object(v_xy_sect);

    rstream.open("test_meta_sectoral_cylinder_volume_xy.svg");
    rstream << rfile_xy_sect;
    rstream.close();
}
