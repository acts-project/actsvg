// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/proto/detector.hpp"

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

TEST(proto, cylindrical_detector) {

    // Create and define a volume
    proto::detector<point3_container> d;

    ASSERT_TRUE(d._name == "unnamed_detector");
    ASSERT_TRUE(d._volumes.empty());

    d._name = "bp_pix";

    // beam pipe
    proto::volume<point3_container> bp;
    bp._index = 0u;
    // Set up the volume parameters
    bp._bound_values = {0., 40., 0., 400., M_PI, 0.};
    bp._type = decltype(bp)::type::e_cylinder;

    proto::volume<point3_container> pix_nec;
    pix_nec._index = 1u;
    // Set up the volume parameters
    pix_nec._bound_values = {40., 100., -350., 50., M_PI, 0.};
    pix_nec._type = decltype(bp)::type::e_cylinder;

    proto::volume<point3_container> pix_b;
    pix_b._index = 2u;
    // Set up the volume parameters
    pix_b._bound_values = {40., 100., 0., 300., M_PI, 0.};
    pix_b._type = decltype(bp)::type::e_cylinder;

    proto::volume<point3_container> pix_pec;
    pix_pec._index = 3u;
    // Set up the volume parameters
    pix_pec._bound_values = {40., 100., 350., 50., M_PI, 0.};
    pix_pec._type = decltype(bp)::type::e_cylinder;

    // Negative endcap portal: nec
    proto::portal<point3_container> bp_nec;
    proto::surface<point3_container> bp_s_nec;
    bp_s_nec._type = proto::surface<point3_container>::type::e_disc;
    bp_s_nec._sf_type = proto::surface<point3_container>::sf_type::e_portal;
    bp_s_nec._radii = {0., 40.};
    bp_s_nec._zparameters = {-400., 0.};

    std::vector<style::color> vColors = {
        style::color({{50, 50, 50}}), style::color({{255, 0, 0}}),
        style::color({{0, 255, 0}}), style::color({{0, 0, 255}})};
    for (auto& c : vColors) {
        c._opacity = 0.1_scalar;
    }

    // The portals
    // (1) joint nec
    proto::portal<point3_container> nec;
    nec._name = "joint_nec";
    proto::surface<point3_container> s_nec;
    s_nec._type = proto::surface<point3_container>::type::e_disc;
    s_nec._sf_type = proto::surface<point3_container>::sf_type::e_portal;
    s_nec._radii = {0., 100.};
    s_nec._zparameters = {-400., 0.};

    // Assign the surface
    nec._surface = s_nec;
    // nec links
    proto::portal<point3_container>::link nec_link_bp;
    nec_link_bp._start = {20., 0., -400.};
    nec_link_bp._end = {20., 0., -380.};
    nec_link_bp._link_index = 0u;
    proto::portal<point3_container>::link nec_link_pix_nec;
    nec_link_pix_nec._start = {70., 0., -400.};
    nec_link_pix_nec._end = {70., 0., -380.};
    nec_link_pix_nec._link_index = 1u;

    nec._volume_links = {nec_link_bp, nec_link_pix_nec};

    // (2) joint inner cylinder
    proto::portal<point3_container> ci;
    ci._name = "joint_c_i";
    proto::surface<point3_container> s_ci;
    s_ci._type = proto::surface<point3_container>::type::e_cylinder;
    s_ci._sf_type = proto::surface<point3_container>::sf_type::e_portal;
    s_ci._radii = {0., 40.};
    s_ci._zparameters = {0., 400.};

    // Assign the surface
    ci._surface = s_ci;

    proto::portal<point3_container>::link ci_link_bp;
    ci_link_bp._start = {40., 0., 0.};
    ci_link_bp._end = {20., 0., 0.};
    ci_link_bp._link_index = 0u;

    proto::portal<point3_container>::link ci_link_pix_nec;
    ci_link_pix_nec._start = {40., 0., -350.};
    ci_link_pix_nec._end = {60., 0., -350.};
    ci_link_pix_nec._link_index = 1u;

    proto::portal<point3_container>::link ci_link_pix_b;
    ci_link_pix_b._start = {40., 0., 0.};
    ci_link_pix_b._end = {60., 0., 0.};
    ci_link_pix_b._link_index = 2u;

    proto::portal<point3_container>::link ci_link_pix_pec;
    ci_link_pix_pec._start = {40., 0., 350.};
    ci_link_pix_pec._end = {60., 0., 350.};
    ci_link_pix_pec._link_index = 3u;

    ci._volume_links = {ci_link_bp, ci_link_pix_nec, ci_link_pix_b,
                        ci_link_pix_pec};

    // (3) joint pec
    proto::portal<point3_container> pec;
    pec._name = "joint_pec";
    proto::surface<point3_container> s_pec;
    s_pec._type = proto::surface<point3_container>::type::e_disc;
    s_pec._radii = {0., 100.};
    s_pec._zparameters = {400., 0.};

    // Assign the surface
    pec._surface = s_pec;
    // pec links
    proto::portal<point3_container>::link pec_link_bp;
    pec_link_bp._start = {20., 0., 400.};
    pec_link_bp._end = {20., 0., 380.};
    pec_link_bp._link_index = 0u;
    proto::portal<point3_container>::link pec_link_pix_pec;
    pec_link_pix_pec._start = {70., 0., 400.};
    pec_link_pix_pec._end = {70., 0., 380.};
    pec_link_pix_pec._link_index = 3u;

    pec._volume_links = {pec_link_bp, pec_link_pix_pec};

    // (4) joint inner cylinder
    proto::portal<point3_container> co;
    co._name = "joint_c_o";
    proto::surface<point3_container> s_co;
    s_co._type = proto::surface<point3_container>::type::e_cylinder;
    s_co._sf_type = proto::surface<point3_container>::sf_type::e_portal;
    s_co._radii = {0., 100.};
    s_co._zparameters = {0., 400.};

    // Assign the surface
    co._surface = s_co;

    proto::portal<point3_container>::link co_link_pix_nec;
    co_link_pix_nec._start = {100., 0., -350.};
    co_link_pix_nec._end = {80., 0., -350.};
    co_link_pix_nec._link_index = 1u;

    proto::portal<point3_container>::link co_link_pix_b;
    co_link_pix_b._start = {100., 0., 0.};
    co_link_pix_b._end = {80., 0., 0.};
    co_link_pix_b._link_index = 2u;

    proto::portal<point3_container>::link co_link_pix_pec;
    co_link_pix_pec._start = {100., 0., 350.};
    co_link_pix_pec._end = {80., 0., 350.};
    co_link_pix_pec._link_index = 3u;

    co._volume_links = {co_link_pix_nec, co_link_pix_b, co_link_pix_pec};

    // (5) pix_nec - pix_b
    proto::portal<point3_container> pix_nec_b;
    pix_nec_b._name = "pix_nec_b";
    proto::surface<point3_container> s_pix_nec_b;
    s_pix_nec_b._type = proto::surface<point3_container>::type::e_disc;
    s_pix_nec_b._sf_type = proto::surface<point3_container>::sf_type::e_portal;
    s_pix_nec_b._radii = {40., 100.};
    s_pix_nec_b._zparameters = {-300., 0.};

    // Assign the surface
    pix_nec_b._surface = s_pix_nec_b;

    proto::portal<point3_container>::link pix_nec_b_link_pix_nec;
    pix_nec_b_link_pix_nec._start = {70., 0., -300.};
    pix_nec_b_link_pix_nec._end = {70., 0., -320.};
    pix_nec_b_link_pix_nec._link_index = 1u;

    proto::portal<point3_container>::link pix_nec_b_link_pix_b;
    pix_nec_b_link_pix_b._start = {70., 0., -300.};
    pix_nec_b_link_pix_b._end = {70., 0., -280.};
    pix_nec_b_link_pix_b._link_index = 2u;

    pix_nec_b._volume_links = {pix_nec_b_link_pix_nec, pix_nec_b_link_pix_b};

    // (6) pix_b - pix_bec
    proto::portal<point3_container> pix_b_pec;
    pix_b_pec._name = "pix_b_pec";
    proto::surface<point3_container> s_pix_b_pec;
    s_pix_b_pec._type = proto::surface<point3_container>::type::e_disc;
    s_pix_b_pec._sf_type = proto::surface<point3_container>::sf_type::e_portal;
    s_pix_b_pec._radii = {40., 100.};
    s_pix_b_pec._zparameters = {300., 0.};

    // Assign the surface
    pix_b_pec._surface = s_pix_b_pec;

    proto::portal<point3_container>::link s_pix_b_pec_link_pix_b;
    s_pix_b_pec_link_pix_b._start = {70., 0., 300.};
    s_pix_b_pec_link_pix_b._end = {70., 0., 280.};
    s_pix_b_pec_link_pix_b._link_index = 2u;

    proto::portal<point3_container>::link s_pix_b_pec_link_pix_pec;
    s_pix_b_pec_link_pix_pec._start = {70., 0., 300.};
    s_pix_b_pec_link_pix_pec._end = {70., 0., 320.};
    s_pix_b_pec_link_pix_pec._link_index = 3u;

    pix_b_pec._volume_links = {s_pix_b_pec_link_pix_b,
                               s_pix_b_pec_link_pix_pec};

    bp._portals = {nec, ci, pec};
    pix_nec._portals = {nec, ci, co, pix_nec_b};
    pix_b._portals = {ci, co, pix_nec_b, pix_b_pec};
    pix_pec._portals = {pec, ci, co, pix_b_pec};

    // Get the volumes and colorize
    d._volumes = {bp, pix_nec, pix_b, pix_pec};
    for (auto& v : d._volumes) {
        v.colorize(vColors);
    }

    // Test the disc in z-r view
    svg::object d_zr = display::detector("cylinder_detector", d, views::z_r{});

    svg::file rfile_zr;
    rfile_zr.add_object(d_zr);

    std::ofstream rstream;
    rstream.open("test_meta_cylinder_detector_zr.svg");
    rstream << rfile_zr;
    rstream.close();

    // Let's adjust the view box as a test
    std::array<scalar, 4u> view_box = rfile_zr.view_box();
    view_box[0] = -10;
    view_box[2] *= 0.5;
    view_box[2] += 100.;
    rfile_zr.set_view_box(view_box);

    rstream.open("test_meta_cylinder_detector_zr_half.svg");
    rstream << rfile_zr;
    rstream.close();
}
