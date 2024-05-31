// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/proto/portal.hpp"

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "actsvg/display/geometry.hpp"
#include "actsvg/proto/surface.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

TEST(proto, cylinder_portal) {

    proto::portal<point3_container> p;

    ASSERT_TRUE(p._name == "unnamed_portal");
    ASSERT_TRUE(p._volume_links.empty());

    proto::surface<point3_container> s;
    s._type = proto::surface<point3_container>::type::e_cylinder;
    s._sf_type = proto::surface<point3_container>::sf_type::e_portal;
    s._radii = {0., 100.};
    s._zparameters = {0., 200.};

    // Assign the surface
    p._surface = s;

    proto::portal<point3_container>::link link_to_self;

    link_to_self._start = {0., 100., 0.};
    link_to_self._end = {0., 90., 0.};

    proto::portal<point3_container>::link link_to_outside;
    link_to_outside._start = {0., 100., 0.};
    link_to_outside._end = {0., 110., 0.};
    link_to_outside._stroke = style::stroke({style::color{{255, 0, 0}}});
    link_to_outside._end_marker._fill =
        style::fill{{link_to_outside._stroke._sc}};
    link_to_outside._end_marker._stroke = link_to_outside._stroke;

    p._volume_links = {link_to_self, link_to_outside};

    // Test the portal in x-y view
    svg::object cp_xy = display::portal("cylinder_portal", p, views::x_y{});

    svg::file rfile_xy;
    rfile_xy.add_object(cp_xy);

    std::ofstream rstream;
    rstream.open("test_meta_portal_cylinder_xy.svg");
    rstream << rfile_xy;
    rstream.close();

    // Test the disc in z-r view
    svg::object cp_zr = display::portal("cylinder_portal", p, views::z_r{});

    svg::file rfile_zr;
    rfile_zr.add_object(cp_zr);

    rstream.open("test_meta_portal_cylinder_zr.svg");
    rstream << rfile_zr;
    rstream.close();
}

TEST(proto, full_disc_portal) {

    proto::portal<point3_container> p;

    ASSERT_TRUE(p._name == "unnamed_portal");
    ASSERT_TRUE(p._volume_links.empty());

    proto::surface<point3_container> s;
    s._type = proto::surface<point3_container>::type::e_disc;
    s._sf_type = proto::surface<point3_container>::sf_type::e_portal;
    s._radii = {10., 100.};
    s._zparameters = {50., 0.};
    s._fill = defaults::__nn_fill;

    proto::portal<point3_container>::link link_to_self;

    link_to_self._start = {0., 55., 50.};
    link_to_self._end = {0., 55., 40.};

    proto::portal<point3_container>::link link_to_outside;
    link_to_outside._start = {0., 55., 50.};
    link_to_outside._end = {0., 55., 60.};
    link_to_outside._stroke = style::stroke({style::color{{255, 0, 0}}});
    link_to_outside._end_marker._fill =
        style::fill{{link_to_outside._stroke._sc}};
    link_to_outside._end_marker._stroke = link_to_outside._stroke;

    p._volume_links = {link_to_self, link_to_outside};

    // Test the portal in x-y view
    svg::object dp_xy = display::portal("full_disc_portal", p, views::x_y{});

    svg::file rfile_xy;
    rfile_xy.add_object(dp_xy);

    std::ofstream rstream;
    rstream.open("test_meta_portal_full_disc_xy.svg");
    rstream << rfile_xy;
    rstream.close();

    // Test the disc in z-r view
    svg::object dp_zr = display::portal("full_disc_portal", p, views::z_r{});

    svg::file rfile_zr;
    rfile_zr.add_object(dp_zr);

    rstream.open("test_meta_portal_full_disc_zr.svg");
    rstream << rfile_zr;
    rstream.close();
}
