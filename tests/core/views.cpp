// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core/views.hpp"

#include <gtest/gtest.h>

#include <array>
#include <fstream>

#include "actsvg/core/defs.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;

TEST(views, xy) {
    std::vector<point3> p3s = {{1., 2., 0.}, {3., 4., 0.}, {7., 8., 0.}};

    views::x_y x_y_view;
    auto c = x_y_view.path(p3s);

    ASSERT_TRUE(c.size() == 3u);

    views::contour exptected = {{1., 2.}, {3., 4.}, {7., 8.}};
    ASSERT_TRUE(exptected == c);
}

TEST(views, zr) {
    std::vector<point3> p3s = {{1., 0., 3.}, {0., 4., 0.}, {2., 2., 2.}};

    views::z_r z_r_view;
    auto c = z_r_view.path(p3s);

    ASSERT_TRUE(c.size() == 3u);

    views::contour exptected = {
        {3., 1.}, {0., 4.}, {2., static_cast<scalar>(std::sqrt(8.))}};
    ASSERT_TRUE(exptected == c);
}

TEST(views, zphi) {
    std::vector<point3> p3s = {{1., 0., 3.}, {0., 4., 0.}, {2., 2., 2.}};

    views::z_phi z_phi_view;
    auto c = z_phi_view.path(p3s);

    ASSERT_TRUE(c.size() == 3u);

    views::contour exptected = {{3., 0.}, {0., 0.5 * pi}, {2., 0.25 * pi}};
    ASSERT_TRUE(exptected == c);
}

TEST(views, zrphi) {
    std::vector<point3> p3s = {{1., 0., 3.}, {0., 4., 0.}, {2., 2., 2.}};

    views::z_rphi z_rphi_view;
    z_rphi_view._fixed_r = 2.;
    auto c = z_rphi_view.path(p3s);

    ASSERT_TRUE(c.size() == 3u);

    views::contour exptected = {{3., 0.}, {0., 1. * pi}, {2., 0.5 * pi}};
    ASSERT_TRUE(exptected == c);
}
