// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include "actsvg/actsvg.hpp"
#include "actsvg/display/sheets.hpp"
#include "actsvg/proto/surface.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

TEST(display, surface_sheet) {

    proto::surface<point3_container> trapezoid;
    trapezoid._vertices = {{-8.5, -34, 0.}, {8.5, -34, 0.}, {14.5, 34., 0.}, {-14.5, 34., 0.}};

    svg::object surface_sheet = display::surface_sheet(trapezoid);
    svg::file surface_file;
    surface_file.add_object(surface_sheet);

    // Write out the file
    std::ofstream sout;
    sout.open("trapezoid_sheet.svg");
    sout << surface_file;
    sout.close();
}