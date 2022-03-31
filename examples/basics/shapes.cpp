// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "actsvg/actsvg.hpp"

namespace po = boost::program_options;
using namespace actsvg;

int main(int argc, char* argv[]) {

    int r = 200;
    int g = 150;
    int b = 50;

    try {
        po::options_description desc("Allowed options");
        // clang-format off
        desc.add_options()
            ("help", "produce help message")
            ("r",po::value<int>(&r),"Fill color R in (R,G,B)")
            ("g",po::value<int>(&r),"Fill color G in (R,G,B)")
            ("b",po::value<int>(&r),"Fill color B in (R,G,B)");

        // clang-format on
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help") != 0u) {
            std::cout << desc << std::endl;
            return 0;
        }
    } catch (std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }

    // Some definitions
    using point3 = std::array<scalar, 3>;
    views::x_y x_y_view;

    style::fill fill_style({r, g, b});
    fill_style._fc._highlight = {"mouseover", "mouseout"};
    fill_style._fc._hl_rgb = {0, 255, 0};
    fill_style._fc._opacity = 0.5;

    style::stroke stroke_style({r, g, b});

    // Make a rectangle shape
    //
    std::vector<point3> rectangle = {{-100., -200., 0.},
                                     {100., -200., 0.},
                                     {100., 200., 0.},
                                     {-100., 200., 0.}};
    auto rectangle_2d = x_y_view(rectangle);

    auto rectangle_svg =
        draw::polygon(rectangle_2d, "r0", fill_style, stroke_style);
    svg::file rectangle_file;
    rectangle_file._objects.push_back(rectangle_svg);

    auto x_y_axes_svg = draw::x_y_axes({-150, 150}, {-250, 250});
    rectangle_file._objects.insert(rectangle_file._objects.begin(),
                                   x_y_axes_svg.begin(), x_y_axes_svg.end());

    std::ofstream rectangle_stream;
    rectangle_stream.open("basic_rectangle.svg");
    rectangle_stream << rectangle_file;
    rectangle_stream.close();
}
