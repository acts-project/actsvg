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

#include "actsvg/core.hpp"

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

    style::fill fill_style{{{r, g, b}}};
    fill_style._fc._highlight = {"mouseover", "mouseout"};
    fill_style._fc._hl_rgb = {0, 255, 0};
    fill_style._fc._opacity = 0.5;

    style::stroke stroke_style{{{r, g, b}}};
    style::stroke stroke_black = style::stroke();

    // Make a rectangle shape
    //
    std::vector<point3> rectangle = {{-100., -200., 0.},
                                     {100., -200., 0.},
                                     {100., 200., 0.},
                                     {-100., 200., 0.}};
    auto rectangle_2d = x_y_view(rectangle);

    // rectangle svg object
    auto rectangle_svg =
        draw::polygon("r0", rectangle_2d, fill_style, stroke_style);
    svg::file rectangle_file;
    rectangle_file.add_object(rectangle_svg);

    auto x_y_a =
        draw::x_y_axes("xy", {-150, 150}, {-250, 250}, stroke_black, "x", "y");
    rectangle_file.add_object(x_y_a);

    // measure labeling
    auto measure_marker = style::marker({"|<"});
    auto measure_hlx =
        draw::measure("hlx", {0, 210}, {100., 210}, stroke_black,
                      measure_marker, style::font(), "hx", {50., 220.});
    auto measure_hly =
        draw::measure("hly", {110, 0}, {110., 200}, stroke_black,
                      measure_marker, style::font(), "hy", {120., 50.});
    rectangle_file.add_object(measure_hly);
    rectangle_file.add_object(measure_hlx);

    std::ofstream rectangle_stream;
    rectangle_stream.open("basic_rectangle.svg");
    rectangle_stream << rectangle_file;
    rectangle_stream.close();

    // Make a trapezoid shape
    std::vector<point3> trapezoid = {{-50., -200., 0.},
                                     {50., -200., 0.},
                                     {100., 200., 0.},
                                     {-100, 200., 0.}};

    auto trapezoid_2d = x_y_view(trapezoid);

    // rectangle svg object
    auto trapezoid_svg =
        draw::polygon("t0", trapezoid_2d, fill_style, stroke_style);
    svg::file trapezoid_file;
    trapezoid_file.add_object(trapezoid_svg);
    trapezoid_file.add_object(x_y_a);

    auto measure_hlx_min =
        draw::measure("hlx_min", {0, -210}, {50., -210}, stroke_black,
                      measure_marker, measure_marker, "hx_min", style::font(), {25., -220.} );

    auto measure_hlx_max =
        draw::measure("hlx_max", {0, 210}, {100., 210}, stroke_black,
                      measure_marker, measure_marker, "hx_max", style::font(), {50.,220.});

    trapezoid_file.add_object(measure_hlx_min);
    trapezoid_file.add_object(measure_hlx_max);
    trapezoid_file.add_object(measure_hly);

    std::ofstream trapezoid_stream;
    trapezoid_stream.open("basic_trapezoid.svg");
    trapezoid_stream << trapezoid_file;
    trapezoid_stream.close();
}
