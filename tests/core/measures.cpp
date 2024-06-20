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

TEST(draw, measure) {

    // Set a playground
    auto pg = test::playground({-400, -400}, {400, 400});

    // Measure
    auto m0 = draw::measure("m0", {200, 10}, {300, 10});

    // Measure with text
    auto m0_t = draw::measure("m0_t", {100, 100}, {400, 400}, style::stroke(),
                              style::marker({"|<<"}), style::marker({"|<<"}),
                              style::font(), "300 mm", {270., 270.});

    scalar phi_min = -0.25;
    scalar phi_max = 0.75;
    scalar r = 75.;

    point2 start = {r * std::cos(phi_min), r * std::sin(phi_min)};
    point2 end = {r * std::cos(phi_max), r * std::sin(phi_max)};

    // Add the line
    style::marker arc_marker = style::marker{"|<<"};
    arc_marker._fill._fc = style::color{{255, 0, 0}};
    arc_marker._size = 10;
    arc_marker._stroke = style::stroke{style::color{{255, 0, 0}}};
    auto m_arc = draw::arc_measure("m_arc", r, start, end,
                                   style::stroke{style::color{{255, 0, 0}}, 2},
                                   style::marker(), arc_marker);

    point2 start2 = {static_cast<scalar>(1.5 * r * std::cos(phi_min)),
                     static_cast<scalar>(1.5 * r * std::sin(phi_min))};
    point2 end2 = {static_cast<scalar>(1.5 * r * std::cos(phi_max)),
                   static_cast<scalar>(1.5 * r * std::sin(phi_max))};

    // Add the line
    style::marker start_marker = style::marker{"|"};
    start_marker._fill._fc = style::color{{255, 0, 0}};
    start_marker._size = 10;
    start_marker._stroke = style::stroke{style::color{{255, 0, 0}}};
    auto m_arc2 = draw::arc_measure(
        "m_arc", static_cast<scalar>(1.5 * r), start2, end2,
        style::stroke{style::color{{255, 0, 0}}, 2}, start_marker, arc_marker);

    svg::file mfile;
    mfile.add_object(pg);
    mfile.add_object(m0);
    mfile.add_object(m0_t);
    mfile.add_object(m_arc);
    mfile.add_object(m_arc2);

    std::ofstream tstream;
    tstream.open("test_core_measures.svg");
    tstream << mfile;
    tstream.close();
}
