// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "actsvg/core/defs.hpp"
#include "actsvg/core/svg.hpp"

namespace actsvg {

namespace test {

/** This method constructs a playground without using
 * the helper classes in the core module
 * in oder to check consistency of the objects and
 * transforms.
 *
 * @param llc_ the left lower corner
 * @param ruc_ the right upper corner
 *
 **/
static inline svg::object playground(
    const point2& llc_, const point2& ruc_,
    const std::vector<scalar>& ticks_x_ = {50., 100.},
    const std::vector<scalar>& ticks_y_ = {50., 100.}) {

    svg::object pg = svg::object::create_group("playground");
    pg._x_range = {llc_[0], ruc_[0]};
    pg._y_range = {-ruc_[1], -llc_[1]};

    /// Add a playground
    scalar half = 0.5;
    scalar cx = half * (llc_[0] + ruc_[0]);
    scalar cy = -half * (llc_[0] + ruc_[0]);
    scalar wx = std::abs(ruc_[0] - llc_[0]);
    scalar hy = std::abs(ruc_[0] - llc_[0]);

    svg::object field;
    field._tag = "rect";
    field._attribute_map["x"] = std::to_string(cx - half * wx);
    field._attribute_map["y"] = std::to_string(cy - half * hy);
    field._attribute_map["width"] = std::to_string(wx);
    field._attribute_map["height"] = std::to_string(hy);
    field._attribute_map["style"] =
        "fill:grey;stroke:grey;stroke-width:1;fill-opacity:0.1;stroke-opacity:"
        "0.5";

    pg._sub_objects.push_back(field);

    // x and y axis
    svg::object x_axis;
    x_axis._tag = "line";
    x_axis._attribute_map["x1"] = std::to_string(0.9 * llc_[0]);
    x_axis._attribute_map["y1"] = "0.";
    x_axis._attribute_map["x2"] = std::to_string(0.9 * ruc_[0]);
    x_axis._attribute_map["y2"] = "0.";
    x_axis._attribute_map["style"] = "stroke:black;stroke-width:1;";
    pg._sub_objects.push_back(x_axis);

    svg::object y_axis;
    y_axis._tag = "line";
    y_axis._attribute_map["x1"] = "0.";
    y_axis._attribute_map["y1"] = std::to_string(-0.9 * llc_[1]);
    y_axis._attribute_map["x2"] = "0.";
    y_axis._attribute_map["y2"] = std::to_string(-0.9 * ruc_[1]);
    y_axis._attribute_map["style"] = "stroke:black;stroke-width:1;";
    pg._sub_objects.push_back(y_axis);

    // ticks with label
    for (const auto& tx : ticks_x_) {
        svg::object t_tick;
        t_tick._tag = "line";
        t_tick._attribute_map["x1"] = std::to_string(tx);
        t_tick._attribute_map["y1"] = std::to_string(0.005 * hy);
        t_tick._attribute_map["x2"] = std::to_string(tx);
        t_tick._attribute_map["y2"] = std::to_string(-0.005 * hy);
        t_tick._attribute_map["style"] = "stroke:black;stroke-width:1;";
        pg._sub_objects.push_back(t_tick);
        // add the label
        svg::object t_text;
        t_text._tag = "text";
        t_text._attribute_map["x"] = std::to_string(tx - 0.01 * wx);
        t_text._attribute_map["y"] = std::to_string(0.025 * hy);
        t_text._attribute_map["font-family"] = "Arial";
        t_text._fill = style::fill{{{0, 0, 0}}};
        t_text._field = {std::to_string(static_cast<int>(tx))};
        pg._sub_objects.push_back(t_text);
    }

    for (const auto& ty : ticks_y_) {
        svg::object t_tick;
        t_tick._tag = "line";
        t_tick._attribute_map["x1"] = std::to_string(-0.005 * wx);
        t_tick._attribute_map["y1"] = std::to_string(-ty);
        t_tick._attribute_map["x2"] = std::to_string(0.005 * wx);
        t_tick._attribute_map["y2"] = std::to_string(-ty);
        t_tick._attribute_map["style"] = "stroke:black;stroke-width:1;";
        pg._sub_objects.push_back(t_tick);
        // add the label
        svg::object t_text;
        t_text._tag = "text";
        t_text._attribute_map["x"] = std::to_string(0.01 * wx);
        t_text._attribute_map["y"] = std::to_string(-ty + 0.005 * hy);
        t_text._attribute_map["font-family"] = "Arial";
        t_text._fill = style::fill{{{0, 0, 0}}};
        t_text._field = {std::to_string(static_cast<int>(ty))};
        pg._sub_objects.push_back(t_text);
    }

    return pg;
}

}  // namespace test

}  // namespace actsvg
