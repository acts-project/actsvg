// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <vector>

#include "svg.hpp"

namespace actsvg {

namespace connectors {

/** Helper method to connect objects
 *
 * @param sources_ are the source objects
 * @param targests_ are the target objects
 * @param s_t_connections are the connections from source to target
 * @param on_off_ are the connection effects
 **/
static void connect_objects(
    std::vector<svg::object> &sources_, std::vector<svg::object> &targets_,
    const std::vector< std::vector<size_t> > &s_t_connections,
    const std::array<std::string, 2u> &on_off_ = {"mouseover", "mouseout"}) {

    for (auto [s, ts] : utils::enumerate(s_t_connections)) {
        if (s < sources_.size()) {
            auto sog = sources_[s];
            // Continue if you do not have a source identification
            if (sog._id.empty()) {
                continue;
            }
            // Make the connections
            for (auto t : ts) {
                if (t < targets_.size()) {
                    auto &tog = targets_[t];

                    // Highlight it
                    svg::object grid_on;
                    grid_on._tag = "animate";
                    grid_on._attribute_map["fill"] = "freeze";
                    grid_on._attribute_map["attributeName"] = "fill";
                    grid_on._attribute_map["begin"] =
                        sog._id + __d + on_off_[0];
                    grid_on._attribute_map["dur"] = "0.05s";
                    grid_on._attribute_map["to"] =
                        style::rgb_attr(tog._fill._fc._hl_rgb);
                    grid_on._attribute_map["from"] =
                        style::rgb_attr(tog._fill._fc._rgb);
                    // Back to normal
                    svg::object grid_off;
                    grid_off._tag = "animate";
                    grid_off._attribute_map["fill"] = "freeze";
                    grid_off._attribute_map["attributeName"] = "fill";
                    grid_off._attribute_map["begin"] =
                        sog._id + __d + on_off_[1];
                    grid_off._attribute_map["dur"] = "0.05s";
                    grid_off._attribute_map["from"] =
                        style::rgb_attr(tog._fill._fc._hl_rgb);
                    grid_off._attribute_map["to"] =
                        style::rgb_attr(tog._fill._fc._rgb);

                    // Record them
                    tog._sub_objects.push_back(grid_on);
                    tog._sub_objects.push_back(grid_off);
                }
            }
        }
    }
}
}  // namespace connectors
}  // namespace actsvg
