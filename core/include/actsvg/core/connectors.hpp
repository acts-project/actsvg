// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <fstream>
#include <map>
#include <vector>

#include "draw.hpp"
#include "svg.hpp"
#include "utils.hpp"

namespace actsvg {

namespace connectors {

// The type of connections that exist
enum type { e_highlight, e_associate_id, e_associate_info };

/** Helper method to connect objects
 *
 * @param sources_ are the source objects
 * @param targests_ are the target objects
 * @param s_t_connections_ are the connections from source to target
 * @param ls_ label the source or not (in case of multiple connection sheets)
 * @param on_off_ are the connection effects
 *
 * In case of e.g. a grid surface connection:
 * - the @param sources_ [in, out]are the grid tiles
 * - the @param targets_ [in, out] are the surfaces
 * - the @param s_t_connections are the indices which type of association
 * - the @param font_ is used for associate info action
 *
 **/
std::vector<svg::object> connect_action(
    std::vector<svg::object> &sources_, std::vector<svg::object> &targets_,
    const std::vector<std::vector<size_t>> &s_t_connections_, bool ls_ = true,
    const std::array<std::string, 2u> &on_off_ = {"mouseover", "mouseout"},
    const std::vector<connectors::type> &ct_ = {e_highlight, e_associate_info},
    const style::font &font_ = style::font());

struct name_extractor {
    /// Strips the name out of an object
    template <typename source_type>
    auto name(const source_type &s_) {
        return s_._name;
    }
};

struct id_extractor {
    /// Strips the id out of an object
    template <typename source_type>
    auto name(const source_type &s_) {
        return s_._id;
    }
};

/** Helper method to connect the appearance of an object to some other
 * objects
 *
 * @tparam source_objects_type the type of the source objects
 * @tparam source_name_getter_type the way the name is provided
 *
 * @param sources_ are the source objects
 * @param targests_ [in, out] are the target objects
 * @param s_t_connections_ an association which source belongs to which target
 *
 **/
template <typename source_objects_type,
          typename source_name_getter_type = name_extractor>
void connect_action(const source_objects_type &sources_,
                    std::vector<svg::object> &targets_,
                    const std::vector<std::vector<size_t>> &s_t_connections_) {

    // Ignore if the s_t_connections and sources don't match up
    if (sources_.size() != s_t_connections_.size()) {
        return;
    }

    // Loop over the sources and and register the connections
    for (auto [is, s] : utils::enumerate(sources_)) {
        // Isolate the name
        std::string sid = source_name_getter_type{}.name(s);

        // Get the target list and create the connection
        auto st_cs = s_t_connections_[is];

        for (auto it : st_cs) {
            // Ignore if not available
            if (it >= targets_.size()) {
                continue;
            }

            svg::object &t = targets_[it];
            t._attribute_map["display"] = "none";

            // Object information to appear
            svg::object on;
            on._tag = "animate";
            on._attribute_map["fill"] = "freeze";
            on._attribute_map["attributeName"] = "display";
            on._attribute_map["from"] = "none";
            on._attribute_map["to"] = "block";
            on._attribute_map["begin"] = sid + __d + "mouseout";

            svg::object off;

            off._tag = "animate";
            off._attribute_map["fill"] = "freeze";
            off._attribute_map["attributeName"] = "display";
            off._attribute_map["to"] = "none";
            off._attribute_map["from"] = "block";
            off._attribute_map["begin"] = sid + __d + "mouseover";

            // Store the animation
            t._sub_objects.push_back(off);
            t._sub_objects.push_back(on);
        }
    }
}

}  // namespace connectors
}  // namespace actsvg
