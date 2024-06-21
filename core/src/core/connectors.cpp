// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core/connectors.hpp"

#include <cmath>
#include <fstream>
#include <map>
#include <vector>

#include "actsvg/core/draw.hpp"
#include "actsvg/core/svg.hpp"
#include "actsvg/core/utils.hpp"

namespace actsvg {

namespace connectors {
std::vector<svg::object> connect_action(
    std::vector<svg::object> &sources_, std::vector<svg::object> &targets_,
    const std::vector<std::vector<size_t>> &s_t_connections_, bool ls_,
    const std::array<std::string, 2u> &on_off_,
    const std::vector<connectors::type> &ct_, const style::font &font_) {

    std::vector<svg::object> additional_connections;

    for (auto [s, ts] : utils::enumerate(s_t_connections_)) {
        if (s < sources_.size()) {
            // The source object group
            auto &sog = sources_[s];
            // A connections string
            std::vector<std::string> id_assoc_text;
            std::vector<std::string> info_assoc_text;
            // Continue if you do not have a source identification
            if (sog._id.empty()) {
                continue;
            }
            // Associate the id tags, i.e. labl the source
            if (ls_) {
                id_assoc_text.push_back("Source: ");
                id_assoc_text.push_back("* " + sog._id);
                // Associate the auxiliary information
                info_assoc_text.push_back("Source: ");
                for (const auto &sai : sog._aux_info) {
                    info_assoc_text.push_back(sai);
                }
            }
            // Enumerate over the connections
            for (auto [it, t] : utils::enumerate(ts)) {
                if (t < targets_.size()) {
                    // Get the target (object group)
                    auto &tog = targets_[t];
                    // remember the associations
                    if (it == 0) {
                        id_assoc_text.push_back("Target:");
                        info_assoc_text.push_back("Target:");
                    }
                    id_assoc_text.push_back("* " + tog._id);
                    for (const auto &tai : tog._aux_info) {
                        info_assoc_text.push_back(tai);
                    }
                    // Highlight connection
                    if (std::find(ct_.begin(), ct_.end(), e_highlight) !=
                        ct_.end()) {
                        svg::object on_off;
                        on_off._tag = "set";
                        on_off._attribute_map["attributeName"] = "fill";
                        on_off._attribute_map["begin"] =
                            sog._id + __d + on_off_[0];
                        on_off._attribute_map["end"] =
                            sog._id + __d + on_off_[1];
                        // Stroke and fill sterile
                        on_off._stroke._sterile = true;
                        on_off._fill._sterile = true;
                        // If the object has a use object, the connection goes
                        // to the use object and not the the top object
                        bool connection_done = false;
                        for (auto &sob_tog : tog._sub_objects) {
                            if (sob_tog._tag == "use") {
                                on_off._attribute_map["to"] =
                                    style::rgb_attr(sob_tog._fill._fc._hl_rgb);
                                sob_tog._sub_objects.push_back(on_off);
                                connection_done = true;
                                break;
                            }
                        }
                        // On off tests
                        if (not connection_done) {
                            on_off._attribute_map["to"] =
                                style::rgb_attr(tog._fill._fc._hl_rgb);
                            tog._sub_objects.push_back(on_off);
                        }
                    }
                }
            }
            // Assocate as text
            if (std::find(ct_.begin(), ct_.end(), e_associate_id) !=
                ct_.end()) {
                // Associate id connection
                auto c_text = draw::connected_text(
                    sog._id + "_id_associations",
                    {sog._barycenter[0], sog._barycenter[1]}, id_assoc_text,
                    font_, style::transform(), sog);
                additional_connections.push_back(c_text);
            } else if (std::find(ct_.begin(), ct_.end(), e_associate_info) !=
                       ct_.end()) {
                // Associate info connection
                auto c_text = draw::connected_text(
                    sog._id + "_info_associations",
                    {sog._barycenter[0], sog._barycenter[1]}, info_assoc_text,
                    font_, style::transform(), sog);
                additional_connections.push_back(c_text);
            }
        }
    }
    return additional_connections;
}
}  // namespace connectors
}  // namespace actsvg
