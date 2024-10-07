// This file is part of the actsvg packge.
//
// Copyright (C) 2024 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string>

#include "actsvg/core.hpp"
#include "actsvg/display/grids.hpp"
#include "actsvg/proto/grid.hpp"
#include "actsvg/proto/material.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {

using namespace defaults;

namespace display {
svg::object surface_material(const std::string& id_,
                             const proto::surface_material& m_) {

    svg::object m = svg::object::create_group(id_);
    m._sterile = true;

    // Create the grid object
    std::string g_tag = "_grid";
    svg::object g = display::grid(id_ + g_tag, m_._grid);

    if (m_._material_matrix.empty()) {
        return g;
    }

    std::vector<svg::object> m_info;

    for (auto& bin : g._sub_objects) {
        std::string bin_str = bin._id;
        // Remove the grid id from the string
        bin_str.erase(0u, g._id.size() + 1);
        // Take the sub string after the last underscore
        auto last_underscore = bin_str.find_last_of("_");
        // Take the two sub strings left and right of the last underscore
        std::size_t b0 = std::stoul(bin_str.substr(0u, last_underscore));
        std::size_t b1 =
            std::stoul(bin_str.substr(last_underscore + 1u, bin_str.size()));
        // Associated material
        const auto& slab = m_._material_matrix[b1][b0];
        // For the moment we support X0 view
        scalar t_X0 = slab._properties[5u] / slab._properties[0u];
        scalar t_X0_rel =
            (t_X0 - m_._material_ranges[0][0]) /
            (m_._material_ranges[0][1] - m_._material_ranges[0][0]);

        // create the fill color from it
        std::vector<std::string> info;
        // Make a connected text box
        std::string info_bin =
            "material, bin " + std::to_string(b0) + " " + std::to_string(b1);
        info.push_back(info_bin);
        // Fill the bin with the color
        if (slab._properties[5] == 0.) {
            bin._fill = __w_fill;
            info.push_back("- No material");
        } else {
            bin._fill = style::fill{
                style::color{m_._gradient.rgb_from_scale(t_X0_rel)}};
            // Add the info to the bin
            std::string info_tX0 = "- t/X0: " + std::to_string(t_X0);
            std::string info_tL0 =
                "- t/L0: " +
                std::to_string(slab._properties[5u] / slab._properties[1u]);
            info.push_back(info_tX0);
            info.push_back(info_tL0);
        }

        // Connect the text info about the bin to it
        svg::object m_info_text =
            draw::connected_text(bin_str + "_material_info", m_._info_pos, info,
                                 m_._info_font, style::transform{}, bin);
        m_info.push_back(m_info_text);
    }
    // Draw the gradient box
    using stop_value = std::tuple<style::gradient::stop, scalar>;
    std::vector<stop_value> stops;

    // Get the min and max of the material
    scalar tX0_min = m_._material_ranges[0][0];
    scalar tX0_max = m_._material_ranges[0][1];

    for (const auto& s : m_._gradient._stops) {
        stops.push_back({s, tX0_min + (tX0_max - tX0_min) * s.first});
    }

    svg::object gbox = draw::gradient_box(
        id_ + "_gradient_box", m_._gradient_pos, m_._gradient_box, stops,
        m_._gradient_label, m_._gradient_stroke, m_._gradient_font);

    m.add_object(g);
    m.add_objects(m_info);
    m.add_object(gbox);

    return m;
}

}  // namespace display
}  // namespace actsvg
