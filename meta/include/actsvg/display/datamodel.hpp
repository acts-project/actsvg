// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/proto/cluster.hpp"

namespace actsvg {

using namespace defaults;

namespace display {

/** Make a cluster  in 2D
 *
 * @param grid_ is the grid in which this cluster sits
 * @param id_ is the identification tag
 * @param cluster_ is the cluster to be drawn
 * @param fill_low_ is the low color band
 * @param fill_high is hte high color band
 * @param fill_m_ is the measurement fill color
 * @param expand_ is the focus expansion
 *
 * @return a cluster display and the focussed grid area
 *
 **/
template <size_t DIM>
static inline std::pair<svg::object, svg::object> cluster(
    const svg::object& grid_, const std::string& id_,
    const proto::cluster<DIM>& cluster_,
    const style::fill& fill_low_ = style::fill{style::color{{255, 255, 0}}},
    const style::fill& fill_high_ = style::fill{style::color{{250, 0, 0}}},
    const style::fill& fill_m_ = style::fill{style::color{{0, 0, 255}}},
    const std::array<unsigned int, 2>& expand_ = {2, 2}) {

    svg::object cluster_group;
    cluster_group._tag = "g";
    cluster_group._id = id_;

    // Grid indices
    unsigned int i_min = std::numeric_limits<unsigned int>::max();
    unsigned int i_max = std::numeric_limits<unsigned int>::min();
    unsigned int j_min = std::numeric_limits<unsigned int>::max();
    unsigned int j_max = std::numeric_limits<unsigned int>::min();

    // Cluster low/high data
    scalar low_data = std::numeric_limits<scalar>::max();
    scalar high_data = std::numeric_limits<scalar>::lowest();

    for (const auto& c : cluster_._channels) {
        low_data = std::min(low_data, c._data);
        high_data = std::max(high_data, c._data);
    }
    // Cluster data scaling
    scalar delta_data = high_data - low_data;

    int low_r = fill_low_._fc._rgb[0];
    int low_g = fill_low_._fc._rgb[1];
    int low_b = fill_low_._fc._rgb[2];
    int delta_r = (fill_high_._fc._rgb[0] - low_r);
    int delta_g = (fill_high_._fc._rgb[1] - low_g);
    int delta_b = (fill_high_._fc._rgb[2] - low_b);

    // Measurement, covariance, truth
    scalar t_x = 0.;
    scalar t_y = 0.;
    scalar t_r = 0.;
    scalar t_phi = 0.;

    scalar m_x = 0.;
    scalar m_y = 0.;
    scalar m_r = 0.;
    scalar m_phi = 0.;

    scalar d_r = 0.;
    scalar d_phi = 0.;

    scalar c_x = 0.;
    scalar c_y = 0.;
    scalar c_r = 0.;

    std::array<scalar, 2> x_range = {std::numeric_limits<scalar>::max(),
                                     std::numeric_limits<scalar>::lowest()};
    std::array<scalar, 2> y_range = x_range;
    std::array<scalar, 2> r_range = x_range;
    std::array<scalar, 2> phi_range = x_range;

    /// Helper function to generate a fill color
    ///
    /// @param channel_ the channel in question
    auto generate_fill =
        [&](const proto::channel<DIM>& channel_) -> style::fill {
        // Scale the colors for the filling
        scalar rel_delta_data = (channel_._data - low_data) / delta_data;

        int ch_r = low_r + rel_delta_data * delta_r;
        int ch_g = low_g + rel_delta_data * delta_g;
        int ch_b = low_b + rel_delta_data * delta_b;

        return style::fill{style::color{{ch_r, ch_g, ch_b}}};
    };

    if constexpr (DIM == 2) {

        if (cluster_._type == proto::cluster<DIM>::e_polar) {
            // Polar case, measurement and truth
            m_r = cluster_._measurement[DIM - 2];
            m_phi = cluster_._measurement[DIM - 1];

            d_r = cluster_._variance[DIM - 2];
            d_phi = cluster_._variance[DIM - 1];

            t_r = cluster_._measurement[DIM - 2];
            t_phi = cluster_._measurement[DIM - 1];

        } else {
            // Cartesian/Fan case
            m_x = cluster_._measurement[DIM - 2];
            m_y = cluster_._measurement[DIM - 1];
            c_x = cluster_._variance[DIM - 2];
            c_y = cluster_._variance[DIM - 1];
            t_x = cluster_._truth[DIM - 2];
            t_y = cluster_._truth[DIM - 1];

            // Correlation between -1 (-45 deg)and 1 (45 deg)
            scalar corr = cluster_._correlation;
            corr = corr < -1. ? -1. : (corr > 1. ? 1. : corr);
            c_r = -corr * 45;
        }
    }

    // Loop over the cluster channels and draw them
    for (const auto& c : cluster_._channels) {
        // Conventional indices
        size_t i = 0;
        size_t j = 0;
        // 2-Dimensional clusters: trivial
        if constexpr (DIM == 2) {
            i = c._cid[DIM - 2];
            j = c._cid[DIM - 1];
        }
        // 1-Dimensional cluster
        if constexpr (DIM == 1) {
            if (cluster_._coords[DIM - 1] == proto::cluster<DIM>::e_r or
                cluster_._coords[DIM - 1] == proto::cluster<DIM>::e_x) {
                i = c._cid[DIM - 1];
            } else {
                j = c._cid[DIM - 1];
            }
        }
        // Remember for the focus point
        i_min = i < i_min ? i : i_min;
        i_max = i_max < i ? i : i_max;
        j_min = j < j_min ? j : j_min;
        j_max = j_max < j ? j : j_max;

        // Get the corresponding grid tile (will be a copy already if
        // successful)
        std::string tile_id =
            grid_._id + "_" + std::to_string(i) + "_" + std::to_string(j);
        std::optional<svg::object> grid_tile = grid_.find_object(tile_id);

        if (grid_tile.has_value()) {
            // Get the tile
            svg::object tile = grid_tile.value();
            // Record the x/y/r/phi area
            x_range[0] = std::min(x_range[0], tile._x_range[0]);
            x_range[1] = std::max(x_range[1], tile._x_range[1]);
            y_range[0] = std::min(y_range[0], tile._y_range[0]);
            y_range[1] = std::max(y_range[1], tile._y_range[1]);
            r_range[0] = std::min(r_range[0], tile._r_range[0]);
            r_range[1] = std::max(r_range[1], tile._r_range[1]);
            phi_range[0] = std::min(phi_range[0], tile._phi_range[0]);
            phi_range[1] = std::max(phi_range[1], tile._phi_range[1]);

            auto points = tile._attribute_map.find("points");
            if (points != tile._attribute_map.end()) {
                svg::object cell;
                cell._tag = "polygon";
                cell._id = id_ + "_cell_" + std::to_string(i) + "_" +
                           std::to_string(j);
                cell._fill = generate_fill(c);
                cell._stroke = style::stroke();
                cell._attribute_map["points"] = points->second;
                // Add the cell object
                cluster_group.add_object(cell);
            }
        }
    }

    // Need to fill the remaining parameters for measurement & truth
    if constexpr (DIM == 1) {
        if (cluster_._coords[DIM - 1] == proto::cluster<DIM>::e_r) {
            m_r = cluster_._measurement[DIM - 1];
            t_r = cluster_._truth[DIM - 1];

            d_r = cluster_._variance[DIM - 1];
            d_phi = std::abs(phi_range[1] - phi_range[0]) / std::sqrt(12.);

            m_phi = 0.5 * (phi_range[0] + phi_range[1]);
            t_phi = m_phi;
        } else if (cluster_._coords[DIM - 1] == proto::cluster<DIM>::e_phi) {
            m_phi = cluster_._measurement[DIM - 1];
            t_phi = cluster_._truth[DIM - 1];

            d_r = std::abs(r_range[1] - r_range[0]) / std::sqrt(12.);
            d_phi = cluster_._variance[DIM - 1];

            m_r = 0.5 * (r_range[0] + r_range[1]);
            t_r = m_r;
        } else if (cluster_._coords[DIM - 1] == proto::cluster<DIM>::e_x) {
            m_x = cluster_._measurement[DIM - 1];
            t_x = cluster_._truth[DIM - 1];
            c_x = cluster_._variance[DIM - 1];
            c_y = std::abs(y_range[1] - y_range[0]) / std::sqrt(12.);
            m_y = 0.5 * (y_range[0] + y_range[1]);
            t_y = m_y;

            // Correlation between -1 (-45 deg)and 1 (45 deg)
            scalar corr = cluster_._correlation;
            corr = corr < -1. ? -1. : (corr > 1. ? 1. : corr);
            c_r = -corr * 45;

        } else {
            m_y = cluster_._measurement[DIM - 1];
            t_y = cluster_._truth[DIM - 1];
            c_y = cluster_._variance[DIM - 1];
            c_x = std::abs(x_range[1] - x_range[0]) / std::sqrt(12.);
            m_x = 0.5 * (x_range[0] + x_range[1]);
            t_x = m_x;

            // Correlation between -1 (-45 deg)and 1 (45 deg)
            scalar corr = cluster_._correlation;
            corr = corr < -1. ? -1. : (corr > 1. ? 1. : corr);
            c_r = corr * 45;
        }
    }

    // Finally, if the cluster is polar, transform
    if (cluster_._type == proto::cluster<DIM>::e_polar) {

        m_x = m_r * std::cos(m_phi);
        m_y = m_r * std::sin(m_phi);

        scalar cos_phi = std::cos(m_phi);
        scalar sin_phi = std::sin(m_phi);

        c_x = std::abs(cos_phi * d_r - m_r * sin_phi * d_phi);
        c_y = std::abs(sin_phi * d_r + m_r * cos_phi * d_phi);

        c_r = m_phi * 180 / M_PI;

        t_x = t_r * std::cos(t_phi);
        t_y = t_r * std::sin(t_phi);
    }

    // Add an error ellipse
    style::fill fill_c = fill_m_;
    fill_c._fc._opacity = 0.5;
    style::stroke stroke_c;
    style::transform t_c;
    t_c._rot = {c_r, m_x, m_y};
    cluster_group.add_object(
        draw::ellipse("c", {m_x, m_y}, {c_x, c_y}, fill_c, stroke_c, t_c));

    // Add the marker
    style::marker marker{"o"};
    marker._fill = fill_m_;
    cluster_group.add_object(draw::marker("m", {m_x, m_y}, marker));

    if (cluster_._mc) {
        style::marker truth_marker{"o"};
        cluster_group.add_object(draw::marker("t", {t_x, t_y}, truth_marker));
    }

    // Get the focussed grid
    svg::object grid_area;
    grid_area._tag = "g";
    grid_area._id = id_ + std::string("_focussed_grid");

    // Expand
    i_min = (i_min >= expand_[0]) ? i_min - expand_[0] : 0;
    i_max += expand_[0];
    j_min = (j_min >= expand_[1]) ? j_min - expand_[1] : 0;
    j_max += expand_[1];

    for (unsigned int i = i_min; i <= i_max; ++i) {
        for (unsigned int j = j_min; j < j_max; ++j) {
            std::string tile_id =
                grid_._id + "_" + std::to_string(i) + "_" + std::to_string(j);
            std::optional<svg::object> grid_tile = grid_.find_object(tile_id);
            if (grid_tile.has_value()) {
                grid_area.add_object(grid_tile.value());
            }
        }
    }

    return {cluster_group, grid_area};
}

}  // namespace display
}  // namespace actsvg
