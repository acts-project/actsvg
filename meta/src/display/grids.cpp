// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/display/grids.hpp"

#include <algorithm>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/proto/grid.hpp"

namespace actsvg {

using namespace defaults;

namespace display {
svg::object grid(const std::string& id_, const proto::grid& g_) {

    // The grid object
    svg::object g = svg::object::create_group(id_);

    // Get the tiles
    std::vector<svg::object> grid_tiles;

    // Generate the grid
    if (g_._type == proto::grid::type::e_r_phi) {
        grid_tiles = draw::tiled_polar_grid(id_, g_._edges_0, g_._edges_1,
                                            g_._fill, g_._stroke)
                         ._sub_objects;

    } else {
        // cartesian view grid
        std::vector<scalar> edges_1 = g_._edges_1;
        if (g_._type == proto::grid::type::e_z_phi) {
            std::for_each(edges_1.begin(), edges_1.end(),
                          [&](auto& e) { e *= g_._reference_r; });
        }
        grid_tiles = draw::tiled_cartesian_grid(id_, g_._edges_0, edges_1,
                                                g_._fill, g_._stroke)
                         ._sub_objects;
    }

    if (not g_._connections.empty()) {
        std::vector<svg::object> local_connections = g_._connections;
        connectors::connect_action<decltype(grid_tiles),
                                   connectors::id_extractor>(
            grid_tiles, local_connections, g_._connection_associations);
        g.add_objects(local_connections);
    } else {
        g.add_objects(g_._connections);
    }
    g.add_objects(grid_tiles);

    return g;
}

}  // namespace display
}  // namespace actsvg
