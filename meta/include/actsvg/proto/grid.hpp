// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>

#include "actsvg/core/defs.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {

namespace proto {

/** A proto grid class to describe the grid setup
 *
 * For convenience, it is forced to be 2-dimensional, for 1-dim
 * grid descriptions, only provide 2 eges
 *
 */
struct grid {

    /** Type of grid, enum defintion */
    enum type { e_x_y = 0, e_r_phi = 1, e_z_phi = 2 };

    /// Grid name
    std::string _name = "grid";

    /// Chose the type
    type _type = e_r_phi;

    /// The edges in the two given directions, loc0
    std::vector<scalar> _edges_0 = {};

    /// The edges in the two given directions, loc1
    std::vector<scalar> _edges_1 = {};

    /// Reference r for drawing (e.g. in case of phi axis)
    scalar _reference_r = 0.;

    /// Bin identifications
    std::vector<std::string> _bin_ids = {};

    /// Connections
    std::vector<svg::object> _connections;
    /// Connection types
    std::vector<connectors::type> _connection_types = {};
    /// Associations
    std::vector<std::vector<size_t>> _connection_associations = {};

    /// Stroke
    style::fill _fill = defaults::__g_fill;

    /// Stroke
    style::stroke _stroke = defaults::__g_stroke;
};

}  // namespace proto

}  // namespace actsvg
