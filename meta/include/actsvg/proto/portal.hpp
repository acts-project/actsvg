// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "actsvg/core/style.hpp"
#include "actsvg/proto/grid.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {

namespace proto {

/** A proto portal class as a simple translation layer
 * from a portal description
 *
 * @tparam point3_container a vertex description of surfaces
 **/

template <typename point3_container>
struct portal {
    /// A nested link type
    struct link {
        /// The position of the volume link
        typename point3_container::value_type _position;
        /// The direciton of the volume link
        typename point3_container::value_type _direction;
        /// The stroke style (ideally synchronized with volume)
        style::stroke _stroke;
        /// The span
        std::optional<std::array<scalar, 2>> _span = std::nullopt;
        /// Binning type
        std::string _binning = "";
        /// Auxiliary information as container map
        std::map<std::string, std::vector<std::string>> _aux_info = {};
    };

    /// Name of the surface
    std::string _name = "unnamed";

    /// Auxiliary information as container map
    std::map<std::string, std::vector<std::string>> _aux_info = {};

    /// The surface representation of this portal
    surface<point3_container> _surface;

    /// The list of volume links
    std::vector<link> _volume_links;
};

}  // namespace proto

}  // namespace actsvg