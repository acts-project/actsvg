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

    // Expose the container type
    using container_type = point3_container;

    /// A nested link type
    struct link {

        // Expose the point3 type
        using point3_type = typename container_type::value_type;

        /// The start of the volume link
        point3_type _start;
        /// The end of the volume link
        point3_type _end;

        /// The stroke style (ideally synchronized with volume)
        style::marker _start_marker = style::marker({});
        style::marker _end_marker = style::marker({"<<"});
        style::stroke _stroke;

        /// @brief  The link index of the volume
        unsigned int _link_index = 0u;

        /// The span
        std::optional<std::array<scalar, 2>> _span = std::nullopt;
        /// Binning type
        std::string _binning = "";
        /// Auxiliary information as container map
        std::map<std::string, std::vector<std::string>> _aux_info = {};
    };

    /// Name of the surface
    std::string _name = "unnamed_portal";

    /// Auxiliary information as container map
    std::map<std::string, std::vector<std::string>> _aux_info = {};

    /// The surface representation of this portal
    surface<point3_container> _surface;

    /// The list of volume links
    std::vector<link> _volume_links;
};

}  // namespace proto

}  // namespace actsvg