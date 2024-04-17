// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "actsvg/core/defs.hpp"
#include "actsvg/proto/cluster.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {
namespace proto {

template <typename point3_type>
struct trajectory {

    using point3 = point3_type;

    /// The origin of the trajectory
    point3_type _origin;

    /// The direction of the trajectory
    point3_type _direction;

    /// The path of the trajectory, first positions, second optional directions
    std::vector<std::pair<point3_type, std::optional<point3_type>>> _path;

    /// The origin style
    scalar _origin_size = 0.;
    style::stroke _origin_stroke;
    style::fill _origin_fill;

    /// The trajectory style
    style::marker _path_arrow;
    style::stroke _path_stroke;
};

template <typename point3_type>
struct seed {

    /// The space pionts
    std::vector<point3_type> _space_points;

    /// The trajectory
    std::optional<trajectory<point3_type>> _trajectory;

    /// The space point style
    scalar _space_point_size = 5.;
    style::fill _space_point_fill;
    style::stroke _space_point_stroke;
};

template <typename surface_type>
struct track_state {
    /// The surface
    surface_type _surface;

    /// The cluster
    std::variant<cluster<1>, cluster<2>> _cluster;
};

}  // namespace proto

}  // namespace actsvg
