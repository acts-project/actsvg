// This file is part of the actsvg packge.
//
// Copyright (C) 2024 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <vector>

#include "actsvg/core/defs.hpp"
#include "actsvg/proto/track.hpp"

namespace actsvg {
namespace test {

/** Generate helix points along a helix with constant magnetic field
 * along a z axis
 *
 * @tparam point3_type The point type
 *
 * @param[in] origin_ The particle origin
 * @param[in] direction_ The particle direction
 * @param[in] radius_ The radius of the helix
 * @param[in] o_ The orientation of the helix (pos, neg)
 * @param[in] phi_range_ The range of the phi angle (phi==0 is the origin)
 * @param[in] step_size_ The step size
 *
 * @return The helix points
 */
template <typename point3_type>
inline static proto::trajectory<point3_type> generate_helix(
    const point3_type& origin_, const point3_type& direction_, scalar radius_,
    scalar o_, std::array<scalar, 2u> phi_range_, std::size_t n_points_) {
    // The helix points to be generated
    proto::trajectory<point3_type> trj;
    trj._origin = origin_;
    trj._direction = direction_;

    // Create the center of the helix
    scalar t_l = std::sqrt(direction_[0] * direction_[0] +
                           direction_[1] * direction_[1]);
    point2 t_dir = {direction_[0] / t_l, direction_[1] / t_l};
    point2 n_dir = {std::copysign(t_dir[1], o_), -std::copysign(t_dir[0], o_)};

    point3_type center = origin_;
    center[0] -= radius_ * n_dir[0];
    center[1] -= radius_ * n_dir[1];

    // Calculate phi0
    scalar phi0 = std::atan2(origin_[1] - center[1], origin_[0] - center[0]);

    // The number of steps
    scalar phi_step = (phi_range_[1] - phi_range_[0]) / n_points_;

    for (std::size_t i = 0; i < n_points_; ++i) {
        scalar phi = phi_range_[0] +
                     static_cast<scalar>(std::copysign(i, o_)) * phi_step;
        scalar x = center[0] + radius_ * std::cos(phi0 + phi);
        scalar y = center[1] + radius_ * std::sin(phi0 + phi);
        scalar z = center[2] + direction_[2] * phi;
        // Create the current position
        point3_type p{x, y, z};
        // Create the direction
        point3_type d{-std::sin(phi + phi0), std::cos(phi + phi0),
                      direction_[2]};
        trj._path.push_back({p, d});
    }
    // Return the trajectory
    return trj;
}

}  // namespace test
}  // namespace actsvg
