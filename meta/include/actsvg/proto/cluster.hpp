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

#include "actsvg/core/style.hpp"
#include "actsvg/styles/defaults.hpp"
#include "grid.hpp"

namespace actsvg {

namespace proto {

template <size_t DIM = 1>
struct channel {
    /// Channel identification
    std::array<unsigned int, DIM> _cid;
    /// The data
    scalar _data;
};

/** A proto surface class as a simple translation layer
 * from a surface description
 *
 * @tparam point3_container a vertex description of surfaces
 **/
template <size_t DIM = 1>
struct cluster {

    enum type { e_cartesian, e_fan, e_polar, e_drift };

    enum coordinate { e_x, e_y, e_r, e_phi };

    /// The type of the cluster
    type _type = e_cartesian;

    /// Measurements coordinates
    std::array<coordinate, DIM> _coords;

    /// The list of channels
    std::vector<channel<DIM>> _channels;

    /// Reconstructed variables
    std::array<scalar, DIM> _measurement;
    std::array<scalar, DIM> _variance;
    scalar _correlation = 0.;

    /// Truth
    std::array<scalar, DIM> _truth;
    bool _mc = true;
};

}  // namespace proto

}  // namespace actsvg