// This file is part of the actsvg packge.
//
// Copyright (C) 2024 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/proto/material.hpp"

#include <array>

#include "actsvg/core/defs.hpp"
#include "actsvg/proto/grid.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {

namespace proto {
std::array<std::array<scalar, 2u>, 3u> material_ranges(
    const matrial_matrix& mm_) {
    std::array<std::array<scalar, 2u>, 3u> material_ranges;
    material_ranges[0u][0u] = std::numeric_limits<scalar>::max();
    material_ranges[0u][1u] = 0.;
    material_ranges[1u][0u] = std::numeric_limits<scalar>::max();
    material_ranges[1u][1u] = 0.;
    material_ranges[2u][0u] = std::numeric_limits<scalar>::max();
    material_ranges[2u][1u] = 0.;
    // Loop over the material matrix and fill min/max ranges
    for (const auto& row : mm_) {
        for (const auto& slab : row) {
            scalar t_X0 = slab._properties[5u] / slab._properties[0u];
            scalar t_L0 = slab._properties[5u] / slab._properties[1u];
            scalar a_z_rho = slab._properties[2] * slab._properties[3] /
                             (slab._properties[4] * slab._properties[4]);
            material_ranges[0u][0u] = std::min(material_ranges[0u][0u], t_X0);
            material_ranges[0u][1u] = std::max(material_ranges[0u][1u], t_X0);
            material_ranges[1u][0u] = std::min(material_ranges[1u][0u], t_L0);
            material_ranges[1u][1u] = std::max(material_ranges[1u][1u], t_L0);
            material_ranges[2u][0u] =
                std::min(material_ranges[2u][0u], a_z_rho);
            material_ranges[2u][1u] =
                std::max(material_ranges[2u][1u], a_z_rho);
        }
    }
    // Return the ranges
    return material_ranges;
}
}  // namespace proto
}  // namespace actsvg
