// This file is part of the actsvg packge.
//
// Copyright (C) 2024 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>

#include "actsvg/core/defs.hpp"
#include "actsvg/proto/grid.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {

namespace proto {

struct material_slab {
    // x0, l0, a, z, rho, thickness
    std::array<scalar, 6u> _properties = {0., 0., 0., 0., 0., 0.};
};

using matrial_matrix = std::vector<std::vector<material_slab> >;

/** A helper function to get the min max ranges from material properties
 *
 * @param mm_ the material slab
 *
 * @return the boundaries in t/X0, t/L0, a*z/rho
 */
std::array<std::array<scalar, 2u>, 3u> material_ranges(
    const matrial_matrix& mm_);

/** A proto surface material class that represents material properties
 *  in a binned and non-binned manner
 **/
struct surface_material {

    /// The material matrix
    std::vector<std::vector<material_slab> > _material_matrix = {};

    /// The bin ids in a grid
    proto::grid _grid = {};

    // The material ranges
    std::array<std::array<scalar, 2u>, 3u> _material_ranges = {};

    /// The gradient
    style::gradient _gradient = defaults::__rgb_gradient;

    /// The gradient pos
    point2 _gradient_pos = {0, 0};

    /// Gradient bos size
    std::array<scalar, 2u> _gradient_box = {0., 0.};

    /// Gradient stroke
    style::stroke _gradient_stroke = style::stroke{};

    /// Gradient font
    style::font _gradient_font = style::font{};

    /// Gradient label
    style::label _gradient_label = style::label{};

    /// The info position
    point2 _info_pos = {0, 0};

    /// The info font
    style::font _info_font = style::font{};

    /// Evaluate the material ranges
    void evaluate_material_ranges() {
        _material_ranges = material_ranges(_material_matrix);
    }
};

}  // namespace proto
}  // namespace actsvg
