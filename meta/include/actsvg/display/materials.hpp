// This file is part of the actsvg packge.
//
// Copyright (C) 2024 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

#include "actsvg/core.hpp"
#include "actsvg/display/grids.hpp"
#include "actsvg/proto/grid.hpp"
#include "actsvg/proto/material.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {

using namespace defaults;

namespace display {

/** Draw a grid
 *
 * @param id_ the identification of this grid
 * @param m_ the surface material
 * @param info_pos_ the position of the info text
 *
 * @return an svg object representing the grid
 */
svg::object surface_material(const std::string& id_,
                             const proto::surface_material& m_);
}  // namespace display
}  // namespace actsvg
