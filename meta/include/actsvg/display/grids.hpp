// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

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

/** Draw a grid
 *
 * @param id_ the identification of this grid
 * @param g_ the grid to be drawn
 *
 * @return an svg object representing the grid
 */
svg::object grid(const std::string& id_, const proto::grid& g_);
}  // namespace display
}  // namespace actsvg
