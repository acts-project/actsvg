// This file is part of the actsvg packge.
//
// Copyright (C) 2022 Andreas Salzburger <asalzburger@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <vector>

#include "actsvg/core/defs.hpp"

namespace actsvg {

/// This file is generated from the Open Data detector and
/// depicts one Pixel endcap layer
namespace data {

using trapezoid = std::vector<std::array<scalar, 3u>>;

std::vector<trapezoid> generate_endcap_modules();
}  // namespace data

}  // namespace actsvg
