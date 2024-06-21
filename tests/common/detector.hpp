// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <vector>

#include "actsvg/core/defs.hpp"
#include "actsvg/data/odd_pixel_barrel.hpp"
#include "actsvg/data/odd_pixel_endcap.hpp"
namespace actsvg {

namespace test {

using rectangle = std::array<std::array<scalar, 3u>, 4u>;

inline static std::vector<rectangle> generate_barrel_modules() {
    std::vector<rectangle> modules;
    size_t number_of_modules = data::odd_pixel_barrel.size() / 4u;
    modules.reserve(number_of_modules);
    for (size_t im = 0; im < number_of_modules; ++im) {
        modules.push_back({data::odd_pixel_barrel[4 * im],
                           data::odd_pixel_barrel[4 * im + 1],
                           data::odd_pixel_barrel[4 * im + 2],
                           data::odd_pixel_barrel[4 * im + 3]});
    }
    return modules;
}

}  // namespace test
}  // namespace actsvg
