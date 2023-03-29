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
#include "actsvg/styles/defaults.hpp"

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
svg::object grid(const std::string& id_, const proto::grid& g_) {

    svg::object g;

    if (g_._type == proto::grid::type::e_r_phi) {
        g = draw::tiled_polar_grid(id_, g_._edges_0, g_._edges_1, __g_fill,
                                   __g_stroke);
    } else {

        std::vector<scalar> edges_1 = g_._edges_1;
        if (g_._type == proto::grid::type::e_z_phi) {
            std::for_each(edges_1.begin(), edges_1.end(),
                          [&](auto& e) { e *= g_._reference_r; });
        }

        g = draw::tiled_cartesian_grid(id_, g_._edges_0, edges_1, __g_fill,
                                       __g_stroke);
    }

    return g;
}

}  // namespace display
}  // namespace actsvg
