// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core/generators.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "actsvg/core/defs.hpp"

namespace actsvg::generators {
std::vector<scalar> phi_values(scalar start_phi, scalar end_phi,
                               unsigned int lseg) {
    std::vector<scalar> values;
    values.reserve(lseg + 1);
    scalar step_phi = (end_phi - start_phi) / lseg;
    for (unsigned int istep = 0; istep <= lseg; ++istep) {
        values.push_back(start_phi + istep * step_phi);
    }
    return values;
}

std::vector<point2> sector_contour(scalar inner_r, scalar outer_r,
                                   scalar start_phi, scalar end_phi,
                                   unsigned int lseg) {

    // Re-bound phi
    if (start_phi > 0. and end_phi < 0.) {
        end_phi += 2_scalar * pi;
    }

    auto inner_phi = phi_values(end_phi, start_phi, lseg);
    auto outer_phi = phi_values(start_phi, end_phi, lseg);

    std::vector<point2> sector_vertices;
    sector_vertices.reserve(inner_phi.size() + outer_phi.size());
    for (auto iphi : inner_phi) {
        sector_vertices.push_back(
            {inner_r * std::cos(iphi), inner_r * std::sin(iphi)});
    }
    for (auto ophi : outer_phi) {
        sector_vertices.push_back(
            {outer_r * std::cos(ophi), outer_r * std::sin(ophi)});
    }
    return sector_vertices;
}
}  // namespace actsvg::generators
