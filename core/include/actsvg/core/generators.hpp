// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "defs.hpp"

namespace actsvg {

namespace generators {

/** Generate phi values, respecting phi boundaries
 *
 * @param start_phi is the start for the arc generation
 * @param end_phi is the end of the arc generation
 * @param lseg is the number of segments used to gnerate the arc
 *
 * @return a vector of phi values for the arc
 */
std::vector<scalar> phi_values(scalar start_phi, scalar end_phi,
                               unsigned int lseg);

/** Generate a contour for a sector
 *
 * @param inner_r is the inner radius
 * @param outer_r is the outer radius
 * @param start_phi is the start for the arc generation
 * @param end_phi is the end of the arc generation
 * @param lseg is the number of segments to approximate the arc
 *
 * @note this contour generation DOES NOT perform the y flip
 * @note start/end are respected in terms of orientation
 *
 **/
std::vector<point2> sector_contour(scalar inner_r, scalar outer_r,
                                   scalar start_phi, scalar end_phi,
                                   unsigned int lseg = 16);
}  // namespace generators

}  // namespace actsvg
