// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "actsvg/core.hpp"

namespace actsvg {

namespace display {

/// Find inner outer radius at edges in STRIP PC
///
/// @note have a look at Acts/Surfaces/AnnulusBounds.hpp
/// for more information
///
static inline point2 annulusCircleIx(scalar O_x, scalar O_y, scalar r,
                                     scalar phi) {
    //                      _____________________________________________
    //                     /      2  2                    2    2  2    2
    //     O_x + O_y*m - \/  - O_x *m  + 2*O_x*O_y*m - O_y  + m *r  + r
    // x =
    // --------------------------------------------------------------
    //                                  2
    //                                 m  + 1
    //
    // y = m*x
    //
    scalar m = std::tan(phi);
    point2 dir = {std::cos(phi), std::sin(phi)};
    scalar x1 = static_cast<scalar>(
        (O_x + O_y * m -
         std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) + 2 * O_x * O_y * m -
                   std::pow(O_y, 2) + std::pow(m, 2) * std::pow(r, 2) +
                   std::pow(r, 2))) /
        (std::pow(m, 2) + 1));
    scalar x2 = static_cast<scalar>(
        (O_x + O_y * m +
         std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) + 2 * O_x * O_y * m -
                   std::pow(O_y, 2) + std::pow(m, 2) * std::pow(r, 2) +
                   std::pow(r, 2))) /
        (std::pow(m, 2) + 1));

    point2 v1 = {x1, m * x1};
    if (v1[0] * dir[0] + v1[1] * dir[1] > 0) {
        return v1;
    }
    return {x2, m * x2};
};

}  // namespace display
}  // namespace actsvg
