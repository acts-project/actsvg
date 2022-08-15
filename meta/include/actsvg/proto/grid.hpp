// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include <string>

#include "actsvg/core/defs.hpp"

namespace actsvg {

namespace proto {

    /** A proto grid class to describe the grid setup
     * 
     * For convenience, it is forced to be 2-dimensional, for 1-dim
     * grid descriptions, only provide 2 eges
     * 
     */
    struct grid {

        /** Type of grid, enum defintion */ 
        enum type { e_x_y = 0, e_r_phi = 1, e_z_phi };

        /// Name the type
        type _type = e_r_phi;

        /// The edges in the two given directions, loc0
        std::vector<scalar> _edges_0 = {};
        
        /// The edges in the two given directions, loc1
        std::vector<scalar> _edges_1 = {};
        
    };

}  // namespace proto

}  // namespace actsvg

