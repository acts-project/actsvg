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

#include "grid.hpp"
#include "actsvg/core/style.hpp"

namespace actsvg {

namespace proto {

/** A proto volume class as a simple translation layer 
 * from a volume description
 *
 * @tparam point3_container a vertex description of surfaces
 **/

template <typename point3_container>
struct volume {

    /// Type enumeration 
    enum type { e_barrel = 0, e_endcap = 1, e_other = 2};

    unsigned int _index = 0u;

    /// Name of the volume
    std::string _name = "unnamed";

    type _type = e_barrel;

    /// Auxiliary information 
    std::vector<std::string> _info = {};

    /// The contained surfaces 
    std::vector<point3_container> _surfaces = {};
    style::fill _surface_fill;
    style::stroke _surface_stroke;
    /// The associated surface grid 
    grid _surface_grid;

    /// The portals 
    std::vector<point3_container> _portals = {};
    style::fill _portals_fill;
    style::stroke _portals_stroke;

};

}  // namespace proto

}  // namespace actsvg