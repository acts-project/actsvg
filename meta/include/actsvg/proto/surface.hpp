// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>

#include "actsvg/core/style.hpp"
#include "actsvg/styles/defaults.hpp"
#include "grid.hpp"

namespace actsvg {

namespace proto {

/** A proto surface class as a simple translation layer
 * from a surface description
 *
 * @tparam point3_container a vertex description of surfaces
 **/

template <typename point3_container>
struct surface {

    enum type { e_annulus, e_cylinder, e_disc, e_rectangle, e_trapez };

    /// Name of the volume
    std::string _name = "unnamed";

    /// Auxiliary information
    std::vector<std::string> _info = {};

    /// The contained surfaces
    point3_container _vertices = {};
    style::fill _fill = defaults::__s_fill;
    style::stroke _stroke = defaults::__s_stroke;

    /// Type of the surfaces
    type _type = e_trapez;
    /// And their measures
    std::vector<scalar> _measures = {};
};

}  // namespace proto

}  // namespace actsvg