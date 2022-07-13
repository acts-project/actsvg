// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <optional>
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

    enum type { e_annulus, e_cylinder, e_disc, e_polygon, e_rectangle, e_trapez };

    enum boolean { e_clipping, e_union, e_subtraction, e_none };

    /// Name of the surface
    std::string _name = "unnamed";

    /// Auxiliary information
    std::vector<std::string> _info = {};

    /// The contained vertices - for polygon surfaces
    point3_container _vertices = {};

    /// Dedicated disc descriptions, simplifies the set
    std::array<scalar, 2> _radii = {0., 0.};
    std::array<scalar, 2> _opening = {-M_PI, M_PI};

    /// Boolean surfaces
    std::vector<surface<point3_container>> _boolean_surface = {};
    boolean _boolean_operation = e_none;

    /// Fill and stroke
    style::fill _fill = defaults::__s_fill;
    style::stroke _stroke = defaults::__s_stroke;
    style::transform _transform = defaults::__t_identity;

    /// Type of the surfaces
    type _type = e_trapez;

    /// And their measures
    std::vector<scalar> _measures = {};

    /// A (potential) template for this surface
    svg::object _template;
};

}  // namespace proto

}  // namespace actsvg