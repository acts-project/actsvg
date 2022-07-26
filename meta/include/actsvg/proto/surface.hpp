// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <optional>
#include <string>
#include <utility>
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

    enum type {
        e_annulus,
        e_cylinder,
        e_disc,
        e_polygon,
        e_rectangle,
        e_trapez
    };

    enum boolean { e_clipping, e_union, e_subtraction, e_none };

    /// Name of the surface
    std::string _name = "unnamed";

    /// Auxiliary information as container map
    std::map<std::string, std::vector<std::string>> _aux_info = {};

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
    svg::object _template_object;

    using point3_type = typename point3_container::value_type;

    /** Static constructor from a templat
     * @param t_ the template
     * @param o_ the tempalte object
     * @param name_ the new name
     **/
    static surface<point3_container> from_template(
        const surface<point3_container>& t_, const svg::object& o_,
        const std::string& name_) {
        surface<point3_container> s;
        s._name = name_;
        s._type = t_._type;
        s._aux_info = t_._aux_info;
        s._vertices = t_._vertices;
        s._measures = t_._measures;
        s._radii = t_._radii;
        s._opening = t_._opening;
        s._fill = t_._fill;
        s._stroke = t_._stroke;
        s._transform = t_._transform;
        s._template_object = o_;
        return s;
    }
};

}  // namespace proto

}  // namespace actsvg