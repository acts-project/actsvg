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
#include "actsvg/proto/material.hpp"
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

    using point3 = typename point3_container::value_type;

    /// A transform3 structure before placement
    struct transform3 {
        // original translation
        point3 _translation = {0., 0., 0.};
        // original rotation
        std::array<point3, 3u> _rotation = {
            point3{1., 0., 0.}, point3{0., 1., 0.}, point3{0., 0., 1.}};
        // Rotate a point to the global frame
        point3 rotate(const point3& p_) const {

            point3 ret{0.f, 0.f, 0.f};

            ret[0] += _rotation[0][0] * p_[0];
            ret[1] += _rotation[1][0] * p_[0];
            ret[2] += _rotation[2][0] * p_[0];

            ret[0] += _rotation[0][1] * p_[1];
            ret[1] += _rotation[1][1] * p_[1];
            ret[2] += _rotation[2][1] * p_[1];

            ret[0] += _rotation[0][2] * p_[2];
            ret[1] += _rotation[1][2] * p_[2];
            ret[2] += _rotation[2][2] * p_[2];

            return ret;
        }

        /// Apply the translation and rotation
        point3 point_to_global(const point3& p_) const {
            // Apply the rotation
            point3 ret = rotate(p_);
            ret[0] += _translation[0];
            ret[1] += _translation[1];
            ret[2] += _translation[2];
            return ret;
        }

        // The identity tranform
        static transform3 identity() { return transform3{}; }
    };

    enum class type {
        e_annulus,
        e_cylinder,
        e_diamond,
        e_disc,
        e_polygon,
        e_rectangle,
        e_straw,
        e_trapez
    };

    enum class sf_type { e_portal, e_sensitive, e_passive };

    enum class boolean { e_clipping, e_union, e_subtraction, e_none };

    /// Name of the surface
    std::string _name = "unnamed_surface";

    /// Auxiliary information as container map
    std::map<std::string, std::vector<std::string>> _aux_info = {};

    /// The contained vertices - for polygon surfaces
    point3_container _vertices = {};

    /// This is the optional surface transform - to be applied at the
    /// vertices before the view is created
    std::optional<transform3> _surface_transform = std::nullopt;

    /// Dedicated regular disc/cylinder descriptions
    /// - if this is not applicable the _vertices view needs to be chosen
    std::array<scalar, 2> _radii = {0., 0.};
    std::array<scalar, 2> _opening = {-pi, pi};
    std::array<scalar, 2> _zparameters = {0., 0.};

    /// Boolean surfaces
    std::vector<surface<point3_container>> _boolean_surface = {};
    boolean _boolean_operation = boolean::e_none;

    /// Fill and stroke
    style::fill _fill = defaults::__s_fill;
    style::stroke _stroke = defaults::__s_stroke;
    style::transform _transform = defaults::__t_identity;

    /// Type of the surfaces
    type _type = type::e_trapez;
    sf_type _sf_type = sf_type::e_sensitive;

    /// And their measures
    std::vector<scalar> _measures = {};

    /// Decoration
    std::map<std::string, svg::object> _decorations;

    /// The surface material
    surface_material _material = {};

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
        s._sf_type = t_._sf_type;
        s._aux_info = t_._aux_info;
        s._vertices = t_._vertices;
        s._measures = t_._measures;
        s._radii = t_._radii;
        s._opening = t_._opening;
        s._fill = t_._fill;
        s._stroke = t_._stroke;
        s._transform = t_._transform;
        s._template_object = o_;
        s._decorations = t_._decorations;
        return s;
    }
};

}  // namespace proto

}  // namespace actsvg
