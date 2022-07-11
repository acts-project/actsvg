// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <sstream>
#include <string>

#include "defs.hpp"
#include "utils.hpp"

namespace actsvg {

namespace style {

using rgb = std::array<int, 3>;

/** Helper method to convert and @return a color into an attribute
 *
 * @param rgb_ is the r,g,b color to be represented
 *
 */
static std::string rgb_attr(const rgb &rgb_) {
    if (rgb_ == rgb{-1, -1, -1}) {
        return "none";
    }

    return std::string("rgb(") + std::to_string(rgb_[0]) + __c +
           std::to_string(rgb_[1]) + __c + std::to_string(rgb_[2]) +
           std::string(")");
}

/// Color specification
struct color {
    /// The color
    rgb _rgb = {255, 255, 255};
    /// The opacity
    scalar _opacity = 1.;
    /// The highlight mode it is assumed on/off
    std::vector<std::string> _highlight = {};
    rgb _hl_rgb = {255, 0, 0};
};

/// Fill type specification
struct fill {

    /// The fill color
    color _fc = color{{0, 0, 0}};
    bool _sterile = false;

    /// A constructor from @param fc_ color
    fill(const color &fc_) : _fc(fc_) {}

    /// A constructor for empty
    fill(bool s_ = false) : _sterile(s_) {
        _fc = color();
        _fc._opacity = 0.;
    }

    /** Attach this fill attribute to an object
     *
     * @tparam object_type the type of the object
     *
     * @param o_ [in,out] the object in question
     **/
    template <typename object_type>
    void attach_attributes(object_type &o_) const {

        if (not _sterile) {
            o_._attribute_map["fill"] = rgb_attr(_fc._rgb);
            o_._attribute_map["fill-opacity"] = std::to_string(_fc._opacity);
        }

        if (_fc._highlight.size() == 2u) {
            object_type on_off;
            on_off._tag = "set";
            on_off._attribute_map["attributeName"] = "fill";
            on_off._attribute_map["begin"] = _fc._highlight[0];
            on_off._attribute_map["end"] = _fc._highlight[1];
            on_off._attribute_map["to"] = rgb_attr(_fc._hl_rgb);
            o_.add_object(on_off);
        }
    }
};

/// Stroke type speficiation
struct stroke {

    /// The stroke color
    color _sc{{0, 0, 0}};
    /// Width definition
    scalar _width = 0.5;
    /// Dashing definition
    std::vector<int> _dasharray = {};
    /// Nothing is written out
    bool _sterile = false;

    stroke(const color &c_, scalar w_ = 0.5, const std::vector<int> &d_ = {})
        : _sc(c_), _width(w_), _dasharray(d_) {}

    stroke(bool s_ = false) : _sterile(s_) {}

    /** Attach this fill attribute to an object
     *
     * @tparam object_type the type of the object
     *
     * @param o_ [in,out] the object in question
     **/
    template <typename object_type>
    void attach_attributes(object_type &o_) const {

        if (not _sterile) {
            o_._attribute_map["stroke"] = rgb_attr(_sc._rgb);
            o_._attribute_map["stroke-opacity"] = std::to_string(_sc._opacity);
            o_._attribute_map["stroke-width"] = std::to_string(_width);
            if (not _dasharray.empty()) {
                std::string da_str;
                for (auto [i, d] : utils::enumerate(_dasharray)) {
                    da_str += std::to_string(d);
                    if (i + 1 < _dasharray.size()) {
                        da_str += __blk;
                    }
                }
                o_._attribute_map["stroke-dasharray"] = da_str;
            }
        }
    }
};

/// Font sytle specification
struct font {
    /// The font color
    color _fc{{0, 0, 0}};

    std::string _family = "Arial";

    unsigned int _size = 12;

    scalar _line_spacing = 1.4;

    std::string _style = "";

    /** Attach this fill attribute to an object
     *
     * @tparam object_type the type of the object
     *
     * @param o_ [in,out] the object in question
     **/
    template <typename object_type>
    void attach_attributes(object_type &o_) const {
        o_._attribute_map["fill"] = rgb_attr(_fc._rgb);
        o_._attribute_map["font-size"] = std::to_string(_size);
        if (not _family.empty()) {
            o_._attribute_map["font-family"] = _family;
        }
        if (not _style.empty()) {
            o_._attribute_map["font-style"] = _style;
        }
    }
};

/// The transform struct
struct transform {

    std::array<scalar, 3> _tr = {0., 0.};
    std::array<scalar, 3> _rot = {0., 0., 0.};
    std::array<scalar, 2> _skew = {0., 0.};
    std::array<scalar, 2> _scale = {1., 1.};

    bool _sterile = false;

    /** Attrribute conversion
     *
     * @note that the scale is directly applied on the objects,
     * in order to control the viewBox boundaries
     **/
    std::string attr() const {
        bool translate = (_tr[0] != 0. or _tr[1] != 0.);
        bool rotate = (_rot[0] != 0.);
        bool scale = (_scale[0] != 1. or _scale[1] != 1.);
        bool skew = (_skew[0] != 0. or _skew[1] != 0.);
        std::stringstream tr_str;
        if (translate) {
            tr_str << "translate(" << _tr[0] << __c << -_tr[1] << ")";
            if (rotate or scale or skew) {
                tr_str << __blk;
            }
        }
        if (rotate) {
            tr_str << "rotate(" << -_rot[0] << __c << _rot[1] << __c << -_rot[2]
                   << ")";
            if (scale or skew) {
                tr_str << __blk;
            }
        }
        if (skew) {
            tr_str << "skewX(" << _skew[0] << ") skewY(" << _skew[1] << ")";
        }
        return tr_str.str();
    }

    /** Attach this fill attribute to an object
     *
     * @tparam object_type the type of the object
     *
     * @param o_ [in,out] the object in question
     *
     * @note is applies the transform to the x range
     **/
    template <typename object_type>
    void attach_attributes(object_type &o_) const {
        auto transform_attribute = attr();
        if (not transform_attribute.empty()) {
            o_._attribute_map["transform"] = transform_attribute;
            scalar tx = _tr[0];
            scalar ty = _tr[1];
            o_._x_range = {o_._x_range[0] + tx, o_._x_range[1] + tx};
            o_._y_range = {o_._y_range[0] + ty, o_._y_range[1] + ty};
        }
    }
};

/// The marker struct
///
/// Allowed types are:
/// - none, <<, <, <|, |<, |<<, o
struct marker {
    std::string _type = "none";

    scalar _size = 4.;

    fill _fill = fill{{{0, 0, 0}}};

    stroke _stroke = stroke();
};

// The axis marker types
template <unsigned int kDIM>
using axis_markers = std::array<std::array<marker, 2u>, kDIM>;

/// Some standard styles to be used as default

}  // namespace style

static style::marker __no_marker = style::marker();
static style::marker __standard_marker = style::marker{{"<<"}};
static std::array<style::marker, 2u> __standard_axis_markers = {
    __no_marker, __standard_marker};

}  // namespace actsvg