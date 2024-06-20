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
#include <vector>

#include "defs.hpp"
#include "utils.hpp"
#include "views.hpp"

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
            o_._attribute_map["fill-opacity"] = utils::to_string(_fc._opacity);
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

// Define a gradient definition
struct gradient {

    /// Define a gradient id
    std::string _id = "";
    /// The gradient direction (x1, y1, 2, y2)
    /// default is left to right
    std::array<scalar, 4> _direction = {0., 0., 1., 0.};
    /// The type of the gradient
    std::string _type = "linear";

    using stop = std::pair<scalar, color>;

    /// The gradient stops
    std::vector<stop> _stops = {};

    /// An optional label string (non-empty will trigger label rendering)
    std::string _label = "";

    /// Get a color from a scale parameter
    ///
    /// @param s_ the scale parameter for the lookup point
    rgb rgb_from_scale(scalar s_) const {
        scalar s_reg = s_ < 0. ? 0. : (s_ > 1. ? 1. : s_);
        // find our stops
        unsigned int is = 1u;
        for (; is <= _stops.size(); ++is) {
            if (_stops[is].first > s_reg) {
                break;
            }
        }
        // Bail out
        if (is >= _stops.size()) {
            return _stops.back().second._rgb;
        }

        // Interpolate between the two stops to get the new color
        scalar s0 = _stops[is - 1].first;
        scalar s1 = _stops[is].first;
        rgb c0 = _stops[is - 1].second._rgb;
        rgb c1 = _stops[is].second._rgb;
        scalar s_diff = s1 - s0;
        scalar s_rel = (s_reg - s0) / s_diff;

        rgb c;
        for (unsigned int ic = 0; ic < 3; ++ic) {
            c[ic] = static_cast<int>(c0[ic] + s_rel * (c1[ic] - c0[ic]));
        }
        return c;
    }
};

/// Stroke type speficiation
struct stroke {

    /// The stroke color
    color _sc{{0, 0, 0}};
    /// Width definition
    scalar _width = 0.5;
    /// hl with
    scalar _hl_width = 0.5;
    /// Dashing definition
    std::vector<int> _dasharray = {};
    /// Nothing is written out
    bool _sterile = false;

    /// @brief Contructor for stroke
    /// @param c_ the color of the stroke
    /// @param w_ the with of the stroke
    /// @param d_ the dashed harray of the stroke
    stroke(const color &c_, scalar w_ = 0.5, const std::vector<int> &d_ = {})
        : _sc(c_), _width(w_), _dasharray(d_) {}

    /// @brief Constructor for sterile stroke
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
            o_._attribute_map["stroke-opacity"] =
                utils::to_string(_sc._opacity);
            o_._attribute_map["stroke-width"] = utils::to_string(_width);
            // dashed array
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
            // Stroke color
            if (_sc._highlight.size() == 2u) {
                object_type on_off;
                on_off._tag = "set";
                on_off._attribute_map["attributeName"] = "stroke";
                on_off._attribute_map["begin"] = _sc._highlight[0];
                on_off._attribute_map["end"] = _sc._highlight[1];
                on_off._attribute_map["to"] = rgb_attr(_sc._hl_rgb);
                o_.add_object(on_off);
            }
            // Stroke width
            if (_width != _hl_width) {
                object_type on_off;
                on_off._tag = "set";
                on_off._attribute_map["attributeName"] = "stroke-width";
                on_off._attribute_map["begin"] = "mouseover";
                on_off._attribute_map["end"] = "mouseout";
                on_off._attribute_map["to"] = utils::to_string(_hl_width);
                o_.add_object(on_off);
            }
        }
    }
};

/// Font sytle specification
struct font {
    /// The font color
    color _fc{{0, 0, 0}};

    std::string _family = "Andale Mono";

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

/// The label struct
///
struct label {
    /// The horizontal alignment enum
    enum class horizontal { left, center, right };

    /// The vertical alignment enum
    enum class vertical { top, center, bottom };

    /// The label text
    std::string _text = "";

    /// The label horizontal alignment
    horizontal _horizontal = horizontal::left;

    /// The vertical horizontal alignment
    vertical _vertical = vertical::bottom;

    /// The font type of the label
    font _font = font{};

    /// The position
    std::array<scalar, 2u> _position = {0., 0.};

    /// Place the label at the position - assuming a bounding box
    ///
    /// @param lhc_ the left hand corner of the bounding box
    /// @param ruc_ the right upper corner of the bounding box
    void place(const std::array<scalar, 2u> &lhc_,
               const std::array<scalar, 2u> &rhc_) {

        scalar x = 0.;
        scalar y = 0.;

        // First determine the y position
        if (_vertical == vertical::top) {
            y = rhc_[1] + 0.6 * _font._size;
        } else if (_vertical == vertical::bottom) {
            y = lhc_[1] - 1.1 * _font._size;
        } else if (_vertical == vertical::center) {
            y = 0.5 * (lhc_[1] + rhc_[1] - _font._size);
        }

        if (_horizontal == horizontal::left) {
            x = lhc_[0];
            if (_vertical == vertical::center) {
                x -= 0.64 * _font._size * _text.size();
            }
        } else if (_horizontal == horizontal::right) {
            x = rhc_[0] - 0.6 * _font._size * _text.size();
            if (_vertical == vertical::center) {
                x += 0.64 * _font._size * _text.size();
            }
        } else if (_horizontal == horizontal::center) {
            x = 0.5 * (lhc_[0] + rhc_[0] - 0.6 * _font._size * _text.size());
        }

        _position = {x, y};
    }
};

/// The transform struct
///
/// This is a style transform, not an object transform,
/// i.e. it applies to where to draw and modify the object
/// while displaying
struct transform {

    std::array<scalar, 3> _tr = {0., 0., 0.};
    std::array<scalar, 3> _rot = {0., 0., 0.};
    std::array<scalar, 2> _skew = {0., 0.};
    std::array<scalar, 2> _scale = {1., 1.};

    bool _sterile = false;

    /** Test if it is a identity/sterile transform */
    bool is_identity() const {
        return _sterile or
               (_tr[0] == 0. and _tr[1] == 0. and _rot[0] == 0. and
                _rot[1] == 0. and _rot[2] == 0. and _skew[0] == 0. and
                _skew[1] == 0. and _scale[0] == 1. and _scale[1] == 1.);
    }

    /** Apply to a point
     *
     * @param p_ the point to be transformed
     * @param flip is the y-axis flipped (for viewing)
     *
     * @return the transformed point
     **/
    template <typename point_type>
    point_type apply(const point_type &p_, bool flip = true) const {
        // Apply the scale and the translation
        point_type p = p_;
        p[0] = _scale[0] * (p[0] + _tr[0]);
        p[1] = _scale[1] * (p[1] + _tr[1]);
        if (flip) {
            p[1] = -p[1];
        }
        return p;
    }

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
            tr_str << "translate(" << utils::to_string(_scale[0] * _tr[0])
                   << __c << utils::to_string(-_scale[0] * _tr[1]) << ")";
            if (rotate or scale or skew) {
                tr_str << __blk;
            }
        }
        if (rotate) {
            tr_str << "rotate(" << utils::to_string(-_rot[0]) << __c
                   << utils::to_string(_rot[1]) << __c
                   << utils::to_string(-_rot[2]) << ")";
            if (scale or skew) {
                tr_str << __blk;
            }
        }
        if (skew) {
            tr_str << "skewX(" << utils::to_string(_skew[0]) << ") skewY("
                   << utils::to_string(_skew[1]) << ")";
        }
        if (scale) {
            tr_str << "scale(" << utils::to_string(_scale[0]) << ","
                   << utils::to_string(_scale[1]) << ")";
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
    std::string _type = "";

    scalar _size = 4.;

    fill _fill = fill{{{0, 0, 0}}};

    stroke _stroke = stroke();
};

// The axis marker types
template <unsigned int kDIM>
using axis_markers = std::array<std::array<marker, 2u>, kDIM>;

/// Some standard styles to be used as default

}  // namespace style

static style::marker __no_marker = style::marker{};
static style::marker __standard_marker = style::marker{{"<<"}};
static std::array<style::marker, 2u> __standard_axis_markers = {
    __no_marker, __standard_marker};

}  // namespace actsvg