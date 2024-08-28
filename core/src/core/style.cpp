// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core/style.hpp"

#include <array>
#include <sstream>
#include <string>
#include <vector>

#include "actsvg/core/defs.hpp"
#include "actsvg/core/utils.hpp"

namespace actsvg {

namespace style {

std::string rgb_attr(const rgb &rgb_) {
    if (rgb_ == rgb{-1, -1, -1}) {
        return "none";
    }

    return std::string("rgb(") + std::to_string(rgb_[0]) + __c +
           std::to_string(rgb_[1]) + __c + std::to_string(rgb_[2]) +
           std::string(")");
}

fill::fill(const color &fc_) : _fc(fc_) {}

fill::fill(bool s_) : _sterile(s_) {
    _fc = color();
    _fc._opacity = 0.;
}

rgb gradient::rgb_from_scale(scalar s_) const {
    scalar s_reg =
        s_ < 0._scalar ? 0._scalar : (s_ > 1._scalar ? 1._scalar : s_);
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

stroke::stroke(const color &c_, scalar w_, const std::vector<int> &d_)
    : _sc(c_), _width(w_), _dasharray(d_) {}

stroke::stroke(bool s_) : _sterile(s_) {}

bool transform::is_identity() const {
    return _sterile or (_tr[0] == 0. and _tr[1] == 0. and _rot[0] == 0. and
                        _rot[1] == 0. and _rot[2] == 0. and _skew[0] == 0. and
                        _skew[1] == 0. and _scale[0] == 1. and _scale[1] == 1.);
}

std::string transform::attr() const {
    bool translate = (_tr[0] != 0. or _tr[1] != 0.);
    bool rotate = (_rot[0] != 0.);
    bool scale = (_scale[0] != 1. or _scale[1] != 1.);
    bool skew = (_skew[0] != 0. or _skew[1] != 0.);
    std::stringstream tr_str;
    if (translate) {
        tr_str << "translate(" << utils::to_string(_scale[0] * _tr[0]) << __c
               << utils::to_string(-_scale[0] * _tr[1]) << ")";
        if (rotate or scale or skew) {
            tr_str << __blk;
        }
    }
    if (rotate) {
        tr_str << "rotate(" << utils::to_string(-_rot[0]) << __c
               << utils::to_string(_rot[1]) << __c << utils::to_string(-_rot[2])
               << ")";
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

void label::place(const std::array<scalar, 2u> &lhc_,
                  const std::array<scalar, 2u> &rhc_) {

    scalar x = 0.;
    scalar y = 0.;

    // First determine the y position
    if (_vertical == vertical::top) {
        y = rhc_[1] + 0.6_scalar * _font._size;
    } else if (_vertical == vertical::bottom) {
        y = lhc_[1] - 1.1_scalar * _font._size;
    } else if (_vertical == vertical::center) {
        y = 0.5_scalar * (lhc_[1] + rhc_[1] - _font._size);
    }

    if (_horizontal == horizontal::left) {
        x = lhc_[0];
        if (_vertical == vertical::center) {
            x -= 0.64_scalar * _font._size * _text.size();
        }
    } else if (_horizontal == horizontal::right) {
        x = rhc_[0] - 0.6_scalar * _font._size * _text.size();
        if (_vertical == vertical::center) {
            x += 0.64_scalar * _font._size * _text.size();
        }
    } else if (_horizontal == horizontal::center) {
        x = 0.5_scalar *
            (lhc_[0] + rhc_[0] - 0.6_scalar * _font._size * _text.size());
    }

    _position = {x, y};
}
}  // namespace style
}  // namespace actsvg
