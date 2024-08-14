// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core/svg.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/core/defs.hpp"
#include "actsvg/core/style.hpp"

namespace actsvg::svg {
bool object::is_defined() const {
    return (not _tag.empty());
}

bool object::is_empty_group() const {
    return (_tag == "g" and _sub_objects.empty());
}
void object::add_object(const svg::object &o_) {
    if (o_._active and not o_.is_empty_group()) {
        // Add the object
        _sub_objects.push_back(o_);
        // Collect eventual definitions
        _definitions.insert(_definitions.end(), o_._definitions.begin(),
                            o_._definitions.end());
        // Re-measure, x/y/r/phi ranges include transforms already
        _x_range = {std::min(_x_range[0], o_._x_range[0]),
                    std::max(_x_range[1], o_._x_range[1])};
        _y_range = {std::min(_y_range[0], o_._y_range[0]),
                    std::max(_y_range[1], o_._y_range[1])};
        _r_range = {std::min(_r_range[0], o_._r_range[0]),
                    std::max(_r_range[1], o_._r_range[1])};
        _phi_range = {std::min(_phi_range[0], o_._phi_range[0]),
                      std::max(_phi_range[1], o_._phi_range[1])};
    }
}

std::optional<svg::object> object::find_object(const std::string &id_) const {

    auto found_object = std::find_if(
        _sub_objects.begin(), _sub_objects.end(),
        [&](const svg::object &candidate_) { return (candidate_._id == id_); });

    if (found_object != _sub_objects.end()) {
        return (*found_object);
    }
    return std::nullopt;
}

/** Stream operator with @param os_ the output stream and @param o_ the streamed
 * object */
std::ostream &operator<<(std::ostream &os_, const object &o_) {

    if (o_.is_empty_group()) {
        return os_;
    }

    // We make a temporary copy for writing, this way we can
    // write the same one with different attributes sets
    object o_copy = o_;

    // Write the file
    os_ << __l << o_copy._tag;
    if (not o_copy._id.empty()) {
        os_ << __blk << "id=\"" << o_copy._id << "\"";
    }

    // Attach the styles: fill, stroke,
    if (not o_._sterile) {
        o_._fill.attach_attributes(o_copy);
        o_._stroke.attach_attributes(o_copy);
    }

    // Attach the transform
    if (not o_._transform.is_identity()) {
        o_._transform.attach_attributes(o_copy);
    }

    // The attribute map
    for (auto [key, value] : o_copy._attribute_map) {
        os_ << __blk << key << "=\"" << value << "\"";
    }

    // This is done return
    if (o_copy._sub_objects.empty() and o_copy._field.empty()) {
        os_ << __er;
        return os_;
    }

    os_ << __r;
    for (const auto &a : o_copy._sub_objects) {
        os_ << a;
    }
    if (not o_copy._field.empty()) {
        if (o_copy._field.size() == 1) {
            os_ << o_copy._field[0];
        } else {
            for (const auto &fl : o_copy._field) {
                os_ << "<tspan x=\"";
                os_ << o_copy._x_range[0] << "\"";
                os_ << " dy=\"" << o_copy._field_span << "\">";
                os_ << fl;
                os_ << "</tspan>" << __nl;
            }
        }
    }
    // Close-out
    os_ << __el << o_copy._tag << __r;
    return os_;
}

void file::add_object(const svg::object &o_) {
    // Add the object
    if (o_._active) {
        _objects.push_back(o_);
    }
}

void file::add_objects(const std::vector<svg::object> &os_) {
    // Add the objects one by one
    for (const auto &o_ : os_)
        if (o_._active) {
            _objects.push_back(o_);
        }
}

std::ostream &operator<<(std::ostream &os_, const file &f_) {
    // Do the viewBox adjustment
    std::array<scalar, 2> x_range = {std::numeric_limits<scalar>::max(),
                                     std::numeric_limits<scalar>::min()};
    std::array<scalar, 2> y_range = {std::numeric_limits<scalar>::max(),
                                     std::numeric_limits<scalar>::min()};

    std::array<scalar, 4> viewBox = {-800, -800, 1600, 1600};

    std::map<std::string, svg::object> definitions;

    for (const auto &o : f_._objects) {
        x_range[0] = std::min(x_range[0], o._x_range[0]);
        x_range[1] = std::max(x_range[1], o._x_range[1]);
        y_range[0] = std::min(y_range[0], o._y_range[0]);
        y_range[1] = std::max(y_range[1], o._y_range[1]);
        for (const auto &d : o._definitions) {
            definitions[d._id] = d;
        }
    }
    // Enlarge the view box by 10 percent
    viewBox[2] = 1.2_scalar * (x_range[1] - x_range[0]);
    viewBox[3] = 1.2_scalar * (y_range[1] - y_range[0]);

    // Include a fixed size border
    viewBox[0] = (x_range[0] - 0.1_scalar * viewBox[2]) - f_._border;
    viewBox[1] = (y_range[0] - 0.1_scalar * viewBox[3]) - f_._border;
    viewBox[2] += f_._border;
    viewBox[3] += f_._border;

    if (f_._add_html) {
        os_ << f_._html_head;
    }
    os_ << f_._svg_head;
    os_ << " width=\"" << f_._width << "\" height=\"" << f_._height << "\"";
    os_ << " viewBox=\"" << viewBox[0] << " " << viewBox[1] << " " << viewBox[2]
        << " " << viewBox[3] << "\"";
    os_ << f_._svg_def_end;
    // Write the definitions first
    if (not definitions.empty()) {
        os_ << "<defs>";
        for (auto [key, value] : definitions) {
            os_ << value;
        }
        os_ << "</defs>";
    }

    // Now write the objects
    for (auto &o : f_._objects) {
        os_ << o;
    }
    os_ << f_._svg_tail;
    if (f_._add_html) {
        os_ << f_._html_tail;
    }
    return os_;
}

}  // namespace actsvg::svg
