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
    add_objects({o_});
}

void file::add_objects(const std::vector<svg::object> &os_) {
    // Add the objects one by one
    for (const auto &o : os_) {
        if (o._active) {
            _objects.push_back(o);
            // Adjust range
            _x_range[0] = std::min(_x_range[0], o._x_range[0]);
            _x_range[1] = std::max(_x_range[1], o._x_range[1]);
            _y_range[0] = std::min(_y_range[0], o._y_range[0]);
            _y_range[1] = std::max(_y_range[1], o._y_range[1]);
        }
    }
}

std::array<scalar, 4> file::view_box() const {

    std::array<scalar, 4> vbox = {-800, -800, 1600, 1600};

    // Enlarge the view box by 10 percent
    vbox[2] = 1.2_scalar * (_x_range[1] - _x_range[0]);
    vbox[3] = 1.2_scalar * (_y_range[1] - _y_range[0]);

    // Include a fixed size border
    vbox[0] = (_x_range[0] - 0.1_scalar * vbox[2]) - _border;
    vbox[1] = (_y_range[0] - 0.1_scalar * vbox[3]) - _border;
    vbox[2] += _border;
    vbox[3] += _border;

    return vbox;
}

void file::set_view_box(const std::array<scalar, 4> &vbox_) {
    _view_box = vbox_;
}

std::ostream &operator<<(std::ostream &os_, const file &f_) {

    std::map<std::string, svg::object> definitions;
    for (const auto &o : f_._objects) {
        if (o._active) {
            for (const auto &d : o._definitions) {
                definitions[d._id] = d;
            }
        }
    }

    std::array<scalar, 4> vbox = f_._view_box.value_or(f_.view_box());

    if (f_._add_html) {
        os_ << f_._html_head;
    }
    os_ << f_._svg_head;
    os_ << " width=\"" << f_._width << "\" height=\"" << f_._height << "\"";
    os_ << " viewBox=\"" << vbox[0] << " " << vbox[1] << " " << vbox[2] << " "
        << vbox[3] << "\"";
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
