// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "defs.hpp"
#include "style.hpp"

namespace actsvg {

namespace svg {
/** An svg object
 *
 * It carries its attributes in a map, but can also have additional
 * animiations.
 *
 * To generate the final svg, it also records the barycenter and the view
 * range
 **/
struct object {

    /// SVG tag of the objec
    std::string _tag = "";

    /// Identification string
    std::string _id = "";

    /// In-field of the svg tag
    std::vector<std::string> _field = {};
    scalar _field_span = 12;

    /// Fill specifications
    style::fill _fill;

    /// Stroke specifications
    style::stroke _stroke;

    /// Transform of the object
    style::transform _transform;

    /// Attribute map for writing
    std::map<std::string, std::string> _attribute_map;

    /// Containes sub objects ( group, animation, etc.)
    std::vector<object> _sub_objects;

    /// Barycenter - Detector frame
    std::array<scalar, 2> _real_barycenter = {0., 0.};

    /// Range in x - Detector frame
    std::array<scalar, 2> _x_range = {std::numeric_limits<scalar>::max(),
                                      std::numeric_limits<scalar>::min()};
    /// Range in y - Detector frame
    std::array<scalar, 2> _y_range = {std::numeric_limits<scalar>::max(),
                                      std::numeric_limits<scalar>::min()};
    /// Range in r  - Detector frame
    std::array<scalar, 2> _r_range = {std::numeric_limits<scalar>::max(),
                                      std::numeric_limits<scalar>::min()};
    /// Range in phi - Detector frame
    std::array<scalar, 2> _phi_range = {std::numeric_limits<scalar>::max(),
                                        std::numeric_limits<scalar>::min()};
    /** Add a sub object and respect the min/max range
     *
     * @param o_ is the object in question
     **/
    void add_object(const svg::object &o_) {
        // Add the object
        _sub_objects.push_back(o_);
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

    /** Find an object by string,
     *
     * @param id_ is the identifier with which the sub object
     * is searched for
     *
     * @note std::optional is used and std::nullopt as @return if
     * the search was not successful
     */
    std::optional<svg::object> find_object(const std::string &id_) const {

        auto found_object =
            std::find_if(_sub_objects.begin(), _sub_objects.end(),
                         [&](const svg::object &candidate_) {
                             return (candidate_._id == id_);
                         });

        if (found_object != _sub_objects.end()) {
            return (*found_object);
        }
        return std::nullopt;
    }

    /** Add a sub object and respect the min/max range
     *
     * @param oc_ is the object container to be added
     **/
    template <typename object_container>
    void add_objects(const object_container &oc_) {
        // Add the object
        for (const auto &o : oc_) {
            add_object(o);
        }
    }

    friend std::ostream &operator<<(std::ostream &os_, const object &o_);
};

/** Stream operator with @param os_ the output stream and @param o_ the streamed
 * object */
inline std::ostream &operator<<(std::ostream &os_, const object &o_) {

    // We make a temporary copy for writing, this way we can
    // write the same one with different attributes sets
    object o_copy = o_;

    // Now write
    os_ << __l << o_copy._tag;
    if (not o_copy._id.empty()) {
        os_ << __blk << "id=\"" << o_copy._id << "\"";
    }

    // Attach the styles: fill, stroke, transform
    o_._fill.attach_attributes(o_copy);
    o_._stroke.attach_attributes(o_copy);
    o_._transform.attach_attributes(o_copy);

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
            for (const auto fl : o_copy._field) {
                os_ << "<tspan x=\"";
                os_ << o_copy._real_barycenter[0] << "\"";
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

/** An svg file scope
 *
 * This contains objects, connections and groups for final writing
 * It can be optionally augmented with an html bracket.
 *
 **/
struct file {

    /// Directive to add html header or not
    bool _add_html = false;

    /// Default header tail definitions
    std::string _html_head = "<html>\n<body>\n";
    std::string _svg_head = "<svg version=\"1.1\"";

    std::string _svg_def_end = " xmlns=\"http://www.w3.org/2000/svg\">\n";

    std::string _svg_tail = "</svg>\n";
    std::string _html_tail = "</body>\n</html>\n";

    scalar _width = 800;
    scalar _height = 800;
    scalar _border = 10;

    std::vector<object> _objects = {};

    /** Add an object and respect the min/max range
     *
     * @param o_ is the object in question
     **/
    void add_object(const svg::object &o_) {
        // Add the object
        _objects.push_back(o_);
    }

    friend std::ostream &operator<<(std::ostream &os_, const file &f_);
};

/** Stream operator with @param os_ the output stream and @param o_ the streamed
 * object */
inline std::ostream &operator<<(std::ostream &os_, const file &f_) {
    // Do the viewBox adjustment
    std::array<scalar, 2> x_range = {std::numeric_limits<scalar>::max(),
                                     std::numeric_limits<scalar>::min()};
    std::array<scalar, 2> y_range = {std::numeric_limits<scalar>::max(),
                                     std::numeric_limits<scalar>::min()};

    std::array<scalar, 4> viewBox = {-800, -800, 1600, 1600};

    for (auto &o : f_._objects) {
        x_range[0] = std::min(x_range[0], o._x_range[0]);
        x_range[1] = std::max(x_range[1], o._x_range[1]);
        y_range[0] = std::min(y_range[0], o._y_range[0]);
        y_range[1] = std::max(y_range[1], o._y_range[1]);
    }
    // Enlarge the view box by 10 percent
    viewBox[2] = 1.2 * (x_range[1] - x_range[0]);
    viewBox[3] = 1.2 * (y_range[1] - y_range[0]);

    // Include a fixed size border
    viewBox[0] = (x_range[0] - 0.1 * viewBox[2]) - f_._border;
    viewBox[1] = (y_range[0] - 0.1 * viewBox[3]) - f_._border;
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

    for (auto &o : f_._objects) {
        os_ << o;
    }
    os_ << f_._svg_tail;
    if (f_._add_html) {
        os_ << f_._html_tail;
    }
    return os_;
}

}  // namespace svg
}  // namespace actsvg
