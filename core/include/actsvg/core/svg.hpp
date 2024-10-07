// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
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

    /** Nested summary struct
     *
     * This struct allows for defining some summary
     * information with which this object has been
     * create.
     *
     **/
    struct summary {
        std::array<scalar, 2> _scale = {1., 1.};
    };

    using bounding_box = std::array<std::array<scalar, 2u>, 2u>;

    /// SVG tag of the objec
    std::string _tag = "";

    /// Identification string
    std::string _id = "";

    /// Auxiliary info for this object
    std::vector<std::string> _aux_info = {};

    /// Sterile - does not write _fill, stroke
    bool _sterile = false;

    /// View is active - is set to false, it is not written out
    bool _active = true;

    /// In-field of the svg tag
    std::vector<std::string> _field = {};
    scalar _field_span = 12;

    /// Fill specifications
    style::fill _fill = style::fill();

    /// Stroke specifications
    style::stroke _stroke = style::stroke();

    /// Transform of the object
    style::transform _transform = style::transform{};

    /// Attribute map for writing
    std::map<std::string, std::string> _attribute_map = {};

    /// Contains sub objects ( group, animation, etc.)
    std::vector<object> _sub_objects = {};

    /// Contains definintions
    std::vector<object> _definitions = {};

    /// Barycenter from vertices
    std::array<scalar, 2> _barycenter = {0., 0.};

    /// Range in x - view frame
    std::array<scalar, 2> _x_range = {std::numeric_limits<scalar>::max(),
                                      std::numeric_limits<scalar>::lowest()};
    /// Range in y - view frame
    std::array<scalar, 2> _y_range = {std::numeric_limits<scalar>::max(),
                                      std::numeric_limits<scalar>::lowest()};
    /// Range in r  - view frame
    std::array<scalar, 2> _r_range = {std::numeric_limits<scalar>::max(),
                                      std::numeric_limits<scalar>::lowest()};
    /// Range in phi - view frame
    std::array<scalar, 2> _phi_range = {std::numeric_limits<scalar>::max(),
                                        std::numeric_limits<scalar>::lowest()};

    /// Summary object
    summary _summary = summary{};

    /// Static constructor for a group
    static object create_group(const std::string &id_) {
        object g;
        g._tag = "g";
        g._id = id_;
        return g;
    }

    /// An object is defined if a tag is set
    bool is_defined() const;

    /// The object is an empty group
    bool is_empty_group() const;

    /// Generate bounding box
    ///
    /// @note this is in the logical frame of the object and not the view frame
    ///
    /// @return a bounding box as left lower corner and upper right corner
    bounding_box generate_bounding_box() const {
        bounding_box bb = {};
        bb[0u] = {_x_range[0], std::min(-_y_range[0], -_y_range[1])};
        bb[1u] = {_x_range[1], std::max(-_y_range[0], -_y_range[1])};
        return bb;
    }

    /** Add a sub object and respect the min/max range
     *
     * @param o_ is the object in question
     **/
    void add_object(const svg::object &o_);

    /** Find an object by string,
     *
     * @param id_ is the identifier with which the sub object
     * is searched for
     *
     * @note std::optional is used and std::nullopt as @return if
     * the search was not successful
     */
    std::optional<svg::object> find_object(const std::string &id_) const;

    /** Add a sub object and respect the min/max range
     *
     * @param oc_ is the object container to be added
     **/
    template <typename object_container>
    void add_objects(const object_container &oc_) {
        // Add the object
        for (const auto &o : oc_) {
            if (o._active) {
                add_object(o);
            }
        }
    }

    friend std::ostream &operator<<(std::ostream &os_, const object &o_);
};

/** Stream operator with @param os_ the output stream and @param o_ the streamed
 * object */
std::ostream &operator<<(std::ostream &os_, const object &o_);

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

    std::string _svg_def_end =
        " xmlns=\"http://www.w3.org/2000/svg\"  "
        "xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

    std::string _svg_tail = "</svg>\n";
    std::string _html_tail = "</body>\n</html>\n";

    scalar _width = 800;
    scalar _height = 800;
    scalar _border = 10;

    std::array<scalar, 2> _x_range = {std::numeric_limits<scalar>::max(),
                                      std::numeric_limits<scalar>::min()};
    std::array<scalar, 2> _y_range = {std::numeric_limits<scalar>::max(),
                                      std::numeric_limits<scalar>::min()};

    std::optional<std::array<scalar, 4>> _view_box = std::nullopt;

    std::vector<object> _objects = {};

    /** Add an object and respect the min/max range
     *
     * @param o_ is the object in question
     **/
    void add_object(const svg::object &o_);

    /** Add an object and respect the min/max range
     *
     * @param os_ is the vector of objects
     **/
    void add_objects(const std::vector<svg::object> &os_);

    /** Get the view box from the x-y range */
    std::array<scalar, 4> view_box() const;

    /** Set a view box, this will overwrite the auto-determined one */
    void set_view_box(const std::array<scalar, 4> &vb_);

    /** Write to ostream */
    friend std::ostream &operator<<(std::ostream &os_, const file &f_);
};

/** Stream operator with @param os_ the output stream and @param o_ the streamed
 * object */
std::ostream &operator<<(std::ostream &os_, const file &f_);
}  // namespace svg
}  // namespace actsvg
