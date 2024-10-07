// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core/draw.hpp"

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "actsvg/core/generators.hpp"
#include "actsvg/core/style.hpp"
#include "actsvg/core/svg.hpp"
#include "actsvg/core/utils.hpp"

namespace actsvg {

namespace detail {
void adapt_range(svg::object &_o_, const std::vector<point2> &vertices_) {
    for (const auto &v : vertices_) {
        _o_._x_range = {std::min(_o_._x_range[0], v[0]),
                        std::max(_o_._x_range[1], v[0])};
        _o_._y_range = {std::min(_o_._y_range[0], v[1]),
                        std::max(_o_._y_range[1], v[1])};

        scalar r = std::sqrt(v[0] * v[0] + v[1] * v[1]);
        scalar phi = std::atan2(-v[1], v[0]);
        _o_._r_range = {std::min(_o_._r_range[0], r),
                        std::max(_o_._r_range[1], r)};
        _o_._phi_range = {std::min(_o_._phi_range[0], phi),
                          std::max(_o_._phi_range[1], phi)};
    }
}

}  // namespace detail

namespace draw {
svg::object line(const std::string &id_, const point2 &start_,
                 const point2 &end_, const style::stroke &stroke_,
                 const style::transform &transform_) {
    svg::object l;
    l._tag = "line";
    l._id = id_;

    // Apply the transform & scale
    // - the line needs the scale direcly applied
    scalar tx = transform_._tr[0];
    scalar ty = transform_._tr[1];
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];

    scalar st_x = sx * (start_[0] + tx);
    scalar st_y = sx * (-start_[1] - ty);
    scalar en_x = sy * (end_[0] + ty);
    scalar en_y = sy * (-end_[1] - ty);

    l._barycenter =
        utils::barycenter<std::array<scalar, 2>>({{st_x, st_y}, {en_x, en_y}});

    // Draw the line, remember the sign flip
    l._attribute_map["x1"] = utils::to_string(st_x);
    l._attribute_map["y1"] = utils::to_string(st_y);
    l._attribute_map["x2"] = utils::to_string(en_x);
    l._attribute_map["y2"] = utils::to_string(en_y);

    // Adapt the range of this object
    detail::adapt_range(l, {{st_x, st_y}, {en_x, en_y}});

    // Remember the stroke attributes and add them
    l._stroke = stroke_;
    return l;
}

svg::object arc(const std::string &id_, scalar r_, const point2 &start_,
                const point2 &end_, const style::fill &fill_,
                const style::stroke &stroke_,
                const style::transform &transform_) {
    svg::object a;
    a._tag = "path";
    a._id = id_;

    // Apply the transform & scale
    // - the arc needs the scale directly attached
    scalar tx = transform_._tr[0];
    scalar ty = transform_._tr[1];
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];

    scalar x_min = sx * (start_[0] + tx);
    scalar y_min = -sy * (start_[1] + ty);
    scalar x_max = sx * (end_[0] + tx);
    scalar y_max = -sy * (end_[1] + ty);

    std::string arc_string =
        "M " + utils::to_string(x_min) + " " + utils::to_string(y_min);
    arc_string +=
        " A " + utils::to_string(sx * r_) + " " + utils::to_string(sy * r_);
    arc_string += " 0 0 0 ";
    arc_string += utils::to_string(x_max) + " " + utils::to_string(y_max);
    a._attribute_map["d"] = arc_string;

    // Adapt the range of this object
    detail::adapt_range(a, {{x_min, y_min}, {x_max, y_max}});

    /// @todo add sagitta
    a._barycenter = utils::barycenter<std::array<scalar, 2>>(
        {{x_min, y_min}, {x_max, y_max}});

    // Remember the stroke attributes and add them
    a._stroke = stroke_;
    a._fill = fill_;
    return a;
}

svg::object polyline(const std::string &id_, const std::vector<point2> &points_,
                     const style::stroke &stroke,
                     const style::transform &transform_) {

    svg::object p;
    p._tag = "polyline";
    p._id = id_;

    // Apply the transform & scale
    scalar tx = transform_._tr[0];
    scalar ty = transform_._tr[1];
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];

    std::string svg_points_string;
    std::vector<point2> display_vertices;
    display_vertices.reserve(points_.size());
    for (auto v : points_) {
        // Add display scaling
        v[0] *= sx;
        v[1] *= sy;
        // Add scaled * translation
        v[0] += sx * tx;
        v[1] += sy * ty;
        v[1] *= -1;
        // Per vertex range estimation
        detail::adapt_range(p, {v});

        // Convert to string attributes, y-switch
        svg_points_string += utils::to_string(v[0]);
        svg_points_string += ",";
        svg_points_string += utils::to_string(v[1]);
        svg_points_string += " ";
    }

    // Barycenter
    p._barycenter = utils::barycenter(display_vertices);
    // Fill the points attributes
    p._attribute_map["points"] = svg_points_string;
    // Attach stroke
    p._stroke = stroke;

    // Return
    return p;
}

svg::object bezier(const std::string &id_,
                   const std::vector<std::array<point2, 2u>> &xds_,
                   const style::stroke &stroke_,
                   const style::transform &transform_) {

    svg::object c = svg::object::create_group(id_);

    point2 lx = {0, 0};
    point2 ld = {0, 0};

    for (const auto &[x, d] : xds_) {
        if (ld == point2{0, 0}) {
            ld = d;
            lx = x;
            continue;
        }
        // Buid the Bezier segments
        // intesect the two lines, apply transform to draw frame
        point2 intersect = transform_.apply(utils::intersect(lx, ld, x, d));
        point2 dlx = transform_.apply(lx);
        point2 dx = transform_.apply(x);

        std::string path_string =
            "M " + utils::to_string(dlx[0]) + " " + utils::to_string(dlx[1]) +
            " C " + utils::to_string(intersect[0]) + " " +
            utils::to_string(intersect[1]) + " " +
            utils::to_string(intersect[0]) + " " +
            utils::to_string(intersect[1]) + " " + utils::to_string(dx[0]) +
            " " + utils::to_string(dx[1]);

        svg::object path;
        path._tag = "path";
        path._id = id_ + "_segment_" + std::to_string(c._sub_objects.size());
        path._attribute_map["d"] = path_string;
        path._stroke = stroke_;
        c.add_object(path);
        ld = d;
        lx = x;
    }

    return c;
}

svg::object circle(const std::string &id_, const point2 &p_, scalar r_,
                   const style::fill &fill_, const style::stroke &stroke_,
                   const style::transform &transform_) {
    // Create the object, tag it, id it (if given)
    svg::object e;
    e._tag = "ellipse";
    e._id = id_;

    // Apply the transform & scale
    point2 c = transform_.apply(p_);

    // Apply the transform & scale
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];

    // Fill the points attributes
    scalar rx = r_ * sx;
    scalar ry = r_ * sy;
    e._attribute_map["cx"] = utils::to_string(c[0]);
    e._attribute_map["cy"] = utils::to_string(c[1]);
    e._attribute_map["rx"] = utils::to_string(rx);
    e._attribute_map["ry"] = utils::to_string(ry);

    // Adapt the range of this object
    detail::adapt_range(e, {{c[0] - rx, c[1] - ry}, {c[0] + rx, c[1] + ry}});
    // Assign it as new barycenter
    e._barycenter = c;

    // Attach fill, stroke & transform attributes and apply
    e._fill = fill_;
    e._stroke = stroke_;

    // The svg object is now set up
    return e;
}

svg::object ellipse(const std::string &id_, const point2 &p_,
                    const std::array<scalar, 2> &rs_, const style::fill &fill_,
                    const style::stroke &stroke_,
                    const style::transform &transform_) {

    // Create the object, tag it, id it (if given)
    svg::object e;
    e._tag = "ellipse";
    e._id = id_;

    // Apply the transform & scale
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];
    scalar cx = sx * (p_[0] + transform_._tr[0]);
    scalar cy = sy * (-p_[1] - transform_._tr[1]);

    // Fill the points attributes
    e._attribute_map["cx"] = utils::to_string(cx);
    e._attribute_map["cy"] = utils::to_string(cy);
    e._attribute_map["rx"] = utils::to_string(rs_[0] * sx);
    e._attribute_map["ry"] = utils::to_string(rs_[1] * sy);

    // Adapt the range of this object
    detail::adapt_range(e, {{cx - sx * rs_[0], cy - sy * rs_[1]},
                            {cx + sx * rs_[0], cy + sy * rs_[1]}});

    // Barycenter
    e._barycenter = {cx, cy};

    // Attach fill, stroke & transform attributes and apply
    e._fill = fill_;
    e._stroke = stroke_;
    // The svg object is now set up
    return e;
}

svg::object polygon(const std::string &id_, const std::vector<point2> &polygon_,
                    const style::fill &fill_, const style::stroke &stroke_,
                    const style::transform &transform_, bool apply_transform_)

{
    // Create the object, tag it, id it (if given)
    svg::object p;
    p._tag = "polygon";
    p._id = id_;
    // Apply the scale
    scalar tx = transform_._tr[0];
    scalar ty = transform_._tr[1];
    scalar sx = transform_._scale[0];
    scalar sy = transform_._scale[1];
    // Write attributes and measure object size, length
    std::string svg_points_string;
    std::vector<point2> display_vertices;
    display_vertices.reserve(polygon_.size());
    for (auto v : polygon_) {
        scalar alpha = transform_._rot[0];
        if (alpha != 0. and apply_transform_) {
            scalar alpha_rad = static_cast<scalar>(alpha / 180. * M_PI);
            v = utils::rotate(v, alpha_rad);
        }

        // Add display scaling
        v[0] *= sx;
        v[1] *= sy;

        point2 vr = v;
        // Add scaled * translation
        if (apply_transform_) {
            v[0] += sx * tx;
            v[1] += sy * ty;
        }
        // range is always adaped
        vr[0] += sx * tx;
        vr[1] += sy * ty;

        // respect y flip
        v[1] *= -1;
        vr[1] *= -1;

        // Per vertex range estimation - always done in display coordinates
        detail::adapt_range(p, {vr});

        // Convert to string attributes, y-switch
        svg_points_string += utils::to_string(v[0]);
        svg_points_string += ",";
        svg_points_string += utils::to_string(v[1]);
        svg_points_string += " ";
    }
    // Barycenter
    p._barycenter = utils::barycenter(display_vertices);
    // Fill the points attributes
    p._attribute_map["points"] = svg_points_string;
    // Attach fill, stroke & transform attributes and apply
    p._fill = fill_;
    p._stroke = stroke_;
    if (not apply_transform_) {
        p._transform = transform_;
    }
    // The svg object is now set up
    return p;
}

svg::object rectangle(const std::string &id_, const point2 &c_, scalar half_x,
                      scalar half_y, const style::fill &fill_,
                      const style::stroke &stroke_,
                      const style::transform &transform_) {

    svg::object r;
    r._id = id_;
    r._tag = "rect";
    // Size attributes
    r._attribute_map["x"] = utils::to_string(c_[0] - half_x);
    r._attribute_map["y"] = utils::to_string(-c_[1] - half_y);
    r._attribute_map["width"] = utils::to_string(half_x * 2);
    r._attribute_map["height"] = utils::to_string(half_y * 2);
    // Attach style attributes
    r._fill = fill_;
    r._stroke = stroke_;
    r._transform = transform_;
    // Set the range
    detail::adapt_range(r, {{c_[0] - half_x, -c_[1] - half_y},
                            {c_[0] + half_x, -c_[1] + half_y}});

    // Return the rectangle objec t
    return r;
}

svg::object text(const std::string &id_, const point2 &p_,
                 const std::vector<std::string> &text_,
                 const style::font &font_, const style::transform &transform_) {
    // Create the object, tag it, id it (if given)
    svg::object t;
    t._tag = "text";
    t._id = id_;
    // Apply the scale
    scalar x = p_[0];
    scalar y = p_[1];

    x *= transform_._scale[0];
    y *= transform_._scale[1];

    t._fill = font_._fc;

    // Fill the field
    t._field = text_;
    t._attribute_map["x"] = utils::to_string(x);
    t._attribute_map["y"] = utils::to_string(-y);
    t._attribute_map["font-family"] = font_._family;
    t._attribute_map["font-size"] = std::to_string(font_._size);
    t._field_span = font_._size * font_._line_spacing;

    // Get the longest text line
    size_t l = 0;
    for (const auto &tl : text_) {
        l = l > tl.size() ? l : tl.size();
    }

    scalar fs = font_._size;

    detail::adapt_range(t, {{x, y - fs * text_.size()},
                            {static_cast<scalar>(x + 0.7 * fs * l), y + l}});

    t._barycenter = utils::barycenter<std::array<scalar, 2>>(
        {{x, y - l}, {static_cast<scalar>(x + 0.7 * fs * l), y + l}});

    return t;
}

svg::object connected_text(const std::string &id_, const point2 &p_,
                           const std::vector<std::string> &text_,
                           const style::font &font_,
                           const style::transform &transform_,
                           const svg::object &object_,
                           const std::vector<std::string> &highlight_) {
    auto t = text(id_, p_, text_, font_, transform_);

    t._attribute_map["display"] = "none";

    svg::object on;
    on._tag = "animate";
    on._attribute_map["fill"] = "freeze";
    on._attribute_map["attributeName"] = "display";
    on._attribute_map["from"] = "none";
    on._attribute_map["to"] = "block";
    on._attribute_map["begin"] = object_._id + __d + highlight_[1];

    svg::object off;

    off._tag = "animate";
    off._attribute_map["fill"] = "freeze";
    off._attribute_map["attributeName"] = "display";
    off._attribute_map["to"] = "none";
    off._attribute_map["from"] = "block";
    off._attribute_map["begin"] = object_._id + __d + highlight_[0];

    // Store the animation
    t._sub_objects.push_back(on);
    t._sub_objects.push_back(off);
    return t;
}

svg::object image_box(const std::string &id_, const std::string &href_,
                      const scalar &height_, const scalar &width_,
                      const scalar &x_, const scalar &y_,
                      const svg::object &object_,
                      const std::vector<std::string> &highlight_,
                      const std::string &onerror_) {
    svg::object i = svg::object::create_group(id_);

    // Image box object
    svg::object imgb;
    imgb._tag = "image";
    imgb._attribute_map["href"] = href_;
    imgb._attribute_map["height"] = utils::to_string(height_);
    imgb._attribute_map["width"] = utils::to_string(width_);
    imgb._attribute_map["x"] = utils::to_string(x_);
    imgb._attribute_map["y"] = utils::to_string(y_);
    imgb._attribute_map["onerror"] = onerror_;

    i.add_object(imgb);
    // Connect it
    if (object_.is_defined()) {
        i._attribute_map["display"] = "none";

        svg::object on;
        on._tag = "animate";
        on._attribute_map["fill"] = "freeze";
        on._attribute_map["attributeName"] = "display";
        on._attribute_map["from"] = "none";
        on._attribute_map["to"] = "block";
        on._attribute_map["begin"] = object_._id + __d + highlight_[1];
        on._attribute_map["fill-opacity"] = 1;

        svg::object off;

        off._tag = "animate";
        off._attribute_map["fill"] = "freeze";
        off._attribute_map["attributeName"] = "display";
        off._attribute_map["to"] = "none";
        off._attribute_map["from"] = "block";
        off._attribute_map["begin"] = object_._id + __d + highlight_[0];

        // Store the animation
        i._sub_objects.push_back(on);
        i._sub_objects.push_back(off);
    }

    detail::adapt_range(i, {{x_, y_}, {x_ + width_, y_ + height_}});

    return i;
}

svg::object connected_info_box(
    const std::string &id_, const point2 &p_, const std::string &title_,
    const style::fill &title_fill_, const style::font &title_font_,
    const std::vector<std::string> &text_, const style::fill &text_fill_,
    const style::font &text_font_, const style::stroke &stroke_,
    const svg::object &object_, const std::vector<std::string> &highlight_) {

    svg::object ib = svg::object::create_group(id_);

    size_t tew = 0;
    for (const auto &t : text_) {
        tew = t.size() > tew ? t.size() : tew;
    }
    scalar tews = (tew + 2) * text_font_._size;
    scalar tiws = (title_.size() + 2) * title_font_._size;
    scalar ws = std::max(tews, tiws);

    scalar tih = 2 * title_font_._size;
    scalar teh =
        static_cast<scalar>(1.5 * text_.size() + 0.5) * text_font_._size;

    std::vector<std::array<scalar, 2>> tic = {p_,
                                              {p_[0], p_[1] - tih},
                                              {p_[0] + ws, p_[1] - tih},
                                              {p_[0] + ws, p_[1]}};

    auto tibox = polygon(id_ + "_title_box", tic, title_fill_, stroke_);
    ib.add_object(tibox);
    auto ti = text(id_ + "_title",
                   {p_[0] + title_font_._size,
                    static_cast<scalar>(p_[1] - 1.5 * title_font_._size)},
                   {title_}, title_font_);
    ib.add_object(ti);

    std::vector<std::array<scalar, 2>> tec = {{p_[0], p_[1] - tih},
                                              {p_[0] + ws, p_[1] - tih},
                                              {p_[0] + ws, p_[1] - tih - teh},
                                              {p_[0], p_[1] - tih - teh}};

    auto tebox = polygon(id_ + "_text_box", tec, text_fill_, stroke_);
    ib.add_object(tebox);

    scalar toffset = text_.size() == 1u ? tih + title_font_._size : tih;
    auto te =
        text(id_ + "_text",
             {p_[0] + title_font_._size, static_cast<scalar>(p_[1] - toffset)},
             text_, text_font_);
    ib.add_object(te);

    // Connect it
    if (object_.is_defined()) {
        ib._attribute_map["display"] = "none";

        svg::object on;
        on._tag = "animate";
        on._attribute_map["fill"] = "freeze";
        on._attribute_map["attributeName"] = "display";
        on._attribute_map["from"] = "none";
        on._attribute_map["to"] = "block";
        on._attribute_map["begin"] = object_._id + __d + highlight_[1];
        on._attribute_map["fill-opacity"] = 1;

        svg::object off;

        off._tag = "animate";
        off._attribute_map["fill"] = "freeze";
        off._attribute_map["attributeName"] = "display";
        off._attribute_map["to"] = "none";
        off._attribute_map["from"] = "block";
        off._attribute_map["begin"] = object_._id + __d + highlight_[0];

        // Store the animation
        ib._sub_objects.push_back(on);
        ib._sub_objects.push_back(off);
    }
    return ib;
}

svg::object cartesian_grid(const std::string &id_,
                           const std::vector<scalar> &l0_edges_,
                           const std::vector<scalar> &l1_edges_,
                           const style::stroke &stroke_,
                           const style::transform &transform_) {
    // The grid group object
    svg::object grid = svg::object::create_group(id_);

    scalar l0_min = l0_edges_[0];
    scalar l0_max = l0_edges_[l0_edges_.size() - 1];

    scalar l1_min = l1_edges_[0];
    scalar l1_max = l1_edges_[l1_edges_.size() - 1];

    for (auto [i0, l0] : utils::enumerate(l0_edges_)) {
        grid.add_object(line(id_ + "_l0_" + std::to_string(i0), {l0, l1_min},
                             {l0, l1_max}, stroke_, transform_));
    }

    for (auto [i1, l1] : utils::enumerate(l1_edges_)) {
        grid.add_object(line(id_ + "_l1_" + std::to_string(i1), {l0_min, l1},
                             {l0_max, l1}, stroke_, transform_));
    }

    return grid;
}

svg::object tiled_cartesian_grid(const std::string &id_,
                                 const std::vector<scalar> &l0_edges_,
                                 const std::vector<scalar> &l1_edges_,
                                 const style::fill &fill_,
                                 const style::stroke &stroke_,
                                 const style::transform &transform_) {
    // The grid group object
    svg::object grid = svg::object::create_group(id_);

    // The list of grid sectors
    for (size_t il0 = 1; il0 < l0_edges_.size(); ++il0) {
        // Grid svg object
        std::string gs = id_ + "_";
        gs += std::to_string(il0 - 1);
        gs += "_";
        for (size_t il1 = 1; il1 < l1_edges_.size(); ++il1) {
            std::array<scalar, 2u> llc = {l0_edges_[il0 - 1],
                                          l1_edges_[il1 - 1]};
            std::array<scalar, 2u> lrc = {l0_edges_[il0], l1_edges_[il1 - 1]};
            std::array<scalar, 2u> rrc = {l0_edges_[il0], l1_edges_[il1]};
            std::array<scalar, 2u> rlc = {l0_edges_[il0 - 1], l1_edges_[il1]};

            std::vector<std::array<scalar, 2u>> tile = {llc, lrc, rrc, rlc};
            // Create the grid tile, add auxiliary info and add it to the grid
            auto grid_tile = polygon(gs + std::to_string(il1 - 1), tile, fill_,
                                     stroke_, transform_);
            grid_tile._aux_info = {std::string("* bin (" +
                                               std::to_string(il0 - 1) + "," +
                                               std::to_string(il1 - 1) + ")")};
            grid.add_object(grid_tile);
        }
    }
    return grid;
}

svg::object fan_grid(const std::string &id_,
                     const std::vector<scalar> &x_low_edges_,
                     const std::vector<scalar> &x_high_edges_,
                     const std::vector<scalar> &y_edges_,
                     const style::stroke &stroke_,
                     const style::transform &transform_) noexcept(false) {
    // The list of grid sectors
    svg::object grid = svg::object::create_group(id_);

    if (x_low_edges_.size() != x_high_edges_.size()) {
        throw std::invalid_argument(
            "fan_grid: mismatch in low/high edge numbers");
    }

    scalar x_low_min = x_low_edges_[0];
    scalar x_low_max = x_low_edges_[x_low_edges_.size() - 1];

    scalar x_high_min = x_high_edges_[0];
    scalar x_high_max = x_high_edges_[x_high_edges_.size() - 1];

    scalar y_min = y_edges_[0];
    scalar y_max = y_edges_[y_edges_.size() - 1];

    // Calculate slopes
    scalar x_low_slope = (x_high_min - x_low_min) / (y_max - y_min);
    scalar x_high_slope = (x_high_max - x_low_max) / (y_max - y_min);

    for (const auto [ix, x] : utils::enumerate(x_low_edges_)) {
        scalar x_min = x;
        scalar x_max = x_high_edges_[ix];
        grid.add_object(line(id_ + "_l0_" + std::to_string(ix), {x_min, y_min},
                             {x_max, y_max}, stroke_, transform_));
    }

    for (const auto &[iy, y] : utils::enumerate(y_edges_)) {
        scalar x_min = x_low_min + (y - y_min) * x_low_slope;
        scalar x_max = x_low_max + (y - y_min) * x_high_slope;
        grid.add_object(line(id_ + "_l1_" + std::to_string(iy), {x_min, y},
                             {x_max, y}, stroke_, transform_));
    }

    return grid;
}

svg::object tiled_fan_grid(const std::string &id_,
                           const std::vector<scalar> &x_low_edges_,
                           const std::vector<scalar> &x_high_edges_,
                           const std::vector<scalar> &y_edges_,
                           const style::fill &fill_,
                           const style::stroke &stroke_,
                           const style::transform &transform_) noexcept(false) {
    svg::object grid = svg::object::create_group(id_);

    // The list of grid sectors
    std::vector<svg::object> grid_tiles;

    if (x_low_edges_.size() != x_high_edges_.size()) {
        throw std::invalid_argument(
            "fan_grid: mismatch in low/high edge numbers");
    }

    scalar y_min = y_edges_[0];
    scalar y_max = y_edges_[y_edges_.size() - 1];

    std::vector<scalar> slopes;
    for (auto [ix, xl] : utils::enumerate(x_low_edges_)) {
        scalar xh = x_high_edges_[ix];
        slopes.push_back((xh - xl) / (y_max - y_min));
    }

    for (auto [iy, yl] : utils::enumerate(y_edges_)) {
        if (iy + 1 == y_edges_.size()) {
            continue;
        }
        scalar yh = y_edges_[iy + 1];
        for (auto [ix, xll] : utils::enumerate(x_low_edges_)) {
            if (ix + 1 == x_low_edges_.size()) {
                continue;
            }
            scalar xlr = x_low_edges_[ix + 1];

            scalar x_l_slope = slopes[ix];
            scalar x_r_slope = slopes[ix + 1];

            // the corrected position given the y values
            scalar xll_c = xll + x_l_slope * (yl - y_min);
            scalar xlr_c = xlr + x_r_slope * (yl - y_min);

            scalar xhl_c = xll + x_l_slope * (yh - y_min);
            scalar xhr_c = xlr + x_r_slope * (yh - y_min);

            // Create the grid tile, add auxiliary info and add it
            std::string g_n =
                id_ + "_" + std::to_string(ix) + "_" + std::to_string(iy);
            auto grid_tile = draw::polygon(
                g_n, {{xll_c, yl}, {xlr_c, yl}, {xhr_c, yh}, {xhl_c, yh}},
                fill_, stroke_, transform_);
            grid_tile._aux_info = {std::string("* bin (" + std::to_string(ix) +
                                               "," + std::to_string(iy) + ")")};
            grid.add_object(grid_tile);
        }
    }

    return grid;
}

svg::object polar_grid(const std::string &id_,
                       const std::vector<scalar> &r_edges_,
                       const std::vector<scalar> &phi_edges_,
                       const style::stroke &stroke_,
                       const style::transform &transform_) {
    // The list of grid sectors
    svg::object grid = svg::object::create_group(id_);

    scalar r_min = r_edges_[0];
    scalar r_max = r_edges_[r_edges_.size() - 1];

    scalar phi_min = phi_edges_[0];
    scalar phi_max = phi_edges_[phi_edges_.size() - 1];

    scalar cos_phi_min = std::cos(phi_min);
    scalar sin_phi_min = std::sin(phi_min);
    scalar cos_phi_max = std::cos(phi_max);
    scalar sin_phi_max = std::sin(phi_max);

    style::fill fill;

    for (const auto [ir, r] : utils::enumerate(r_edges_)) {
        grid.add_object(draw::arc(id_ + "_r_" + std::to_string(ir), r,
                                  {r * cos_phi_min, r * sin_phi_min},
                                  {r * cos_phi_max, r * sin_phi_max}, fill,
                                  stroke_, transform_));
    }

    for (const auto &[iphi, phi] : utils::enumerate(phi_edges_)) {
        scalar cos_phi = std::cos(phi);
        scalar sin_phi = std::sin(phi);
        grid.add_object(draw::line(id_ + "_l_" + std::to_string(iphi),
                                   {r_min * cos_phi, r_min * sin_phi},
                                   {r_max * cos_phi, r_max * sin_phi}, stroke_,
                                   transform_));
    }
    return grid;
}

svg::object tiled_polar_grid(const std::string &id_,
                             const std::vector<scalar> &r_edges_,
                             const std::vector<scalar> &phi_edges_,
                             const style::fill &fill_,
                             const style::stroke &stroke_,
                             const style::transform &transform_) {
    svg::object grid = svg::object::create_group(id_);

    for (size_t ir = 1; ir < r_edges_.size(); ++ir) {
        // Grid svg object
        std::string gs = id_ + "_";
        gs += std::to_string(ir - 1);
        gs += "_";
        for (size_t iphi = 1; iphi < phi_edges_.size(); ++iphi) {
            auto sector_contour = generators::sector_contour(
                r_edges_[ir - 1], r_edges_[ir], phi_edges_[iphi - 1],
                phi_edges_[iphi]);
            // Create the grid tile, add auxiliary info and add it
            auto grid_tile =
                polygon(gs + std::to_string(iphi - 1), sector_contour, fill_,
                        stroke_, transform_);
            grid_tile._aux_info = {std::string("* bin (" +
                                               std::to_string(ir - 1) + "," +
                                               std::to_string(iphi - 1) + ")")};
            grid.add_object(grid_tile);
        }
    }
    return grid;
}

svg::object marker(const std::string &id_, const point2 &at_,
                   const style::marker &marker_, scalar rot_) {
    svg::object marker_group = svg::object::create_group(id_);

    std::vector<point2> arrow_head;
    auto size = marker_._size;

    // Offset due to measureing
    scalar m_offset = 0.;

    // It's a measure type
    if (marker_._type.substr(0u, 1u) == "|") {
        auto measure_line = line(
            id_ + "_line", {at_[0], static_cast<scalar>(at_[1] - 2 * size)},
            {at_[0], static_cast<scalar>(at_[1] + 2 * size)}, marker_._stroke);
        marker_group.add_object(measure_line);
        m_offset = -size;
    }
    // Still an arrow to draw
    if (marker_._type.find("<") != std::string::npos) {
        arrow_head = {{at_[0] - size + m_offset, at_[1] - size},
                      {at_[0] + size + m_offset, at_[1]},
                      {at_[0] - size + m_offset, at_[1] + size}};

        // Modify the arrow-type marker
        if (marker_._type.find("<<") != std::string::npos) {
            // Filled arrow-head
            arrow_head.push_back(
                {static_cast<scalar>(at_[0] - 0.25 * size + m_offset), at_[1]});
        } else if (marker_._type.substr(1u, marker_._type.size()).find("|") ==
                   std::string::npos) {
            arrow_head.push_back(
                {static_cast<scalar>(at_[0] + 1 * size + m_offset), at_[1]});
        }
    } else if (marker_._type.find("o") != std::string::npos) {
        // A dot marker
        svg::object dot =
            circle(id_, at_, 0.5_scalar * size, marker_._fill, marker_._stroke);
        marker_group.add_object(dot);
    } else if (marker_._type.find("x") != std::string::npos) {
        scalar a_x = at_[0];
        scalar a_y = at_[1];
        scalar h_s = 0.5_scalar * size;
        marker_group.add_object(line(id_ + "_ml0", {a_x - h_s, a_y - h_s},
                                     {a_x + h_s, a_y + h_s}, marker_._stroke));
        marker_group.add_object(line(id_ + "_ml1", {a_x - h_s, a_y + h_s},
                                     {a_x + h_s, a_y - h_s}, marker_._stroke));
    }

    // Plot the arrow if not empty
    if (not arrow_head.empty()) {
        auto arrow = polygon(id_, arrow_head, marker_._fill, marker_._stroke);
        marker_group.add_object(arrow);
    }

    // We need to rotate the marker group
    if (rot_ != 0.) {
        style::transform group_transform = style::transform{};
        group_transform._rot = {static_cast<scalar>(rot_ * 180. / M_PI), at_[0],
                                at_[1]};
        marker_group._transform = group_transform;
    }

    return marker_group;
}

svg::object measure(const std::string &id_, const point2 &start_,
                    const point2 &end_, const style::stroke &stroke_,
                    const style::marker &start_marker_,
                    const style::marker &end_marker_, const style::font &font_,
                    const std::string &label_, const point2 &label_pos_) {
    // Measure group here we go
    svg::object measure_group = svg::object::create_group(id_);

    auto mline = line(id_ + "_line", start_, end_, stroke_);
    measure_group.add_object(mline);

    // Calculate the rotation
    scalar theta = std::atan2(end_[1] - start_[1], end_[0] - start_[0]);
    if (std::abs(end_[1] - start_[1]) <
        std::numeric_limits<scalar>::epsilon()) {
        theta = (end_[0] > start_[0]) ? 0.0_scalar : pi;
    }

    measure_group.add_object(marker(id_ + "_start_tag", {start_[0], start_[1]},
                                    start_marker_, pi + theta));
    measure_group.add_object(marker(id_ + "_end_tag", {end_[0], end_[1]},
                                    end_marker_, static_cast<scalar>(theta)));

    if (not label_.empty()) {
        auto ltext = text(id_ + "_label", label_pos_, {label_}, font_);
        measure_group.add_object(ltext);
    }
    return measure_group;
}

svg::object arc_measure(const std::string &id_, scalar r_, const point2 &start_,
                        const point2 &end_, const style::stroke &stroke_,
                        const style::marker &start_marker_,
                        const style::marker &end_marker_,
                        const style::font &font_, const std::string &label_,
                        const point2 &label_pos_) {
    // Measure group here we go
    svg::object measure_group = svg::object::create_group(id_);

    measure_group.add_object(
        arc((id_ + "_arc"), r_, start_, end_, style::fill(), stroke_));

    // Arrow is at end point
    if (not start_marker_._type.empty() and
        start_marker_._type != std::string("none")) {
        scalar theta_start = static_cast<scalar>(atan2(start_[1], start_[0]));
        measure_group.add_object(
            marker(id_ + "_start_tag", {start_[0], start_[1]}, start_marker_,
                   static_cast<scalar>(theta_start - 0.5 * M_PI)));
    }

    scalar theta_end = static_cast<scalar>(atan2(end_[1], end_[0]));
    measure_group.add_object(
        marker(id_ + "_end_tag", {end_[0], end_[1]}, end_marker_,
               static_cast<scalar>(theta_end + 0.5_scalar * pi)));

    if (not label_.empty()) {
        auto ltext = text(id_ + "_label", label_pos_, {label_}, font_);
        measure_group.add_object(ltext);
    }

    return measure_group;
}

svg::object arrow(const std::string &id_, const point2 &start_,
                  const point2 &end_, const style::stroke &stroke_,
                  const style::marker &start_marker_,
                  const style::marker &end_marker_) {

    return measure(id_, start_, end_, stroke_, start_marker_, end_marker_);
}

svg::object x_y_axes(const std::string &id_,
                     const std::array<scalar, 2> &x_range_,
                     const std::array<scalar, 2> &y_range_,
                     const style::stroke &stroke_, const std::string &x_label_,
                     const std::string &y_label_, const style::font &font_,
                     const style::axis_markers<2u> &markers_) {
    svg::object axes = svg::object::create_group(id_);

    auto x =
        line(id_ + "_x_axis", {x_range_[0], 0.}, {x_range_[1], 0.}, stroke_);
    auto y =
        line(id_ + "_y_axis", {0., y_range_[0]}, {0., y_range_[1]}, stroke_);

    axes.add_object(x);
    axes.add_object(y);

    /** Helper method to add marker heads
     *
     * @param p_ is the position of the marker
     * @param b0_ and @param b1_ are the accessors into the marker styles
     * @param rot_ is the rotation parameter
     * @param mid_ is the marker identification
     *
     * */
    auto add_marker = [&](const point2 &p_, unsigned int b0_, unsigned int b1_,
                          scalar rot_, const std::string &mid_) -> void {
        auto lmarker = markers_[b0_][b1_];
        if (not lmarker._type.empty()) {
            axes.add_object(marker(mid_, {p_[0], p_[1]}, lmarker, rot_));
        }
    };

    // Add the markers to the arrows
    add_marker({x_range_[0], 0.}, 0, 0, pi, id_ + "_neg_x_head");
    add_marker({x_range_[1], 0.}, 0, 1, 0., id_ + "pos_x_head");
    add_marker({0., y_range_[0]}, 1, 0, -0.5_scalar * pi, id_ + "neg_y_head");
    add_marker({0., y_range_[1]}, 1, 1, 0.5_scalar * pi, id_ + "pos_y_head");

    // Add the labels: x
    if (not x_label_.empty()) {
        scalar size = markers_[0][1]._size;
        auto xlab = text(id_ + "_x_label", {x_range_[1] + 2 * size, size},
                         {x_label_}, font_);
        axes.add_object(xlab);
    }
    // Add the labels: y
    if (not y_label_.empty()) {
        scalar size = markers_[1][1]._size;
        auto ylab = text(id_ + "_y_label", {-size, y_range_[1] + 2 * size},
                         {y_label_}, font_);
        axes.add_object(ylab);
    }

    return axes;
}

svg::object gradient_as_object(const style::gradient &g_) {
    svg::object grad;
    grad._id = g_._id;
    grad._tag = g_._type + "Gradient";
    grad._sterile = true;
    grad._attribute_map["x1"] = utils::to_string(g_._direction[0] * 100) + "%";
    grad._attribute_map["y1"] = utils::to_string(g_._direction[1] * 100) + "%";
    grad._attribute_map["x2"] = utils::to_string(g_._direction[2] * 100) + "%";
    grad._attribute_map["y2"] = utils::to_string(g_._direction[3] * 100) + "%";
    for (const auto &[s, c] : g_._stops) {
        svg::object stop;
        stop._tag = "stop";
        stop._sterile = true;
        stop._attribute_map["offset"] = utils::to_string(s * 100) + "%";
        stop._attribute_map["stop-color"] = style::rgb_attr(c._rgb);
        grad._sub_objects.push_back(stop);
    }
    return grad;
}

svg::object gradient_box(
    const std::string &id_, const point2 &p_, const std::array<scalar, 2> &w_h_,
    const std::vector<std::tuple<style::gradient::stop, scalar>> &stops_,
    const style::label &label_, const style::stroke &stroke_,
    const style::font &font_, const style::transform t_) {

    svg::object g = svg::object::create_group(id_);
    g._sterile = true;

    auto [w_, h_] = w_h_;

    style::gradient gradient;
    gradient._id = id_ + "_gradient";
    gradient._type = "linear";

    std::vector<scalar> gradient_values;
    for (auto [s, o] : stops_) {
        gradient._stops.push_back(s);
    }
    bool vertical = w_ < h_;
    // vertical
    if (vertical) {
        gradient._direction = {0., 1., 0., 0.};
    } else {
        // horizontal
        gradient._direction = {0., 0., 1., 0.};
    }
    g.add_object(gradient_as_object(gradient));

    bool label_on = not label_._text.empty();

    // Create the box - by hand this time
    svg::object gbox = draw::rectangle(
        id_ + "_box", {p_[0u] + 0.5_scalar * w_, p_[1u] + 0.5_scalar * h_},
        0.5_scalar * w_, 0.5_scalar * h_);
    gbox._sterile = true;
    gbox._attribute_map["fill"] = "url(#" + gradient._id + ")";
    gbox._transform = t_;

    g.add_object(gbox);

    if (label_on) {
        style::label draw_label = label_;
        auto [llc, urc] = gbox.generate_bounding_box();
        draw_label.place(llc, urc);
        g.add_object(draw::label(id_ + "_label", draw_label));
    }

    // Create the tick objects
    for (auto [is, s] : utils::enumerate(gradient._stops)) {
        // Create the tick
        svg::object tick =
            vertical ? line(id_ + "_tick_" + std::to_string(is),
                            {p_[0], p_[1] + h_ * s.first},
                            {p_[0] + w_, p_[1] + h_ * s.first}, stroke_)
                     : line(id_ + "_tick_" + std::to_string(is),
                            {p_[0] + w_ * s.first, p_[1]},
                            {p_[0] + w_ * s.first, p_[1] + h_}, stroke_);
        g.add_object(tick);

        // Create the tick value
        scalar value = std::get<1>(stops_[is]);

        svg::object tval =
            vertical
                ? draw::text(
                      id_ + "_tick_val_" + std::to_string(is),
                      {static_cast<scalar>(p_[0] + w_ + 0.2 * font_._size),
                       p_[1] + h_ * s.first},
                      {utils::to_string(value)}, font_)
                : draw::text(id_ + "_tick_val_" + std::to_string(is),
                             {p_[0] + w_ * s.first,
                              static_cast<scalar>(p_[1] - 1.2 * font_._size)},
                             {utils::to_string(value)}, font_);
        tval._sterile = true;
        tval._attribute_map["alignment-baseline"] = "middle";
        if (not vertical) {
            tval._attribute_map["text-anchor"] = "middle";
        }
        font_.attach_attributes(tval);
        tval._field = {utils::to_string(value)};
        g.add_object(tval);
    }

    return g;
}

svg::object from_template(const std::string &id_, const svg::object &ro_,
                          const style::fill &f_, const style::stroke &s_,
                          const style::transform t_) {
    // Create new svg object
    svg::object nsvg = svg::object::create_group(id_);
    nsvg._sterile = true;

    // Refer to as the linker object
    svg::object use_obj;
    use_obj._tag = "use";
    use_obj._id = id_ + "_use";
    use_obj._attribute_map["xlink:href"] = "#" + ro_._id;
    use_obj._definitions.push_back(ro_);

    // Barycenter is shifted
    use_obj._barycenter = {(ro_._barycenter[0] + t_._tr[0]),
                           (ro_._barycenter[1] + t_._tr[1])};

    // Set the fill attributes
    use_obj._fill = f_;
    use_obj._stroke = s_;
    use_obj._transform = t_;

    // Add it to the main object
    nsvg.add_object(use_obj);

    return nsvg;
}

svg::object label(const std::string &id_, const style::label &l_) {
    auto t_o = text(id_, l_._position, {l_._text}, l_._font);
    return t_o;
}

}  // namespace draw

}  // namespace actsvg
