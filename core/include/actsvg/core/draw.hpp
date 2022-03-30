// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "defs.hpp"
#include "svg.hpp"
#include "style.hpp"
#include "utils.hpp"
#include "generators.hpp"

#include <vector>
#include <string>
#include <map>

namespace actsvg
{

    namespace draw
    {
        /** Draw a polygon object
         *
         * @param polygon_ the polygon points
         * @param id_ is the identification
         * @param fill_ is the fill style
         * @param stroke_ is the stroke style
         * @param transform_ is the optional transform
         */
        static svg::object polygon(const std::vector<point2> &polygon_,
                            const std::string &id_ = "",
                            const style::fill &fill_ = style::fill(),
                            const style::stroke &stroke_ = style::stroke(),
                            const style::transform &transform_ = style::transform())

        {
            // Create the object, tag it, id it (if given)
            svg::object p;
            p._tag = "polygon";
            p._id = id_;
            // Apply the scale
            scalar sx = transform_._scale[0];
            scalar sy = transform_._scale[1];
            // Write attributes and measure object size, length
            std::string svg_points_string;
            for (auto [x, y] : polygon_)
            {
                x *= sx;
                y *= sy;
                // Record min/max for the view point
                p._x_range = {std::min(p._x_range[0], x),
                              std::max(p._x_range[1], x)};
                p._y_range = {std::min(p._y_range[0], y),
                              std::max(p._y_range[1], y)};
                // Add them up
                p._barycenter[0] += x;
                p._barycenter[1] += y;
                // Convert to string attributes
                svg_points_string += std::to_string(x);
                svg_points_string += ",";
                svg_points_string += std::to_string(y);
                svg_points_string += " ";
            }
            // Re-normalize the barycenter
            p._barycenter[0] /= polygon_.size();
            p._barycenter[1] /= polygon_.size();
            // Fill the points attributes
            p._attribute_map["points"] = svg_points_string;
            // Attach fill, stroke & transform attributes and apply
            p._fill = fill_;
            p._stroke = stroke_;
            p._transform = transform_;
            fill_.attach_attributes(p);
            stroke_.attach_attributes(p);
            transform_.attach_attributes(p);
            // The svg object is now set up
            return p;
        }

        /** Draw a text object - unconnected
         *
         * @param x_ is the x position
         * @param y_ is sthe y position
         * @param tid_ is the text object id
         * @param text_ is the actual text to be drawn
         * @param font_ is the font sytle specification
         * @param transform_ defines the text transform
         **/
        static svg::object text(scalar x_, scalar y_,
                         const std::string &tid_,
                         const std::vector<std::string> &text_,
                         const style::font &font_ = style::font(),
                         const style::transform &transform_ = style::transform())
        {

            // Create the object, tag it, id it (if given)
            svg::object t;
            t._tag = "text";
            t._id = tid_;
            // Apply the scale
            x_ *= transform_._scale[0];
            y_ *= transform_._scale[1];
            // Fill the field
            t._field = text_;
            t._attribute_map["x"] = std::to_string(x_);
            t._attribute_map["y"] = std::to_string(y_);
            // barycenter
            t._barycenter = {x_, y_};

            font_.attach_attributes(t);
            transform_.attach_attributes(t);
            return t;
        }

        /** Draw a text object - unconnected
         *
         * @param x_ is the x position
         * @param y_ is sthe y position
         * @param tid_ is the text object id
         * @param text_ is the actual text to be drawn
         * @param font_ is the font sytle specification
         * @param transform_ defines the text transform
         * @param object_ is the connected object
         * @param highlight_ are the hightlighting options
         **/
        static svg::object connected_text(scalar x_, scalar y_,
                                   const std::string &tid_,
                                   const std::vector<std::string> &text_,
                                   const style::font &font_,
                                   const style::transform &transform_,
                                   const svg::object &object_,
                                   const std::vector<std::string> &highlight_ = {"mouseover", "mouseout"})
        {

            auto t = text(x_, y_, tid_, text_, font_, transform_);

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
            t._animations.push_back(on);
            t._animations.push_back(off);
            return t;
        }

        /** Draw a connected r-phi-grid in x_y
         *
         * @param r_edges_ are the edges in r
         * @param phi_edges_ are the edges in phi
         * @param fill_ is the fill style
         * @param stroke_ is the stroke style
         * @param transform_ is the optional transform
         */
        static std::vector<svg::object> r_phi_grid(const std::vector<scalar> &r_edges_,
                                            const std::vector<scalar> &phi_edges_,
                                            const style::fill &fill_ = style::fill(),
                                            const style::stroke &stroke_ = style::stroke(),
                                            const style::transform &transform_ = style::transform())
        {
            // Apply the scale
            scalar sx = transform_._scale[0];
            scalar sy = transform_._scale[1];
            // The list of grid sectors
            std::vector<svg::object> grid_sectors;
            for (size_t ir = 1; ir < r_edges_.size(); ++ir)
            {
                // Grid svg object
                std::string gs = "g_r";
                gs += std::to_string(ir - 1);
                gs += "_phi";
                for (size_t iphi = 1; iphi < phi_edges_.size(); ++iphi)
                {
                    auto sector_contour = generators::sector_contour(r_edges_[ir - 1], r_edges_[ir],
                                                                     phi_edges_[iphi - 1], phi_edges_[iphi]);
                    if (sx != 0. and sy != 0.)
                    {
                        for (auto sc : sector_contour)
                        {
                            sc[0] *= sx;
                            sc[1] *= sy;
                        }
                    }

                    auto grid_sector = polygon(sector_contour, gs + std::to_string(iphi - 1), fill_, stroke_, transform_);
                    scalar r = r_edges_[ir - 1], r_edges_[ir];
                    scalar phi = phi_edges_[iphi - 1], phi_edges_[iphi];
                    grid_sector._barycenter = {sx * r * std::cos(phi), sy * r * std::sin(phi)};
                    grid_sectors.push_back(grid_sector);
                }
            }
            return grid_sectors;
        }

        /** Draw a connected z-phi-grid in x_y
         *
         * @param z_edges_ are the edges in r
         * @param phi_edges_ are the edges in phi
         * @param fill_ is the fill style
         * @param stroke_ is the stroke style
         * @param transform_ is the optional transform
         */
        static std::vector<svg::object> z_phi_grid(const std::vector<scalar> &z_edges_,
                                            const std::vector<scalar> &phi_edges_,
                                            const style::fill &fill_ = style::fill(),
                                            const style::stroke &stroke_ = style::stroke(),
                                            const style::transform &transform_ = style::transform())
        {

            // Apply the scale
            scalar sx = transform_._scale[0];
            scalar sy = transform_._scale[1];
            // The list of grid sectors
            std::vector<svg::object> grid_tiles;
            for (size_t iz = 1; iz < z_edges_.size(); ++iz)
            {
                // Grid svg object
                std::string gs = "g_z";
                gs += std::to_string(iz - 1);
                gs += "_phi";
                for (size_t iphi = 1; iphi < phi_edges_.size(); ++iphi)
                {
                    std::array<scalar, 2u> llc = {sx * z_edges_[iz - 1], sy * phi_edges_[iphi - 1]};
                    std::array<scalar, 2u> lrc = {sx * z_edges_[iz], sy * phi_edges_[iphi - 1]};
                    std::array<scalar, 2u> rrc = {sx * z_edges_[iz], sy * phi_edges_[iphi]};
                    std::array<scalar, 2u> rlc = {sx * z_edges_[iz - 1], sy * phi_edges_[iphi]};

                    std::vector<std::array<scalar, 2u>> tile = {llc, lrc, rrc, rlc};

                    auto grid_tile = polygon(tile, gs + std::to_string(iphi - 1), fill_, stroke_, transform_);
                    grid_tiles.push_back(grid_tile);
                }
            }
            return grid_tiles;
        }

    } // namespace draw

} // namespace actsvg
