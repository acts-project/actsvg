// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "defs.hpp"
#include "generators.hpp"
#include "style.hpp"
#include "svg.hpp"
#include "utils.hpp"

namespace actsvg {

namespace detail {
/** Helper method to estimate ranges of an object
 *
 * @param _o_ is the svg object in question (to be adapted)
 * @param vertices_ are the input vertices
 **/
void adapt_range(svg::object &_o_, const std::vector<point2> &vertices_);

}  // namespace detail

/** The draw namespace encapsulates the CORRECT
 * left handed x-y system from the AWKWARD SVG draw system, and applies
 * the transform to it if necessary
 *
 * That is its main purpose (togehter with given allowing to
 * give an identifier an connect objects.
 *
 * */

namespace draw {

/** Method to draw a simple line
 *
 * @note will perform the y switch
 *
 * @param id_ is the identification tag of this line
 * @param start_ is the start point of the line
 * @param end_ is the end point of the line
 * @param stroke_ are the stroke parameters
 * @param transform_ is an optional transform of the object
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the line
 */
svg::object line(const std::string &id_, const point2 &start_,
                 const point2 &end_,
                 const style::stroke &stroke_ = style::stroke(),
                 const style::transform &transform_ = style::transform());

/** Method to draw an arc
 *
 * @param id_ is the identification tag of this line
 * @param r_ the radius
 * @param start_ is the start point of the line
 * @param end_ is the end point of the line
 * @param stroke_ are the stroke parameters
 * @param transform_ is an optional transform of the object
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the line
 */
svg::object arc(const std::string &id_, scalar r_, const point2 &start_,
                const point2 &end_, const style::fill &fill_ = style::fill(),
                const style::stroke &stroke_ = style::stroke(),
                const style::transform &transform_ = style::transform());

/** Draw a Polyline path
 *
 * @param id_ is the identification
 * @param points_ are the points of the polyline
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @return an object representing the polyline
 */
svg::object polyline(const std::string &id_, const std::vector<point2> &points_,
                     const style::stroke &stroke = style::stroke(),
                     const style::transform &transform_ = style::transform());

/** Draw a Bezier path
 * input is defined as pairs of position, momentum
 *
 * @param id_ is the identification
 * @param xds_ is the collection of point, direction tuples
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @return a group object representing the curve
 */
svg::object bezier(const std::string &id_,
                   const std::vector<std::array<point2, 2u>> &xds_,
                   const style::stroke &stroke_ = style::stroke(),
                   const style::transform &transform_ = style::transform());

/** Draw a circle object
 *  - will translate into ellipse to allow for a scale
 *
 * @param id_ is the identification
 * @param p_ the position
 * @param r_ the radius
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the circle
 */
svg::object circle(const std::string &id_, const point2 &p_, scalar r_,
                   const style::fill &fill_ = style::fill(),
                   const style::stroke &stroke_ = style::stroke(),
                   const style::transform &transform_ = style::transform());

/** Draw an ellipse object
 *
 * @param id_ is the identification
 * @param p_ the position
 * @param rs_ the radii
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the ellipse
 */
svg::object ellipse(const std::string &id_, const point2 &p_,
                    const std::array<scalar, 2> &rs_,
                    const style::fill &fill_ = style::fill(),
                    const style::stroke &stroke_ = style::stroke(),
                    const style::transform &transform_ = style::transform());

/** Draw a polygon object
 *
 * @param id_ is the identification
 * @param polygon_ the polygon points
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 * @param apply_transform_ is the option to either apply the transform
 * or not
 *
 * @note transform is directly applied and not attached as property
 * for raw drawing objects
 *
 * @return an svg object for the polygon
 */
svg::object polygon(const std::string &id_, const std::vector<point2> &polygon_,
                    const style::fill &fill_ = style::fill(),
                    const style::stroke &stroke_ = style::stroke(),
                    const style::transform &transform_ = style::transform(),
                    bool apply_transform_ = true);

/** Draw a rectangle object
 *
 * @param id_ is the rectangle object id
 * @param c_ is the center position
 * @param half_x is the halflength in x
 * @param half_y is the halflength in y
 * @param fill_ are the fill parameters
 * @param stroke_ are the stroke parameters
 * @param transform_ defines the rectangle transform
 *
 * @return an svg object for the rectangle
 */
svg::object rectangle(const std::string &id_, const point2 &c_, scalar half_x,
                      scalar half_y, const style::fill &fill_ = style::fill{},
                      const style::stroke &stroke_ = style::stroke(),
                      const style::transform &transform_ = style::transform());

/** Draw a text object - unconnected
 *
 * @param id_ is the text object id
 * @param p_ is the text position
 * @param text_ is the actual text to be drawn
 * @param font_ is the font sytle specification
 * @param transform_ defines the text transform
 *
 * @return an svg object for the text
 *
 **/
svg::object text(const std::string &id_, const point2 &p_,
                 const std::vector<std::string> &text_,
                 const style::font &font_ = style::font(),
                 const style::transform &transform_ = style::transform());

/** Draw a label
 *
 * @param id_ is the label object id
 *
 * @return a label object
 **/
svg::object label(const std::string &id_, const style::label &l_);

/** Draw a text object - connected
 *
 * @param id_ is the text object id
 * @param p_ is the text position
 * @param text_ is the actual text to be drawn
 * @param font_ is the font sytle specification
 * @param transform_ defines the text transform
 * @param object_ is the connected object
 * @param highlight_ are the hightlighting options
 *
 * @return an svg object with highlight connection
 *
 **/
svg::object connected_text(
    const std::string &id_, const point2 &p_,
    const std::vector<std::string> &text_, const style::font &font_,
    const style::transform &transform_, const svg::object &object_,
    const std::vector<std::string> &highlight_ = {"mouseover", "mouseout"});

/** Draw an image object - connected
 *
 * @param id_ is the image object id
 * @param href_ is the image object href field
 * @param height_ is the image object height field
 * @param width_ is the image object width field
 * @param x_ is the image object x field
 * @param y_ is the image object y field
 * @param object_ is the connected object
 * @param highlight_ are the hightlighting options
 * @param onerror_ is the image object onerror field
 *
 * @return an svg object with highlight connection
 *
 **/
svg::object image_box(const std::string &id_, const std::string &href_,
                      const scalar &height_, const scalar &width_,
                      const scalar &x_, const scalar &y_,
                      const svg::object &object_,
                      const std::vector<std::string> &highlight_,
                      const std::string &onerror_);

/** Draw a text object - connected
 *
 * @param id_ is the text object id
 * @param p_ is the position of the info box
 * @param title_ is the title of the info box
 * @param title_fill_ is the fill color of the title
 * @param title_font_ is the font style of the title
 * @param text_ is the actual text to be drawn
 * @param text_fill_ is the fill color of the text box
 * @param text_font_ is the font style of the text box
 * @param stroke_ is the stroke
 * @param object_ is the connected object
 * @param highlight_ are the hightlighting options
 *
 * @return an svg object with highlight connection
 *
 **/
svg::object connected_info_box(
    const std::string &id_, const point2 &p_, const std::string &title_,
    const style::fill &title_fill_, const style::font &title_font_,
    const std::vector<std::string> &text_, const style::fill &text_fill_,
    const style::font &text_font_, const style::stroke &stroke_,
    const svg::object &object_,
    const std::vector<std::string> &highlight_ = {"mouseover", "mouseout"});

/** Draw a tiled cartesian grid - ready for connecting
 *
 * @param id_ the grid identification
 * @param l0_edges_ are the edges in l0
 * @param l1_edges_ are the edges in l1
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @return a simple cartesian grid
 */
svg::object cartesian_grid(
    const std::string &id_, const std::vector<scalar> &l0_edges_,
    const std::vector<scalar> &l1_edges_,
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform());

/** Draw a tiled cartesian grid - ready for connecting
 *
 * @param id_ the grid identification
 * @param l0_edges_ are the edges in l0
 * @param l1_edges_ are the edges in l1
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note - the single objects of the grid can be found by
 * their unique name "id_+_X_Y" when X is the bin in the
 * first local and j the bin in the second local coordinate
 *
 * @return a tiled grid with internal objects
 */
svg::object tiled_cartesian_grid(
    const std::string &id_, const std::vector<scalar> &l0_edges_,
    const std::vector<scalar> &l1_edges_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform());

/** Draw a simple fan grid
 *
 * @param id_ the grid identification
 * @param x_low_edges_ are the edges in x at low y
 * @param x_high_edges_ are the edges in x at high y
 * @param y_edges_ are the edges in phi
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @return a simple faned grid structure
 */
svg::object fan_grid(
    const std::string &id_, const std::vector<scalar> &x_low_edges_,
    const std::vector<scalar> &x_high_edges_,
    const std::vector<scalar> &y_edges_,
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) noexcept(false);

/** Draw a simple fan grid
 *
 * @param id_ the grid identification
 * @param x_low_edges_ are the edges in x at low y
 * @param x_high_edges_ are the edges in x at high y
 * @param y_edges_ are the edges in phi
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note - the single objects of the grid can be found by
 * their unique name "id_+_X_Y" when X is the bin in the
 * first local and j the bin in the second local coordinate
 *
 * @return a tiled grid with contained individual cells
 */
svg::object tiled_fan_grid(
    const std::string &id_, const std::vector<scalar> &x_low_edges_,
    const std::vector<scalar> &x_high_edges_,
    const std::vector<scalar> &y_edges_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform()) noexcept(false);

/** Draw a simple polar grid
 *
 * @param id_ the grid identification
 * @param r_edges_ are the edges in r
 * @param phi_edges_ are the edges in phi
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @return a simple polar grid object
 */
svg::object polar_grid(const std::string &id_,
                       const std::vector<scalar> &r_edges_,
                       const std::vector<scalar> &phi_edges_,
                       const style::stroke &stroke_ = style::stroke(),
                       const style::transform &transform_ = style::transform());

/** Draw a connected polar grid
 *
 * @param id_ the grid identification
 * @param r_edges_ are the edges in r
 * @param phi_edges_ are the edges in phi
 * @param fill_ is the fill style
 * @param stroke_ is the stroke style
 * @param transform_ is the optional transform
 *
 * @note - the single objects of the grid can be found by
 * their unique name "id_+_X_Y" when X is the bin in the
 * first local and j the bin in the second local coordinate
 *
 * @return a tiled polar grid in individual objects
 */
svg::object tiled_polar_grid(
    const std::string &id_, const std::vector<scalar> &r_edges_,
    const std::vector<scalar> &phi_edges_,
    const style::fill &fill_ = style::fill(),
    const style::stroke &stroke_ = style::stroke(),
    const style::transform &transform_ = style::transform());

/** Marker definition
 *
 *  Arrorws types are: <, <<, <|, |<, |<<, |<|, o, x, *
 * @param id_ is the marker identification
 * @param at_ is the poistion of the marker
 * @param marker_ is the marker style
 * @param rot_ is the rotation in [pi,phi)]
 *
 * @return an svg object for the marker group
 **/
svg::object marker(const std::string &id_, const point2 &at_,
                   const style::marker &marker_, scalar rot_ = 0.);

/** Draw a measure
 *
 * @param id_ is the identification tag of this object
 * @param start_ is the start point of the line
 * @param end_ is the end point of the line
 * @param stroke_ are the stroke parameters
 * @param start_marker_ are the marker parmeters at start
 * @param end_marker_ are the marker parmeters at start
 * @param font_ are the font parameters
 * @param label_ is the label associated
 * @param label_pos_ is the label position
 *
 * @return an svg object for the measure
 */
svg::object measure(const std::string &id_, const point2 &start_,
                    const point2 &end_,
                    const style::stroke &stroke_ = style::stroke(),
                    const style::marker &start_marker_ = style::marker({"|<"}),
                    const style::marker &end_marker_ = style::marker({"|<"}),
                    const style::font &font_ = style::font(),
                    const std::string &label_ = "",
                    const point2 &label_pos_ = {0., 0.});

/** Draw an arc measure
 *
 * @param id_ is the identification tag of this object
 * @param r_ the radius
 * @param start_ is the start point of the line, per definition
 * with smaller phi
 * @param end_ is the end point of the line, defines the marker
 * @param stroke_ are the stroke parameters
 * @param start_marker_ are the marker parmeters
 * @param end_marker_ are the marker parmeters
 * @param font_ are the font parameters
 * @param label_ is the label associated
 * @param label_pos_ is the label position
 *
 * @return an svg object for the measure
 */
svg::object arc_measure(
    const std::string &id_, scalar r_, const point2 &start_, const point2 &end_,
    const style::stroke &stroke_ = style::stroke(),
    const style::marker &start_marker_ = style::marker(),
    const style::marker &end_marker_ = style::marker({"|<"}),
    const style::font &font_ = style::font(), const std::string &label_ = "",
    const point2 &label_pos_ = {0., 0.});

/** Draw an arrow
 *
 * @param id_ is the identification tag of this object
 * @param start_ is the start point of the line
 * @param end_ is the end point of the line
 * @param stroke_ are the stroke parameters
 * @param start_marker_ are the marker parmeters at start
 * @param end_marker_ are the marker parmeters at start
 *
 * @return an svg object for the arrow
 */
svg::object arrow(const std::string &id_, const point2 &start_,
                  const point2 &end_,
                  const style::stroke &stroke_ = style::stroke(),
                  const style::marker &start_marker_ = style::marker({"|<"}),
                  const style::marker &end_marker_ = style::marker({"|<"}));

/** Draw an x-y axes system
 *
 * @param id_ is the id tag of this object
 * @param x_range_ is the x range of the axes to be drawn
 * @param y_range_ is the y range of the axes to be drawn
 * @param stroke_ are the stroke parameters
 * @param x_label_ is the x label of the axis system
 * @param y_label_ is the y label of the axis system
 * @param font_ are the font parameters
 * @param markers_ are the 4 markers on each axis end
 *
 * @return an svg object representing the axes
 */
svg::object x_y_axes(const std::string &id_,
                     const std::array<scalar, 2> &x_range_,
                     const std::array<scalar, 2> &y_range_,
                     const style::stroke &stroke_ = style::stroke(),
                     const std::string &x_label_ = "",
                     const std::string &y_label_ = "",
                     const style::font &font_ = style::font(),
                     const style::axis_markers<2u> &markers_ = {
                         __standard_axis_markers, __standard_axis_markers});

/** Helper method to create a gradient object
 *
 * @param g_ is the gradient definition
 *
 * @return the gradient as an object
 **/
svg::object gradient_as_object(const style::gradient &g_);

/** Place a gradient box
 *
 * @param id_ is the identification of this surface
 * @param p_ is the position of the box
 * @param w_h_ is the width and height of the box
 * @param stops_ is the gradient definition
 * @param label_ is the label string (optional)
 * @param font_ is the font of the box
 * @param transform_ is the transform of the box
 *
 * @note make sure the gradient is added to the file definition
 * @note this will flip the x-axis
 *
 * @return the create object
 */
svg::object gradient_box(
    const std::string &id_, const point2 &p_, const std::array<scalar, 2> &w_h_,
    const std::vector<std::tuple<style::gradient::stop, scalar>> &stops_,
    const style::label &label_ = style::label{},
    const style::stroke &stroke_ = style::stroke(),
    const style::font &font_ = style::font{},
    const style::transform t_ = style::transform());

/** Place a copy with different attributes via xling
 *
 * @param id_ the identification of this surface
 * @param ro_ the reference object
 * @param f_ the fill of the new object
 * @param s_ the stroke of the new object
 * @param t_ the transform of the new object
 *
 * @return the create object from reference
 */
svg::object from_template(const std::string &id_, const svg::object &ro_,
                          const style::fill &f_ = style::fill(),
                          const style::stroke &s_ = style::stroke(),
                          const style::transform t_ = style::transform());

}  // namespace draw

}  // namespace actsvg
