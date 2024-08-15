// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/styles/defaults.hpp"

#include "actsvg/core.hpp"

namespace actsvg {

namespace defaults {
svg::object __e_object;
style::font __t_font;
style::fill __s_fill;
style::stroke __s_stroke;
style::fill __bs_fill;
style::stroke __bs_stroke;
style::fill __ss_fill;
style::stroke __ss_stroke;
style::marker __m_marker;
style::stroke __m_stroke;
style::stroke __m_stroke_guide;
style::font __m_font;
style::fill __g_fill;
style::stroke __g_stroke;
style::stroke __a_stroke;
style::font __a_font;
style::fill __bg_fill;
style::stroke __bg_stroke;
style::transform __t_identity;
style::fill __w_fill;
style::stroke __w_stroke;
style::fill __r_fill;
style::stroke __r_stroke;
style::fill __bl_fill;
style::stroke __bl_stroke;
style::stroke __bl_dashed_stroke;
style::stroke __bl_dotted_stroke;
style::fill __nn_fill;
style::stroke __nn_stroke;
style::gradient __rgb_gradient;
std::array<style::marker, 2u> __a_markers = __standard_axis_markers;

/** Static method to create the defaults in situ */
static bool create_defaults() {

    // Title font
    __t_font._size = 14;

    // Sensitive fill and stroke
    __s_fill._fc._rgb = {66, 182, 245};
    __s_fill._fc._opacity = 0.75;
    __s_fill._fc._highlight = {"mouseover", "mouseout"};
    __s_fill._fc._hl_rgb = {245, 182, 66};

    __bs_fill._fc._rgb = {40, 83, 237};
    __bs_fill._fc._opacity = 0.75;
    __bs_fill._fc._highlight = {"mouseover", "mouseout"};
    __bs_fill._fc._hl_rgb = {237, 83, 40};

    __ss_fill._fc._rgb = {86, 90, 112};
    __ss_fill._fc._opacity = 0.75;
    __ss_fill._fc._highlight = {"mouseover", "mouseout"};
    __ss_fill._fc._hl_rgb = {76, 153, 84};

    // Stroke definition
    __s_stroke._sc._opacity = 1.;
    __s_stroke._width = 0.75;

    // Measurement
    __m_stroke = style::stroke();
    __m_stroke_guide = style::stroke();
    __m_stroke_guide._dasharray = {1, 1};
    __m_marker = style::marker({"|<"});

    // Grid
    __g_fill._fc._rgb = {200, 200, 200};
    __g_fill._fc._opacity = 0.;
    __g_stroke._sc._rgb = {255, 0, 0};
    __g_stroke._width = 0.5;
    __g_stroke._hl_width = 4;

    // Background panel
    __bg_fill._fc._rgb = {235, 235, 235};
    __bg_fill._fc._opacity = 90.;
    __bg_stroke._sc._rgb = {205, 205, 205};
    __bg_stroke._width = 0.5;

    // Transform
    __t_identity = style::transform();

    // Whites
    __w_fill._fc._rgb = {255, 255, 255};
    __w_stroke._sc._rgb = {255, 255, 255};

    // Reds
    __r_fill._fc._rgb = {1, 0, 0};
    __r_stroke._sc._rgb = {1, 0, 0};

    // Blacks
    __bl_fill._fc._rgb = {0, 0, 0};
    __bl_stroke._sc._rgb = {0, 0, 0};

    // Black dashed
    __bl_dashed_stroke._sc._rgb = {0, 0, 0};
    __bl_dashed_stroke._dasharray = {3, 3};

    // Black dashed
    __bl_dotted_stroke._sc._rgb = {0, 0, 0};
    __bl_dotted_stroke._dasharray = {1, 1};

    // Nulls
    __nn_fill = style::fill();

    __nn_stroke = style::stroke();
    __nn_stroke._sterile = true;
    __nn_stroke._width = 0.;

    __rgb_gradient = style::gradient();
    __rgb_gradient._stops = {{0., style::color{style::rgb{255, 0, 0}}},
                             {0.5, style::color{style::rgb{0, 255, 0}}},
                             {1., style::color{style::rgb{0, 0, 255}}}};

    return true;
}

bool __defaults_set = create_defaults();

}  // namespace defaults

}  // namespace actsvg
