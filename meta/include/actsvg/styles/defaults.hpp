// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "actsvg/core.hpp"

namespace actsvg {

namespace defaults {

// Empty object for displaying
static svg::object __e_object;

// Title style
static style::font __t_font;

// Sensitive surface section
static style::fill __s_fill;
static style::stroke __s_stroke;

// Measure
static style::marker __m_marker;
static style::stroke __m_stroke;
static style::stroke __m_stroke_guide;
static style::font __m_font;

// Grid information
static style::fill __g_fill;
static style::stroke __g_stroke;

// Axis section
static style::stroke __a_stroke;
static style::font __a_font;
static std::array<style::marker, 2u> __a_markers = __standard_axis_markers;

// Transform section
static style::transform __t_identity;

/** Static method to create the defaults in situ */
static bool create_defaults() {

    // Title font
    __t_font._size = 14;

    // Sensitive fill and stroke
    __s_fill._fc._rgb = {66, 182, 245};
    __s_fill._fc._opacity = 0.75;
    __s_fill._fc._highlight = {"mouseover", "mouseout"};
    __s_fill._fc._hl_rgb = {245, 182, 66};

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
    __g_fill._fc._opacity = 0.25;
    __g_stroke._sc._rgb = {255, 0, 0};
    __g_stroke._width = 0.5;
    __g_stroke._hl_width = 4;

    // Transform
    __t_identity = style::transform();

    return true;
}

static bool __defaults_set = create_defaults();

}  // namespace defaults

}  // namespace actsvg
