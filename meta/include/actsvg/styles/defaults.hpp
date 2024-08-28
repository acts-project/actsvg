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
extern svg::object __e_object;

// Title style
extern style::font __t_font;

// Sensitive surface section, including backside
extern style::fill __s_fill;
extern style::stroke __s_stroke;

extern style::fill __bs_fill;
extern style::stroke __bs_stroke;

// Support surface section
extern style::fill __ss_fill;
extern style::stroke __ss_stroke;

// Measure
extern style::marker __m_marker;
extern style::stroke __m_stroke;
extern style::stroke __m_stroke_guide;
extern style::font __m_font;

// Grid information
extern style::fill __g_fill;
extern style::stroke __g_stroke;

// Axis section
extern style::stroke __a_stroke;
extern style::font __a_font;
extern std::array<style::marker, 2u> __a_markers;

// Background panel
extern style::fill __bg_fill;
extern style::stroke __bg_stroke;

// Transform section
extern style::transform __t_identity;

// White fill, white stroke
extern style::fill __w_fill;
extern style::stroke __w_stroke;

// Red fill, red stroke
extern style::fill __r_fill;
extern style::stroke __r_stroke;

// Black fill, black stroke
extern style::fill __bl_fill;
extern style::stroke __bl_stroke;

// Black, dashed stroke
extern style::stroke __bl_dashed_stroke;
// Black, dashed stroke
extern style::stroke __bl_dotted_stroke;

// No fill, no stroke
extern style::fill __nn_fill;
extern style::stroke __nn_stroke;

// A default gradient: r,g,b
extern style::gradient __rgb_gradient;

extern bool __defaults_set;

}  // namespace defaults

}  // namespace actsvg
