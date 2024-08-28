// This file is part of the actsvg packge.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core.hpp"

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <fstream>
#include <memory>

#include "utilities.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

namespace actsvg {
namespace python {

/// @brief Add the core module to the
/// @param ctx
void add_core_module(context& ctx) {

    auto& m = ctx.get("main");
    auto c = m.def_submodule("core");

    {
        // The object python class
        py::class_<svg::object, std::shared_ptr<svg::object> >(c, "object")
            .def(py::init<>())
            .def(py::init([](const std::string& tag, const std::string& id) {
                svg::object o;
                o._tag = tag;
                o._id = id;
                return o;
            }))
            .def("add_object", &svg::object::add_object);
    }

    auto d = m.def_submodule("draw");

    {
        /// Draw a line - name the arguments
        d.def("line", &draw::line);

        /// Draw an arc
        d.def("arc", &draw::arc);

        /// Draw bezier curve
        d.def("bezier", &draw::bezier);

        /// Draw a circle
        d.def("circle", &draw::circle);

        /// Draw an ellipse
        d.def("ellipse", &draw::ellipse);

        /// Draw a polygon
        d.def("polygon", &draw::polygon);

        /// Draw a rectangle
        d.def("rectangle", &draw::rectangle);

        /// Draw text
        d.def("text", &draw::text);

        /// Draw connectected text
        d.def("connected_text", &draw::connected_text);

        /// Draw a connected info box
        d.def("connected_info_box", &draw::connected_info_box);

        d.def("image_box", &draw::image_box);

        /// Draw a cartesian grid
        d.def("cartesian_grid", &draw::cartesian_grid);

        /// Draw a tiled cartesian grid
        d.def("tiled_cartesian_grid", &draw::tiled_cartesian_grid);

        /// Draw a fan grid
        d.def("fan_grid", &draw::fan_grid);

        /// Draw a tiled fan grid
        d.def("tiled_fan_grid", &draw::tiled_fan_grid);

        /// Draw a polar grid
        d.def("polar_grid", &draw::polar_grid);

        /// Draw a tiled polar gird
        d.def("tiled_polar_grid", &draw::tiled_polar_grid);

        /// Draw a marker
        d.def("marker", &draw::marker);

        /// Draw a measure
        d.def("measure", &draw::measure);

        /// Draw an arrow
        d.def("arrow", &draw::arrow);

        /// Draw x-y axes
        d.def("x_y_axes", &draw::x_y_axes);

        /// Draw an arc measure
        d.def("arc_measure", &draw::arc_measure);

        /// Draw a gradient box
        d.def("gradient_box", &draw::gradient_box);
    }
}
}  // namespace python
}  // namespace actsvg