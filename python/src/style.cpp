// This file is part of the actsvg packge.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core/style.hpp"

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "actsvg/styles/defaults.hpp"
#include "utilities.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

namespace actsvg {
namespace python {

void add_style_module(context& ctx) {

    auto& m = ctx.get("main");
    auto s = m.def_submodule("style");

    {
        // The color python class
        /// @TODO add set_highlight() with self setting
        py::class_<style::color, std::shared_ptr<style::color>>(s, "color")
            .def(py::init<>())
            .def(py::init(
                [](const std::array<int, 3>& rgb, scalar opacity = 1.) {
                    style::color c;
                    c._rgb = rgb;
                    c._opacity = opacity;
                    return c;
                }));
    }

    {
        // The fill python class
        py::class_<style::fill, std::shared_ptr<style::fill>>(s, "fill")
            .def(py::init<>())
            .def(py::init<style::color>(), py::arg("color"))
            .def(py::init<bool>(), py::arg("sterile"))
            .def(py::init([](const std::array<int, 3>& rgb, scalar opacity) {
                     auto color = style::color{rgb};
                     color._opacity = opacity;
                     return style::fill{color};
                 }),
                 py::arg("rgb"), py::arg("opacity") = 1.);
    }

    {
        // The stroke python class
        py::class_<style::stroke, std::shared_ptr<style::stroke>>(s, "stroke")
            .def(py::init<>())
            .def(py::init<style::color, scalar, std::vector<int>>(),
                 py::arg("color"), py::arg("width"), py::arg("dash"))
            .def(py::init<bool>(), py::arg("sterile"))
            .def(py::init([](const std::array<int, 3>& rgb, scalar width,
                             std::vector<int> dash) {
                     return style::stroke{style::color{rgb}, width, dash};
                 }),
                 py::arg("rgb"), py::arg("width") = 1.,
                 py::arg("dash") = std::vector<int>{});
    }

    {
        // The marker python class
        py::class_<style::marker, std::shared_ptr<style::marker>>(s, "marker")
            .def(py::init<>())
            .def(py::init<>([](const std::string& t, scalar sz,
                               const style::fill& f, const style::stroke& s) {
                     return style::marker{t, sz, f, s};
                 }),
                 py::arg("type"), py::arg("size"),
                 py::arg("fill") = style::fill{},
                 py::arg("stroke") = style::stroke{});
    }

    {

        // The font python class
        py::class_<style::font, std::shared_ptr<style::font>>(s, "font")
            .def(py::init<>())
            .def(py::init<>([](const style::color& c, const std::string& f,
                               unsigned int s, scalar l_s,
                               const std::string st) {
                     return style::font{c, f, s, l_s, st};
                 }),
                 py::arg("color"), py::arg("font_family"), py::arg("size"),
                 py::arg("line_spacing"), py::arg("style"))
            .def(py::init([](const std::array<int, 3>& rgb, unsigned int size) {
                     return style::font{style::color{rgb}, "Andale Mono", size,
                                        1.4, ""};
                 }),
                 py::arg("rgb"), py::arg("size"));
    }

    {
        // The gradient python class
        py::class_<style::gradient, std::shared_ptr<style::gradient>>(
            s, "gradient")
            .def(py::init<>())
            .def_readwrite("id", &style::gradient::_id)
            .def_readwrite("direction", &style::gradient::_direction)
            .def_readwrite("type", &style::gradient::_type)
            .def_readwrite("stops", &style::gradient::_stops);
    }

    {
        // The transform python class
        py::class_<style::transform, std::shared_ptr<style::transform>>(
            s, "transform")
            .def(py::init<>())
            .def_readwrite("translate", &style::transform::_tr)
            .def_readwrite("rotate", &style::transform::_rot)
            .def_readwrite("scale", &style::transform::_scale)
            .def_readwrite("skew", &style::transform::_skew);
    }
    {
        // The style defaults
        auto d = s.def_submodule("defaults");

        d.def("sensitive_fill", []() { return defaults::__s_fill; });
        d.def("sensitive_stroke", []() { return defaults::__s_stroke; });

        // Some handy defaults
        d.def("black_fill", []() { return defaults::__bl_fill; });
        d.def("black_stroke", []() { return defaults::__bl_stroke; });

        d.def("red_fill", []() { return defaults::__r_fill; });
        d.def("red_stroke", []() { return defaults::__r_stroke; });

        d.def("white_fill", []() { return defaults::__w_fill; });
        d.def("white_stroke", []() { return defaults::__w_stroke; });

        d.def("font", []() { return defaults::__t_font; });

        d.def("axis_markers", []() { return defaults::__a_markers; });
    }
}
}  // namespace python
}  // namespace actsvg
