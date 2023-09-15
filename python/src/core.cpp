// This file is part of the actsvg packge.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <fstream>
#include <memory>

#include "actsvg/core/svg.hpp"
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

}
}  // namespace python
}  // namespace actsvg