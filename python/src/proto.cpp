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

#include <memory>

#include "actsvg/meta.hpp"
#include "utilities.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

namespace actsvg {
namespace python {

using surface = proto::surface<point3_collection>;

/// @brief  Adding the proto module to the
/// python bindings
///
/// @param ctx the python context
void add_proto_module(context& ctx) {

    auto& m = ctx.get("main");
    auto p = m.def_submodule("proto");

    {
        // The python surface class
        py::class_<surface, std::shared_ptr<surface>>(p, "surface")
            .def(py::init<>());
    }

    {

        p.def("create_polygon",
              [](const std::string& name, const point3_collection& pcs,
                 const style::fill& f, const style::stroke& s) {
                  // Create the surface
                  surface sf{};
                  sf._name = name;
                  sf._type = surface::type::e_polygon;
                  sf._vertices = pcs;
                  sf._fill = f;
                  sf._stroke = s;                  
                  return sf;
              }, py::arg("name"), py::arg("points"), py::arg("fill"), py::arg("stroke"));
    }
}

}  // namespace python
}  // namespace actsvg