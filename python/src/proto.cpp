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
        /// Create a polygon surface from a point3 collection
        ///
        /// @param name is the name of the created object
        /// @param pcs is the point3 collection
        /// @param f is the fill style
        /// @param s is the stroke style
        p.def(
            "create_polygon",
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
            },
            py::arg("name"), py::arg("points"), py::arg("fill"),
            py::arg("stroke"));
    }

    {
        /// The grid class
        py::class_<proto::grid>(p, "grid")
            .def(py::init<>())
            .def_readwrite("name", &proto::grid::_name)
            .def_readwrite("type", &proto::grid::_type)
            .def_readwrite("edges_0", &proto::grid::_edges_0)
            .def_readwrite("edges_1", &proto::grid::_edges_1)
            .def_readwrite("reference_r", &proto::grid::_reference_r)
            .def_readwrite("bin_ids", &proto::grid::_bin_ids)
            .def_readwrite("connections", &proto::grid::_connections)
            .def_readwrite("connection_types", &proto::grid::_connection_types)
            .def_readwrite("connection_associations",
                           &proto::grid::_connection_associations)
            .def_readwrite("fill", &proto::grid::_fill)
            .def_readwrite("stroke", &proto::grid::_stroke);
    }

    {
        // The material class
        py::class_<proto::surface_material>(p, "surface_material")
            .def(py::init<>())
            .def_readwrite("material_matrix", &proto::surface_material::_material_matrix)
            .def_readwrite("material_ranges", &proto::surface_material::_material_ranges)
            .def_readwrite("grid", &proto::surface_material::_grid);

/**
    style::gradient _gradient = style::gradient{};

    /// The gradient pos
    point2 _gradient_pos = {0, 0};

    /// Gradient bos size
    std::array<scalar, 2u> _gradient_box = {0., 0.};

    /// Gradient stroke
    style::stroke _gradient_stroke = style::stroke{};

    /// Gradient font
    style::font _gradient_font = style::font{};

    /// The info position
    point2 _info_pos_ = {0, 0};

    /// The info font
    style::font _info_font_ = style::font{};
*/
    }
}

}  // namespace python
}  // namespace actsvg