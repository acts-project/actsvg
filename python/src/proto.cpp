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
        auto s =
            py::class_<surface, std::shared_ptr<surface>>(p, "surface")
                .def(py::init<>())
                .def_readwrite("name", &surface::_name)
                .def_readwrite("type", &surface::_type)
                .def_readwrite("vertices", &surface::_vertices)
                .def_readwrite("fill", &surface::_fill)
                .def_readwrite("stroke", &surface::_stroke)
                .def_static(
                    "from_polygon",
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
        auto g = py::class_<proto::grid>(p, "grid")
                     .def(py::init<>())
                     .def_readwrite("name", &proto::grid::_name)
                     .def_readwrite("type", &proto::grid::_type)
                     .def_readwrite("edges_0", &proto::grid::_edges_0)
                     .def_readwrite("edges_1", &proto::grid::_edges_1)
                     .def_readwrite("reference_r", &proto::grid::_reference_r)
                     .def_readwrite("bin_ids", &proto::grid::_bin_ids)
                     .def_readwrite("connections", &proto::grid::_connections)
                     .def_readwrite("connection_types",
                                    &proto::grid::_connection_types)
                     .def_readwrite("connection_associations",
                                    &proto::grid::_connection_associations)
                     .def_readwrite("fill", &proto::grid::_fill)
                     .def_readwrite("stroke", &proto::grid::_stroke);

        /// The grid type enum
        py::enum_<proto::grid::type>(g, "grid_type")
            .value("x_y", proto::grid::type::e_x_y)
            .value("r_phi", proto::grid::type::e_r_phi)
            .value("z_phi", proto::grid::type::e_z_phi)
            .export_values();
    }

    {
        // The material slab class
        py::class_<proto::material_slab>(p, "material_slab")
            .def(py::init<>())
            .def(py::init([](const std::array<scalar, 6u>& properties) {
                proto::material_slab slab;
                slab._properties = properties;
                return slab;
            }))
            .def(py::init(
                [](const std::array<scalar, 5u>& material, scalar thickness) {
                    proto::material_slab slab;
                    slab._properties[0] = material[0];
                    slab._properties[1] = material[1];
                    slab._properties[2] = material[2];
                    slab._properties[3] = material[3];
                    slab._properties[4] = material[4];
                    slab._properties[5] = thickness;
                    return slab;
                }))
            .def_readwrite("properties", &proto::material_slab::_properties);

        // The material class
        auto sm =
            py::class_<proto::surface_material>(p, "surface_material")
                .def(py::init<>())
                .def("evaluate_material_ranges",
                     &proto::surface_material::evaluate_material_ranges)
                .def_readwrite("material_matrix",
                               &proto::surface_material::_material_matrix)
                .def_readwrite("material_ranges",
                               &proto::surface_material::_material_ranges)
                .def_readwrite("grid", &proto::surface_material::_grid)
                .def_readwrite("gradient", &proto::surface_material::_gradient)
                .def_readwrite("gradient_pos",
                               &proto::surface_material::_gradient_pos)
                .def_readwrite("gradient_box",
                               &proto::surface_material::_gradient_box)
                .def_readwrite("gradient_stroke",
                               &proto::surface_material::_gradient_stroke)
                .def_readwrite("gradient_font",
                               &proto::surface_material::_gradient_font)
                .def_readwrite("info_pos", &proto::surface_material::_info_pos)
                .def_readwrite("info_font",
                               &proto::surface_material::_info_font);

        using grid_index = std::pair<unsigned int, unsigned int>;
        using local_material_index = std::pair<grid_index, unsigned int>;

        sm.def("fill_indexed_material",
               [](proto::surface_material& psm,
                  const std::vector<proto::material_slab>& material_vector,
                  const std::vector<local_material_index>&
                      local_material_indices) {
                   for (const auto& [li, mi] : local_material_indices) {
                       unsigned int b0 = li.first - 1;
                       unsigned int b1 = li.second - 1;
                       psm._material_matrix[b1][b0] = material_vector[mi];
                   }
                   return psm;
               });
    }
}

}  // namespace python
}  // namespace actsvg