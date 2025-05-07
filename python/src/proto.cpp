// This file is part of the actsvg package.
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

using volume = proto::volume<point3_collection>;

using portal = proto::portal<point3_collection>;

using trajectory = proto::trajectory<point3>;

using seed = proto::seed<point3>;

using track_state = proto::track_state<surface>;

using channel1 = proto::channel<1u>;

using channel2 = proto::channel<2u>;

using cluster1 = proto::cluster<1u>;

using cluster2 = proto::cluster<2u>;

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
            .def(py::init<>())
            .def_readwrite("_name", &surface::_name)
            .def_readwrite("_type", &surface::_type)
            .def_readwrite("_vertices", &surface::_vertices)
            .def_readwrite("_fill", &surface::_fill)
            .def_readwrite("_stroke", &surface::_stroke)
            .def_static(
                "polygon_from_vertices",
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
                py::arg("name"), py::arg("vertices"), py::arg("fill"),
                py::arg("stroke"))
            .def_static(
                "polygon_from_vertices_and_transform",
                [](const std::string& name, const point3_collection& pcs,
                   const point3& translation,
                   const std::array<point3, 3u>& rotation, const style::fill& f,
                   const style::stroke& s) {
                    // Create the surface
                    surface sf{};
                    sf._name = name;
                    auto identity = surface::transform3::identity();
                    if (translation != identity._translation ||
                        rotation != identity._rotation) {
                        // Set the predefined transform
                        sf._surface_transform =
                            surface::transform3{translation, rotation};
                    }
                    sf._vertices = pcs;
                    sf._fill = f;
                    sf._stroke = s;
                    return sf;
                },
                py::arg("name"), py::arg("vertices"), py::arg("translation"),
                py::arg("rotation"), py::arg("fill"), py::arg("stroke"))
            .def_static(
                "annulus_from_bounds_and_transform",
                [](const std::string& name, const std::vector<scalar>& bounds,
                   const point3& translation,
                   const std::array<point3, 3u>& rotation, const style::fill& f,
                   const style::stroke& s) {
                    // Create the surface
                    surface sf{};
                    sf._name = name;
                    sf._type = surface::type::e_annulus;
                    sf._measures = bounds;
                    auto identity = surface::transform3::identity();
                    if (translation != identity._translation ||
                        rotation != identity._rotation) {
                        // Set the predefined transform
                        sf._surface_transform =
                            surface::transform3{translation, rotation};
                    }
                    sf._fill = f;
                    sf._stroke = s;

                    return sf;
                },
                py::arg("name"), py::arg("bounds"), py::arg("translation"),
                py::arg("rotation"), py::arg("fill"), py::arg("stroke"));
    }

    {
        // The python portal definition
        py::class_<portal, std::shared_ptr<portal>>(p, "portal")
            .def(py::init<>())
            .def_readwrite("_name", &portal::_name)
            .def_readwrite("_surface", &portal::_surface);
    }

    {
        // The pythond volume definition
        py::class_<volume, std::shared_ptr<volume>>(p, "volume")
            .def(py::init<>())
            .def_readwrite("_name", &volume::_name);
    }

    {
        /// The grid class
        auto g = py::class_<proto::grid>(p, "grid")
                     .def(py::init<>())
                     .def_readwrite("_name", &proto::grid::_name)
                     .def_readwrite("_type", &proto::grid::_type)
                     .def_readwrite("_edges_0", &proto::grid::_edges_0)
                     .def_readwrite("_edges_1", &proto::grid::_edges_1)
                     .def_readwrite("_reference_r", &proto::grid::_reference_r)
                     .def_readwrite("_bin_ids", &proto::grid::_bin_ids)
                     .def_readwrite("_connections", &proto::grid::_connections)
                     .def_readwrite("_connection_types",
                                    &proto::grid::_connection_types)
                     .def_readwrite("_connection_associations",
                                    &proto::grid::_connection_associations)
                     .def_readwrite("_fill", &proto::grid::_fill)
                     .def_readwrite("_stroke", &proto::grid::_stroke);

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
            .def_readwrite("_properties", &proto::material_slab::_properties);

        // The material class
        auto sm =
            py::class_<proto::surface_material>(p, "surface_material")
                .def(py::init<>())
                .def("evaluate_material_ranges",
                     &proto::surface_material::evaluate_material_ranges)
                .def_readwrite("_material_matrix",
                               &proto::surface_material::_material_matrix)
                .def_readwrite("_material_ranges",
                               &proto::surface_material::_material_ranges)
                .def_readwrite("_grid", &proto::surface_material::_grid)
                .def_readwrite("_gradient", &proto::surface_material::_gradient)
                .def_readwrite("_gradient_pos",
                               &proto::surface_material::_gradient_pos)
                .def_readwrite("_gradient_box",
                               &proto::surface_material::_gradient_box)
                .def_readwrite("_gradient_stroke",
                               &proto::surface_material::_gradient_stroke)
                .def_readwrite("_gradient_font",
                               &proto::surface_material::_gradient_font)
                .def_readwrite("_info_pos", &proto::surface_material::_info_pos)
                .def_readwrite("_info_font",
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

    {
        // The channel class: 1D
        py::class_<channel1>(p, "channel1")
            .def(py::init<>())
            .def_readwrite("_cid", &channel1::_cid)
            .def_readwrite("_data", &channel1::_data);

        // The cluster class: 1D
        py::class_<cluster1>(p, "cluster1")
            .def(py::init<>())
            .def_readwrite("_type", &cluster1::_type)
            .def_readwrite("_coords", &cluster1::_coords)
            .def_readwrite("_channels", &cluster1::_channels)
            .def_readwrite("_measurement", &cluster1::_measurement)
            .def_readwrite("_variance", &cluster1::_variance)
            .def_readwrite("_correlation", &cluster1::_correlation)
            .def_readwrite("_truth", &cluster1::_truth)
            .def_readwrite("_mc", &cluster1::_mc);

        // The channel class: 2D
        py::class_<channel2>(p, "channel2")
            .def(py::init<>())
            .def_readwrite("_cid", &channel2::_cid)
            .def_readwrite("_data", &channel2::_data);

        // The cluster class: 2D
        py::class_<cluster2>(p, "cluster2")
            .def(py::init<>())
            .def_readwrite("_type", &cluster2::_type)
            .def_readwrite("_coords", &cluster2::_coords)
            .def_readwrite("_channels", &cluster2::_channels)
            .def_readwrite("_measurement", &cluster2::_measurement)
            .def_readwrite("_variance", &cluster2::_variance)
            .def_readwrite("_correlation", &cluster2::_correlation)
            .def_readwrite("_truth", &cluster2::_truth)
            .def_readwrite("_mc", &cluster2::_mc);
    }

    {
        // The trajectory class
        py::class_<trajectory, std::shared_ptr<trajectory>>(p, "trajectory")
            .def(py::init<>())
            .def_readwrite("_origin", &trajectory::_origin)
            .def_readwrite("_direction", &trajectory::_direction)
            .def_readwrite("_path", &trajectory::_path)
            .def_readwrite("_origin_size", &trajectory::_origin_size)
            .def_readwrite("_origin_stroke", &trajectory::_origin_stroke)
            .def_readwrite("_origin_fill", &trajectory::_origin_fill)
            .def_readwrite("_path_arrow", &trajectory::_path_arrow)
            .def_readwrite("_path_stroke", &trajectory::_path_stroke);
    }

    {
        // The seed class
        py::class_<seed, std::shared_ptr<seed>>(p, "seed")
            .def(py::init<>())
            .def_readwrite("_space_points", &seed::_space_points)
            .def_readwrite("_trajectory", &seed::_trajectory)
            .def_readwrite("_space_point_size", &seed::_space_point_size)
            .def_readwrite("_space_point_fill", &seed::_space_point_fill)
            .def_readwrite("_space_point_stroke", &seed::_space_point_stroke);
    }
}

}  // namespace python
}  // namespace actsvg
