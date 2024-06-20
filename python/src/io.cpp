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

#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"
#include "utilities.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

namespace actsvg {
namespace python {

using surface = proto::surface<point3_collection>;

/// @brief  Adding the proto module to the python bindings
///
/// @param ctx the python context
void add_io_module(context& ctx) {

    auto& m = ctx.get("main");
    // The io module: actsvg.io
    auto io = m.def_submodule("io");

    {
        // The file python class
        py::class_<svg::file, std::shared_ptr<svg::file>>(io, "file")
            .def(py::init<>())
            .def("add_object", &svg::file::add_object)
            .def("add_objects", &svg::file::add_objects)
            .def("write", [](svg::file& self, const std::string& fn) {
                std::ofstream fout;
                fout.open(fn);
                fout << self;
                fout.close();
            });
    }

    // The obj submodule: actsvg.io.obj
    auto obj = io.def_submodule("obj");

    // Adding a reader function from obj files for polygons
    {
        /// This function reads polygons from obj files and colors from mtl
        /// files
        ///
        /// @param objfiles the '.obj' files describing the geometry
        /// @param mtlfile the '.mtl' files for the color
        /// @param read_as_template read as template, object + transform
        obj.def("read_polygons", [](const std::vector<std::string>& objfiles,
                                    const std::vector<std::string>& mtlfiles) {
            // The return polygon surfaces
            std::vector<surface> polygons = {};

            // Read the material definitions from the mtl files
            std::map<std::string, style::fill> fill_colors;
            for (const auto& mn : mtlfiles) {
                // Input file
                std::ifstream file(mn.c_str());
                if (file.is_open()) {
                    std::string line;
                    std::string current_material;
                    char del = ' ';
                    while (std::getline(file, line)) {
                        // Check for material
                        if (line.substr(0u, 6u) == std::string("newmtl")) {
                            auto materials =
                                split_string<to_string>(line, del, {0});
                            if (not materials.empty()) {
                                current_material = materials.front();
                                if (fill_colors.find(current_material) ==
                                    fill_colors.end()) {
                                    fill_colors[current_material] =
                                        style::fill{};
                                }
                            }
                        }
                        // Add R/G/B
                        if (line.substr(0, 2u) == std::string("Kd")) {
                            auto& fill_color = fill_colors[current_material];
                            auto rgb_rel =
                                split_string<to_scalar>(line, del, {0});
                            if (rgb_rel.size() == 3u) {
                                int r = int(rgb_rel[0] * 256);
                                int g = int(rgb_rel[1] * 256);
                                int b = int(rgb_rel[2] * 256);
                                fill_color._fc._rgb = {r, g, b};
                            }
                        }
                    }
                    file.close();
                }
            }

            for (const auto& objn : objfiles) {
                // Input file
                std::ifstream file(objn.c_str());
                if (file.is_open()) {
                    std::string line;
                    // Peak loop for allocation
                    unsigned int vCounter = 0;
                    unsigned int fCounter = 0;
                    // Material strings
                    while (std::getline(file, line)) {
                        // Check for vertices
                        if (line.substr(0u, 1u) == std::string("v")) {
                            ++vCounter;
                        }
                        // Check for faces
                        if (line.substr(0u, 1u) == std::string("f")) {
                            ++fCounter;
                        }
                    }
                    file.clear();
                    file.seekg(0, std::ios::beg);

                    // The vertices to be read in
                    std::vector<point3> vertices;
                    vertices.reserve(vCounter + 1);
                    // one offset at obj counts from 1
                    vertices.push_back({0, 0, 0});

                    // The poligons to be contructed - proxies first
                    std::vector<std::vector<int>> polygon_proxies;
                    polygon_proxies.reserve(fCounter);
                    polygons.reserve(fCounter);

                    while (std::getline(file, line)) {

                        char del = ' ';
                        if (line.substr(0u, 1u) == std::string("v")) {
                            auto vertex =
                                split_string<to_scalar>(line, del, {0});
                            point3 p;
                            std::copy_n(vertex.begin(), 3, p.begin());
                            vertices.push_back(p);
                        }
                        if (line.substr(0u, 1u) == std::string("f")) {
                            auto pp = split_string<to_int>(line, del, {0});
                            polygon_proxies.push_back(pp);
                        }
                    }
                    // Construct the surfaces
                    for (auto [ip, pp] : utils::enumerate(polygon_proxies)) {

                        point3_collection surface_vertices;
                        surface_vertices.reserve(pp.size());
                        for (const auto& pv : pp) {
                            surface_vertices.push_back(vertices[pv]);
                        }
                        // add to the polygons
                        surface s;
                        s._name = "polygon_obj_" + std::to_string(ip);
                        s._type = surface::type::e_polygon;
                        s._vertices = surface_vertices;
                        // Add to the polygons
                        polygons.push_back(std::move(s));
                    }
                    file.close();
                }
            }
            // Return the surfaces
            return polygons;
        });
    }

    {
        obj.def("read_steps", [](const std::vector<std::string>& objfiles) {
            // The return polygon surfaces
            std::vector<std::vector<point3>> tracks = {};

            for (const auto& objn : objfiles) {
                // Input file
                std::ifstream file(objn.c_str());
                if (file.is_open()) {
                    std::string line;
                    char del = ' ';
                    // Material strings
                    std::vector<point3> track_steps;
                    while (std::getline(file, line)) {
                        // Check for vertices
                        if (line.substr(0u, 1u) == std::string("v")) {
                            auto vertex =
                                split_string<to_scalar>(line, del, {0});
                            point3 p;
                            std::copy_n(vertex.begin(), 3, p.begin());
                            track_steps.push_back(p);
                        } else if (not track_steps.empty()) {
                            // flush
                            tracks.push_back(track_steps);
                            track_steps.clear();
                        }
                    }
                }
            }
            return tracks;
        });
    }
}

}  // namespace python
}  // namespace actsvg
