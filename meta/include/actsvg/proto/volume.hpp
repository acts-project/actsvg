// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>

#include "actsvg/core/style.hpp"
#include "actsvg/proto/grid.hpp"
#include "actsvg/proto/portal.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {

namespace proto {

/** A proto volume class as a simple translation layer
 * from a volume description
 *
 * @tparam point3_container a vertex description of surfaces
 **/

template <typename point3_container>
struct volume {

    using surface_type = surface<point3_container>;

    using portal_type = portal<point3_container>;

    /// Type enumeration
    enum type { e_cylinder = 0, e_cuboid = 1, e_other = 2 };

    // The index of the volume
    unsigned int _index = 0u;

    // The depth level of the volume
    unsigned int _depth_level = 0u;

    /// Name of the volume
    std::string _name = "unnamed_volume";

    // The type of the volume & its parameters
    type _type = e_cylinder;
    std::vector<scalar> _bound_values;
    point3_container _vertices = {};

    // Style parameters
    style::fill _fill;
    style::stroke _stroke = defaults::__nn_stroke;
    style::transform _transform;

    /// Auxiliary information
    std::vector<std::string> _info = {};

    /// The contained surfaces as batches - optimised for layer view
    std::vector<std::vector<surface<point3_container>>> _surfaces = {};

    /// The container surfaces as flat - opimised for volume view
    std::vector<surface<point3_container>> _v_surfaces = {};

    /// The portals
    std::vector<portal<point3_container>> _portals = {};

    /// The associated surface grid
    grid _surface_grid;

    /** These are the grid associations, when iterating through
     * the associataions are also grouped by batch in order to
     * split the association view if necessary
     *
     *  for (auto : batch){
     *    size_t i = 0;
     *    for (auto: _edges_1){
     *      for (auto: _edges_0){
     *        auto assoc = _grid_associations[i++];
     *        }
     *     }
     *   }
     *
     * The entries of the association point to the index in the surface
     * container of the volume
     **/
    std::vector<std::vector<std::vector<size_t>>> _grid_associations = {};

    /// Colorize method
    ///
    /// @param colors_ are the indexed colors
    ///
    void colorize(std::vector<style::color>& colors_) {
        if (_index < colors_.size()) {
            _fill._fc = colors_[_index];
        }
        // Colorize the portals
        for (auto& p : _portals) {
            // The portal itself has no fill color
            p._surface._fill._fc._rgb = {255, 255, 255};
            p._surface._fill._fc._opacity = 1.;
            // The links are filled
            for (auto& vl : p._volume_links) {
                if (vl._link_index < colors_.size()) {
                    vl._stroke._sc = colors_[vl._link_index];
                    vl._stroke._sc._opacity = 1.;
                    vl._end_marker._fill._fc = colors_[vl._link_index];
                    vl._end_marker._fill._fc._opacity = 1.;
                    vl._end_marker._stroke = vl._stroke;
                }
            }
        }
        // Colorize the surfaces
        for (auto& sbatch : _surfaces) {
            for (auto& s : sbatch) {
                s._fill._fc = colors_[_index];
                s._fill._fc._opacity = 0.25;
            }
        }
    }
};

}  // namespace proto

}  // namespace actsvg
