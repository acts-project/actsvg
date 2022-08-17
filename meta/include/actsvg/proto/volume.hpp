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
#include "actsvg/proto/surface.hpp"
#include "actsvg/styles/defaults.hpp"
#include "grid.hpp"

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

    /// Type enumeration
    enum type { e_barrel = 0, e_endcap = 1, e_other = 2 };

    unsigned int _index = 0u;

    /// Name of the volume
    std::string _name = "unnamed";

    type _type = e_barrel;

    /// Auxiliary information
    std::vector<std::string> _info = {};

    /// The contained surfaces as batches
    std::vector<std::vector<surface<point3_container>>> _surfaces = {};

    /// The portals
    std::vector<surface<point3_container>> _portals = {};

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
};

}  // namespace proto

}  // namespace actsvg