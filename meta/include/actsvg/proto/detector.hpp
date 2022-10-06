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

#include "volume.hpp"

namespace actsvg {

namespace proto {

/** A proto volume class as a simple translation layer
 * from a volume description
 *
 * @tparam point3_container a vertex description of surfaces
 **/

template <typename point3_container>
struct detector {

    using volume_type = volume<point3_container>;

    /// Name of the volume
    std::string _name = "unnamed_detector";

    /// Auxiliary information
    std::vector<volume<point3_container>> _volumes = {};

    /// Colorize method
    ///
    /// @param colors_ are the indexed colors
    ///
    void colorize(std::vector<style::color>& colors_) {
        for (auto& v : _volumes) {
            v.colorize(colors_);
        }
    }
};

}  // namespace proto

}  // namespace actsvg