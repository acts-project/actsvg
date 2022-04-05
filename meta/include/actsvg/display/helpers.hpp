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

#include "actsvg/actsvg.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {

using namespace defaults;

namespace display {

/** Helper method to prepare axis for a view point
 *
 * @param first_ is the first axis
 * @param second_ is the second axis
 * @param sx_ is the first axis scale
 * @param sy_ is the second axis scale
 * @param ax_ is the first axis addon
 * @param ay_ is the second axis addon
 *
 * @return the marker size as 1 percent of the range
 **/
static void prepare_axes(std::array<scalar, 2>& first_,
                         std::array<scalar, 2>& second_, scalar sx_, scalar sy_,
                         scalar ax_ = 0., scalar ay_ = 0.) {
    // Add some extra space for the axis
    first_[0] *= sx_;
    first_[1] *= sx_;
    first_[0] -= ax_;
    first_[1] += ax_;

    second_[0] *= sy_;
    second_[1] *= sy_;
    second_[0] -= ay_;
    second_[1] += ay_;
}

/** Helper method to scale the axis accordingly
 *
 * @param cc_ is an iterable container of contours
 *
 * @return a view range
 **/
template <typename contour_container>
static std::array<std::array<scalar, 2>, 2> view_range(
    const contour_container& cc_) {

    std::array<scalar, 2> x_range = __e_object._x_range;
    std::array<scalar, 2> y_range = __e_object._y_range;

    for (auto& c : cc_) {
        for (auto& v : c) {
            x_range[0] = std::min(v[0], x_range[0]);
            x_range[1] = std::max(v[0], x_range[1]);
            y_range[0] = std::min(v[1], y_range[0]);
            y_range[1] = std::max(v[1], y_range[1]);
        }
    }
    return {x_range, y_range};
}

/** Helper method to connect the surface sheets to the 
 * surfaces of the layer_sheets
 * 
 * @tparam volume_type the type of volume (templated on point3_container)
 * 
 * @param v_ the input volume
 * @param templates_ the given module templtes
 * @param o_ the object to which they are attached
 * 
 **/
template <typename volume_type>
void connect_surface_sheets(const volume_type& v_,
                            std::vector<svg::object>& templates_,
                            svg::object& o_) {
    // Now create an item per surface
    for (auto [is, s] : utils::enumerate(v_._surfaces)) {
        std::string sid = s._name;

        // Template copy
        size_t it = v_._templates[is];
        if (it >= v_._templates.size()){
            continue;
        }
        svg::object s_sheet_s = templates_[v_._templates[is]];

        s_sheet_s._attribute_map["display"] = "none";

        // Get some text
        scalar yl = s_sheet_s._y_range[1];
        auto surface_info = draw::text({0,yl}, "info_"+s._name, s._info);
        s_sheet_s.add_object(surface_info);

        // Object information to appear
        svg::object on;
        on._tag = "animate";
        on._attribute_map["fill"] = "freeze";
        on._attribute_map["attributeName"] = "display";
        on._attribute_map["from"] = "none";
        on._attribute_map["to"] = "block";
        on._attribute_map["begin"] = sid + __d + "mouseout";

        svg::object off;

        off._tag = "animate";
        off._attribute_map["fill"] = "freeze";
        off._attribute_map["attributeName"] = "display";
        off._attribute_map["to"] = "none";
        off._attribute_map["from"] = "block";
        off._attribute_map["begin"] = sid + __d + "mouseover";

        
        // Store the animation
        s_sheet_s._sub_objects.push_back(off);
        s_sheet_s._sub_objects.push_back(on);

        o_.add_object(s_sheet_s);
    }
}

}  // namespace display

}  // namespace actsvg
