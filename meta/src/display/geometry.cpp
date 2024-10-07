// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/display/geometry.hpp"

#include <algorithm>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/proto/surface.hpp"

namespace actsvg {

using namespace defaults;

namespace display {
svg::object eta_lines(
    const std::string& id_, scalar zr_, scalar rr_,
    const std::vector<std::tuple<std::vector<scalar>, style::stroke, bool,
                                 style::font>>& els_,
    const style::transform& tr_) {

    svg::object e = svg::object::create_group(id_);
    e._transform = tr_;

    auto theta_from_eta = [](scalar eta) -> scalar {
        return static_cast<scalar>(2. * std::atan(std::exp(-eta)));
    };

    scalar theta_cut = std::atan2(rr_, zr_);

    for (auto [iet, elt] : utils::enumerate(els_)) {
        auto stroke = std::get<style::stroke>(elt);
        for (auto [ie, eta] :
             utils::enumerate(std::get<std::vector<scalar>>(elt))) {
            scalar theta = theta_from_eta(eta);
            std::array<scalar, 2> start = {0., 0.};
            std::array<scalar, 2> end;
            if (theta < theta_cut) {
                end = {zr_, static_cast<scalar>(zr_ * std::tan(theta))};
            } else if (theta > M_PI - theta_cut) {
                end = {-zr_, static_cast<scalar>(-zr_ * std::tan(theta))};
            } else {
                end = {static_cast<scalar>(rr_ * 1 / std::tan(theta)), rr_};
            }
            // Draw the line
            std::string uid = std::to_string(iet) + "_" + std::to_string(ie);
            auto e_line =
                draw::line(id_ + "eta_line_" + uid, start, end, stroke);
            e.add_object(e_line);
            // Label it if told to do so
            if (std::get<bool>(elt)) {
                auto font = std::get<style::font>(elt);
                end[0] +=
                    static_cast<scalar>(std::cos(theta) * 0.5 * font._size);
                end[1] +=
                    static_cast<scalar>(std::sin(theta) * 0.5 * font._size);
                if (eta == 0.) {
                    end[0] -= static_cast<scalar>(0.5 * font._size);
                }
                const std::string eta_code = "\u03B7";
                auto e_text = eta_code + " = " + utils::to_string(eta);
                auto e_label =
                    draw::text(id_ + "eta_label_" + uid, end, {e_text}, font);
                e.add_object(e_label);
            }
        }
    }
    return e;
}

}  // namespace display

}  // namespace actsvg
