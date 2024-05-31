// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"
#include "actsvg/styles/defaults.hpp"

namespace actsvg {

using namespace defaults;

namespace display {
void prepare_axes(std::array<scalar, 2>& first_, std::array<scalar, 2>& second_,
                  scalar sx_, scalar sy_, scalar ax_, scalar ay_) {
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
}  // namespace display

}  // namespace actsvg
