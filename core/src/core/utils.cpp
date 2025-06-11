// This file is part of the actsvg package.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/core/utils.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <tuple>

namespace actsvg {

namespace utils {
std::string id_to_url(const std::string &id_) {
    return std::string("url(#") + id_ + ")";
}

std::string to_string(const scalar &s_, size_t pr_, value_format f_) {
    std::stringstream sstream;
    if (f_ == value_format::e_scientific) {
        sstream << std::scientific;
    }
    sstream << std::setprecision(pr_) << s_;
    return sstream.str();
}
}  // namespace utils

}  // namespace actsvg
