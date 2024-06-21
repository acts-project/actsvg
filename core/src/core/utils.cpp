// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <tuple>

#include "actsvg/core/defs.hpp"

namespace actsvg {

namespace utils {
std::string id_to_url(const std::string &id_) {
    return std::string("url(#") + id_ + ")";
}

std::string to_string(const scalar &s_, size_t pr_) {
    std::stringstream sstream;
    sstream << std::setprecision(pr_) << s_;
    return sstream.str();
}
}  // namespace utils

}  // namespace actsvg
