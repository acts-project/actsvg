// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <pybind11/pybind11.h>

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "actsvg/core.hpp"

namespace actsvg {
namespace python {

using point3 = std::array<scalar, 3u>;
using point3_collection = std::vector<point3>;

struct context {
    std::unordered_map<std::string, pybind11::module_> modules;

    pybind11::module_& get(const std::string& name) { return modules.at(name); }

    template <typename... Args,
              typename = std::enable_if_t<sizeof...(Args) >= 2>>
    auto get(Args&&... args) {
        return std::make_tuple((modules.at(args))...);
    }
};

/// Translate to int
struct to_int {

    using return_type = int;

    return_type operator()(const std::string& s) const {
        return std::atoi(s.c_str());
    }
};

/// Translate to scalar
struct to_scalar {

    using return_type = scalar;

    return_type operator()(const std::string& s) const {
        return std::atof(s.c_str());
    }
};

/// Translate to string
struct to_string {

    using return_type = std::string;

    return_type operator()(const std::string& s) const { return s; }
};

/// @brief split an input string
///
/// @param source the source string
/// @param sc the split caracter
/// @param ignore indices e.g. to split indicator off
///
/// @return a vector of valut_type
template <typename translator_type>
std::vector<typename translator_type::return_type> split_string(
    const std::string& source, char sc, std::vector<unsigned int> ignore = {}) {

    // Temporary values
    std::string temp;
    std::stringstream stringstream{source};
    std::vector<typename translator_type::return_type> result;

    translator_type tt;

    unsigned int split = 0;
    while (std::getline(stringstream, temp, sc)) {
        if (std::find(ignore.begin(), ignore.end(), split) == ignore.end()) {
            result.push_back(tt(temp));
        }
        ++split;
    }

    return result;
}

}  // namespace python
}  // namespace actsvg
