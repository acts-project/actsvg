// This file is part of the actsvg package.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../utilities.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

namespace actsvg {
namespace python {

/// @brief Add the core module to the
/// @param ctx
void add_examples_module(context& /*unused*/) {}
}  // namespace python
}  // namespace actsvg
