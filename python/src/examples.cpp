// This file is part of the actsvg packge.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "actsvg/data/odd_pixel_barrel.hpp"
#include "actsvg/data/odd_pixel_endcap.hpp"
#include "utilities.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

namespace actsvg {
namespace python {

/// @brief Add the core module to the
/// @param ctx
void add_examples_module(context& ctx) {

    auto& m = ctx.get("main");
    auto e = m.def_submodule("examples");

    // Generate modules ressembling a barrel or endcap
    {
        e.def("generate_barrel_modules",
              []() { return data::generate_barrel_modules(); });

        e.def("generate_endcap_modules",
              []() { return data::generate_endcap_modules(); });
    }
}
}  // namespace python
}  // namespace actsvg