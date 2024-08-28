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

#include "utilities.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

namespace actsvg::python {

/// The modules, for each
/// @param ctx is the python context
void add_core_module(context& ctx);
void add_io_module(context& ctx);
void add_style_module(context& ctx);
void add_proto_module(context& ctx);
void add_display_module(context& ctx);
void add_examples_module(context& ctx);
}  // namespace actsvg::python

PYBIND11_MODULE(pyactsvg, m) {

    actsvg::python::context ctx;
    ctx.modules["main"] = m;
    // These modules are always present
    add_core_module(ctx);
    add_io_module(ctx);
    add_style_module(ctx);
    add_proto_module(ctx);
    add_display_module(ctx);
    // Conditionally: could be dummy module
    add_examples_module(ctx);
}
