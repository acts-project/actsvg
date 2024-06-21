// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <string>

namespace actsvg {
const std::string __l = "<";
const std::string __r = ">\n";
const std::string __el = "</";
const std::string __er = "/>\n";
const std::string __fs = "/";
const std::string __tab = "\t";
const std::string __nl = "\n";
const std::string __blk = " ";
const std::string __c = ",";
const std::string __d = ".";

/// @todo make configurable via compile time
using scalar = float;

using point2 = std::array<scalar, 2>;

}  // namespace actsvg
