// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include "actsvg/core/draw.hpp"

using namespace actsvg;

TEST(text, unconnected_text) {

    svg::object t = draw::text("t0", {10, 10}, {"text"});

    std::cout << t << std::endl;
}
