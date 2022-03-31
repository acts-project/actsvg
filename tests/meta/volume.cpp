// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "actsvg/core/defs.hpp"
#include "actsvg/proto/volume.hpp"

using namespace actsvg;

using point3 = std::array<scalar,3>;
using point3_container = std::vector<point3>;

TEST(proto, volume) {

    // Create and define a volume
    proto::volume<point3_container> v;
    
    ASSERT_TRUE(v._name == "unnamed");
    ASSERT_TRUE(v._surfaces.empty());
    ASSERT_TRUE(v._portals.empty());
    
}
