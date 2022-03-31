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
#include "actsvg/proto/detector.hpp"

using namespace actsvg;

using point3 = std::array<scalar,3>;
using point3_container = std::vector<point3>;

TEST(proto, detector) {

    // Create and define a volume
    proto::detector<point3_container> d;
    
    ASSERT_TRUE(d._name == "unnamed");
    ASSERT_TRUE(d._volumes.empty());
    
}
