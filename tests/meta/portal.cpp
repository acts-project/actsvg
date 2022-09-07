// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/proto/portal.hpp"

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "actsvg/core/defs.hpp"

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;

TEST(proto, portal) {

    // Create and define a volume
    proto::portal<point3_container> p;

    ASSERT_TRUE(p._name == "unnamed");
    ASSERT_TRUE(p._volume_links.empty());
}
