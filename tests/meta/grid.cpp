// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include "actsvg/proto/grid.hpp"

using namespace actsvg;

TEST(proto, grid) {

    proto::grid g;
    ASSERT_TRUE(g._type == proto::grid::type::e_r_phi);

}
