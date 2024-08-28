// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#include <fstream>

#include "actsvg/display/datamodel.hpp"
#include "actsvg/proto/cluster.hpp"

using namespace actsvg;

namespace {

// Poor man's cluster position estimation
//
/// Calculates the @param cluster_ position given a bunch
/// of channels in their @param boundaries_ and stored data
/// used for charge interpolation
///
template <size_t DIM, size_t BOUNDS = DIM>
void cluster_position(
    proto::cluster<DIM>& cluster_,
    const std::array<std::vector<scalar>, BOUNDS>& boundaries_) {

    scalar total_data = 0;

    // straight forward 2D case - cartesian/cartesian & polar/polar
    if constexpr (DIM == 2 and BOUNDS == 2) {
        cluster_._measurement[DIM - 2] = 0.;
        cluster_._measurement[DIM - 1] = 0.;
        for (const auto& c : cluster_._channels) {
            // Total weight
            scalar data = c._data;
            total_data += data;
            scalar channel_p0 = 0.;
            scalar channel_p1 = 0.;
            scalar low_p0 = boundaries_[BOUNDS - 2][c._cid[DIM - 2]];
            scalar high_p0 = boundaries_[BOUNDS - 2][c._cid[DIM - 2] + 1];
            scalar low_p1 = boundaries_[BOUNDS - 1][c._cid[DIM - 1]];
            scalar high_p1 = boundaries_[BOUNDS - 1][c._cid[DIM - 1] + 1];
            channel_p0 = 0.5_scalar * (low_p0 + high_p0);
            channel_p1 = 0.5_scalar * (low_p1 + high_p1);
            cluster_._measurement[DIM - 2] += data * channel_p0;
            cluster_._measurement[DIM - 1] += data * channel_p1;
        }
    }
    // same for 1-d case
    if constexpr (DIM == 1 and BOUNDS == 1) {

        cluster_._measurement[DIM - 1] = 0.;
        for (const auto& c : cluster_._channels) {
            // Total weight
            scalar data = c._data;
            total_data += data;
            scalar channel_p = 0.;
            scalar low_p = boundaries_[BOUNDS - 1][c._cid[DIM - 1]];
            scalar high_p = boundaries_[BOUNDS - 1][c._cid[DIM - 1] + 1];
            channel_p = 0.5_scalar * (low_p + high_p);
            cluster_._measurement[DIM - 1] += data * channel_p;
        }
    }
    // cartesian / fan - in 1D or 3D
    if constexpr (BOUNDS == 3) {

        scalar y_min = boundaries_[BOUNDS - 1][0];
        scalar y_max =
            boundaries_[BOUNDS - 1][boundaries_[BOUNDS - 1].size() - 1];

        std::vector<scalar> slopes;
        for (auto [ix, xl] : utils::enumerate(boundaries_[BOUNDS - 3])) {
            scalar xh = boundaries_[BOUNDS - 2][ix];
            slopes.push_back((xh - xl) / (y_max - y_min));
        }

        // Fill a dimension offset
        unsigned int dim_x_offset = (DIM == 1) ? 1 : 0;

        cluster_._measurement[DIM - 2 + dim_x_offset] = 0.;
        cluster_._measurement[DIM - 1] = 0.;

        for (const auto& c : cluster_._channels) {
            // Total weight
            scalar data = c._data;
            total_data += data;

            scalar low_y = 0.;
            scalar high_y = 0.;
            if constexpr (DIM > 1) {
                low_y = boundaries_[BOUNDS - 1][c._cid[DIM - 1]];
                high_y = boundaries_[BOUNDS - 1][c._cid[DIM - 1] + 1];
            } else {
                low_y = boundaries_[BOUNDS - 1][0];
                high_y = boundaries_[BOUNDS - 1][1];
            }

            scalar y = 0.5_scalar * (low_y + high_y);

            scalar low_x_base =
                boundaries_[BOUNDS - 3][c._cid[DIM - 2 + dim_x_offset]];
            scalar high_x_base =
                boundaries_[BOUNDS - 3][c._cid[DIM - 2 + dim_x_offset] + 1];

            scalar low_x = low_x_base +
                           slopes[c._cid[DIM - 2 + dim_x_offset]] * (y - y_min);
            scalar high_x =
                high_x_base +
                slopes[c._cid[DIM - 2 + dim_x_offset] + 1] * (y - y_min);

            scalar x = 0.5_scalar * (low_x + high_x);

            cluster_._measurement[DIM - 2 + dim_x_offset] += data * x;

            // Protection against DIM=1 case
            if constexpr (DIM > 1) {
                cluster_._measurement[DIM - 1] += data * y;
            }
            cluster_._correlation +=
                data *
                (0.5_scalar * (slopes[c._cid[DIM - 2 + dim_x_offset]] +
                               slopes[c._cid[DIM - 2 + dim_x_offset] + 1]));
        }
    }

    // Re-scale the entire thing
    scalar weight = 1._scalar / total_data;

    for (unsigned int id = 0; id < DIM; ++id) {
        cluster_._measurement[id] *= weight;
    }
    cluster_._correlation *= weight;
}

}  // namespace

TEST(proto, cluster1D) {

    // One way to construct it
    proto::channel<1> c1{{3}, 0.5};

    // Or decide to overwrite
    c1._cid = {4};
    c1._data = 0.6_scalar;

    proto::cluster<1> cl;
    cl._channels = {c1};
}

TEST(proto, cluster2D) {

    proto::channel<2> c2a{{3, 1}, 0.5_scalar};
    proto::channel<2> c2b{{3, 2}, 0.2_scalar};
    proto::cluster<2> cl;
    cl._channels = {c2a, c2b};
}

TEST(display, cluster1D_x_cart) {

    proto::channel<1> c1a{{2}, 0.3_scalar};
    proto::channel<1> c1b{{3}, 0.5_scalar};
    proto::channel<1> c1c{{4}, 0.2_scalar};
    proto::cluster<1> cl;

    cl._channels = {c1a, c1b, c1c};
    cl._coords = {proto::cluster<1>::e_x};
    cl._truth = {-5.};
    cl._variance = {12.};

    std::vector<scalar> x_boundaries = {-100, -75, -50, -25, 0,
                                        25,   50,  75,  100};
    std::vector<scalar> y_boundaries = {-200, 200};

    cluster_position(cl, {x_boundaries});

    style::fill no_fill = style::fill();

    auto segmentation = draw::tiled_cartesian_grid("seg_x", x_boundaries,
                                                   y_boundaries, no_fill);
    auto [cluster, focus] = display::cluster(segmentation, "c1x", cl);

    svg::file cluster_file;
    cluster_file.add_object(segmentation);
    cluster_file.add_object(cluster);

    // Write out the file
    std::ofstream clout;
    clout.open("test_meta_cluster_cart_1D_x.svg");
    clout << cluster_file;
    clout.close();
}

TEST(display, cluster1D_y_cart) {

    proto::channel<1> c1a{{2}, 0.3_scalar};
    proto::channel<1> c1b{{3}, 0.5_scalar};
    proto::channel<1> c1c{{4}, 0.2_scalar};
    proto::cluster<1> cl;

    cl._channels = {c1a, c1b, c1c};
    cl._coords = {proto::cluster<1>::e_y};
    cl._truth = {-8.};
    cl._variance = {12.};

    std::vector<scalar> x_boundaries = {-200, 200};
    std::vector<scalar> y_boundaries = {-100, -75, -50, -25, 0,
                                        25,   50,  75,  100};

    cluster_position(cl, {y_boundaries});

    style::fill no_fill = style::fill();

    auto segmentation = draw::tiled_cartesian_grid("seg_y", x_boundaries,
                                                   y_boundaries, no_fill);
    auto [cluster, focus] = display::cluster(segmentation, "c1y", cl);

    svg::file cluster_file;
    cluster_file.add_object(segmentation);
    cluster_file.add_object(cluster);

    // Write out the file
    std::ofstream clout;
    clout.open("test_meta_cluster_cart_1D_y.svg");
    clout << cluster_file;
    clout.close();
}

TEST(display, cluster2D_cart) {

    proto::channel<2> c2a{{3, 1}, 0.3_scalar};
    proto::channel<2> c2b{{3, 2}, 0.5_scalar};
    proto::channel<2> c2c{{3, 3}, 0.2_scalar};
    proto::channel<2> c2d{{4, 3}, 0.1_scalar};
    proto::cluster<2> cl;

    cl._channels = {c2a, c2b, c2c, c2d};
    cl._coords = {proto::cluster<2>::e_x, proto::cluster<2>::e_y};
    cl._truth = {-35., -30.};
    cl._variance = {6., 10.};

    std::vector<scalar> x_boundaries = {-125, -100, -75, -50, -25, 0,
                                        25,   50,   75,  100, 125};

    std::vector<scalar> y_boundaries = {-100, -75, -50, -25, 0,
                                        25,   50,  75,  100};

    cluster_position(cl, {x_boundaries, y_boundaries});

    style::fill no_fill = style::fill();
    auto segmentation = draw::tiled_cartesian_grid("seg_x_y", x_boundaries,
                                                   y_boundaries, no_fill);

    auto [cluster, focus] = display::cluster(segmentation, "c2", cl);

    svg::file cluster_file;

    cluster_file.add_object(segmentation);
    cluster_file.add_object(cluster);

    // Write out the file
    std::ofstream clout;
    clout.open("test_meta_cluster_cart_2D.svg");
    clout << cluster_file;
    clout.close();
}

TEST(display, cluster1D_polar_r) {

    proto::channel<1> ca{{1}, 0.3_scalar};
    proto::channel<1> cb{{2}, 0.5_scalar};
    proto::channel<1> cc{{3}, 0.2_scalar};

    proto::cluster<1> cl;

    cl._channels = {ca, cb, cc};

    cl._type = proto::cluster<1>::e_polar;
    cl._coords = {proto::cluster<1>::e_r};
    cl._truth = {375.};
    cl._variance = {20.};

    std::vector<scalar> r_boundaries = {100, 200, 300, 400, 500, 600};
    std::vector<scalar> phi_boundaries = {0.25, 0.75};

    cluster_position(cl, {r_boundaries});

    style::fill no_fill = style::fill();
    auto segmentation =
        draw::tiled_polar_grid("seg_r", r_boundaries, phi_boundaries, no_fill);
    auto [cluster, focus] = display::cluster(segmentation, "cl", cl);

    svg::file cluster_file;
    cluster_file.add_object(segmentation);
    cluster_file.add_object(cluster);

    // Write out the file
    std::ofstream clout;
    clout.open("test_meta_cluster_polar_1D_r.svg");
    clout << cluster_file;
    clout.close();
}

TEST(display, cluster1D_polar_phi) {

    proto::channel<1> ca{{4}, 0.3_scalar};
    proto::channel<1> cb{{5}, 0.5_scalar};
    proto::channel<1> cc{{6}, 0.2_scalar};

    proto::cluster<1> cl;

    cl._channels = {ca, cb, cc};

    cl._type = proto::cluster<1>::e_polar;
    cl._coords = {proto::cluster<1>::e_phi};
    cl._truth = {0.51_scalar};
    cl._variance = {0.01_scalar};

    std::vector<scalar> r_boundaries = {100, 600};
    std::vector<scalar> phi_boundaries = {0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
                                          0.55, 0.6, 0.65, 0.7, 0.75};

    cluster_position(cl, {phi_boundaries});

    style::fill no_fill = style::fill();
    auto segmentation = draw::tiled_polar_grid("seg_phi", r_boundaries,
                                               phi_boundaries, no_fill);
    auto [cluster, focus] = display::cluster(segmentation, "cl", cl);

    svg::file cluster_file;
    cluster_file.add_object(segmentation);
    cluster_file.add_object(cluster);

    // Write out the file
    std::ofstream clout;
    clout.open("test_meta_cluster_polar_1D_phi.svg");
    clout << cluster_file;
    clout.close();
}

TEST(display, cluster2D_polar_r_phi) {

    proto::channel<2> ca{{1, 2}, 0.3_scalar};
    proto::channel<2> cb{{2, 2}, 0.5_scalar};
    proto::channel<2> cc{{3, 3}, 0.2_scalar};

    proto::cluster<2> cl;

    cl._channels = {ca, cb, cc};

    cl._type = proto::cluster<2>::e_polar;
    cl._coords = {proto::cluster<2>::e_r, proto::cluster<2>::e_phi};
    cl._truth = {375, 0.45_scalar};
    cl._variance = {20., 0.01_scalar};

    std::vector<scalar> r_boundaries = {100, 200, 300, 400, 500, 600};
    std::vector<scalar> phi_boundaries = {0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
                                          0.55, 0.6, 0.65, 0.7, 0.75};
    cluster_position(cl, {r_boundaries, phi_boundaries});

    style::fill no_fill = style::fill();
    auto segmentation = draw::tiled_polar_grid("seg_r_phi", r_boundaries,
                                               phi_boundaries, no_fill);
    auto [cluster, focus] = display::cluster(segmentation, "cl", cl);

    svg::file cluster_file;
    cluster_file.add_object(segmentation);
    cluster_file.add_object(cluster);

    // Write out the file
    std::ofstream clout;
    clout.open("test_meta_cluster_polar_2D_r_phi.svg");
    clout << cluster_file;
    clout.close();
}

TEST(display, cluster1D_fan_x) {

    proto::channel<1> ca{{1}, 0.3_scalar};
    proto::channel<1> cb{{2}, 0.5_scalar};
    proto::channel<1> cc{{3}, 0.2_scalar};

    proto::cluster<1> cl;

    cl._channels = {ca, cb, cc};

    cl._type = proto::cluster<1>::e_fan;
    cl._coords = {proto::cluster<1>::e_x};
    cl._truth = {-37};
    cl._variance = {20.};

    std::vector<scalar> x_low_boundaries = {-40, -30, -20, -10, 0,
                                            10,  20,  30,  40};
    std::vector<scalar> x_high_boundaries = {-96, -72, -48, -24, 0,
                                             24,  48,  72,  96};
    std::vector<scalar> y_boundaries = {-150, 150};

    cluster_position<1, 3>(cl,
                           {x_low_boundaries, x_high_boundaries, y_boundaries});

    style::fill no_fill = style::fill();
    auto segmentation =
        draw::tiled_fan_grid("seg_fan_x", x_low_boundaries, x_high_boundaries,
                             y_boundaries, no_fill);
    auto [cluster, focus] = display::cluster(segmentation, "cl", cl);

    svg::file cluster_file;
    cluster_file.add_object(segmentation);
    cluster_file.add_object(cluster);

    // Write out the file
    std::ofstream clout;
    clout.open("test_meta_cluster_fan_1D_x.svg");
    clout << cluster_file;
    clout.close();
}

TEST(display, cluster1D_fan_y) {

    proto::channel<1> ca{{1}, 0.3_scalar};
    proto::channel<1> cb{{2}, 0.5_scalar};
    proto::channel<1> cc{{3}, 0.2_scalar};

    proto::cluster<1> cl;

    cl._channels = {ca, cb, cc};

    cl._type = proto::cluster<1>::e_fan;
    cl._coords = {proto::cluster<1>::e_y};
    cl._truth = {-10.};
    cl._variance = {20.};

    std::vector<scalar> x_low_boundaries = {-80, 80};
    std::vector<scalar> x_high_boundaries = {-96, 96};
    std::vector<scalar> y_boundaries = {-150, -100, -50, 0, 50, 100, 150};

    cluster_position<1, 1>(cl, {y_boundaries});

    style::fill no_fill = style::fill();
    auto segmentation =
        draw::tiled_fan_grid("seg_fan_y", x_low_boundaries, x_high_boundaries,
                             y_boundaries, no_fill);
    auto [cluster, focus] = display::cluster(segmentation, "cl", cl);

    svg::file cluster_file;
    cluster_file.add_object(segmentation);
    cluster_file.add_object(cluster);

    // Write out the file
    std::ofstream clout;
    clout.open("test_meta_cluster_fan_1D_y.svg");
    clout << cluster_file;
    clout.close();
}

TEST(display, cluster2D_fan_x_y) {

    proto::channel<2> ca{{1, 2}, 0.3_scalar};
    proto::channel<2> cb{{2, 2}, 0.5_scalar};
    proto::channel<2> cc{{3, 3}, 0.2_scalar};

    proto::cluster<2> cl;

    cl._channels = {ca, cb, cc};

    cl._type = proto::cluster<2>::e_fan;
    cl._coords = {proto::cluster<2>::e_x, proto::cluster<2>::e_y};
    cl._truth = {-25., -12.};
    cl._variance = {5., 25.};

    std::vector<scalar> x_low_boundaries = {-40, -30, -20, -10, 0,
                                            10,  20,  30,  40};
    std::vector<scalar> x_high_boundaries = {-96, -72, -48, -24, 0,
                                             24,  48,  72,  96};
    std::vector<scalar> y_boundaries = {-150, -100, -50, 0, 50, 100, 150};

    cluster_position<2, 3>(cl,
                           {x_low_boundaries, x_high_boundaries, y_boundaries});

    style::fill no_fill = style::fill();
    auto segmentation =
        draw::tiled_fan_grid("seg_fan_x_y", x_low_boundaries, x_high_boundaries,
                             y_boundaries, no_fill);
    auto [cluster, focus] = display::cluster(segmentation, "cl", cl);

    svg::file cluster_file;
    cluster_file.add_object(segmentation);
    cluster_file.add_object(cluster);

    // Write out the file
    std::ofstream clout;
    clout.open("test_meta_cluster_fan_2D_x_y.svg");
    clout << cluster_file;
    clout.close();
}

TEST(display, cluster2D_fan_focus) {

    proto::channel<2> ca{{40, 22}, 0.3_scalar};
    proto::channel<2> cb{{41, 22}, 0.5_scalar};
    proto::channel<2> cc{{42, 23}, 0.2_scalar};

    proto::cluster<2> cl;

    cl._channels = {ca, cb, cc};

    cl._type = proto::cluster<2>::e_polar;
    cl._coords = {proto::cluster<2>::e_x, proto::cluster<2>::e_y};
    cl._truth = {375, 0.22_scalar};
    cl._variance = {0.05_scalar, 0.25_scalar};

    std::vector<scalar> x_low_boundaries = {0.};
    x_low_boundaries.reserve(100);
    std::vector<scalar> x_high_boundaries = {0.};
    x_high_boundaries.reserve(100);
    for (int ix = 1; ix < 51; ++ix) {
        x_low_boundaries.push_back(ix * 3);
        x_low_boundaries.push_back(-ix * 3);
        x_high_boundaries.push_back(ix * 5);
        x_high_boundaries.push_back(-ix * 5);
    }
    std::sort(x_low_boundaries.begin(), x_low_boundaries.end());
    std::sort(x_high_boundaries.begin(), x_high_boundaries.end());

    std::vector<scalar> y_boundaries;
    y_boundaries.reserve(50);
    for (int iy = 0; iy < 51; ++iy) {
        y_boundaries.push_back(-250._scalar + iy * 10._scalar);
    }

    cluster_position<2, 3>(cl,
                           {x_low_boundaries, x_high_boundaries, y_boundaries});

    style::fill no_fill = style::fill();
    auto segmentation = draw::tiled_fan_grid(
        "seg_fan", x_low_boundaries, x_high_boundaries, y_boundaries, no_fill);
    auto [cluster, focus] = display::cluster(segmentation, "cl", cl);

    svg::file cluster_file_nf;
    cluster_file_nf.add_object(segmentation);
    cluster_file_nf.add_object(cluster);

    svg::file cluster_file_f;
    cluster_file_f.add_object(focus);
    cluster_file_f.add_object(cluster);

    // Write out the files
    std::ofstream clout;
    clout.open("test_meta_cluster_fan_no_focus.svg");
    clout << cluster_file_nf;
    clout.close();

    clout.open("test_meta_cluster_fan_focus.svg");
    clout << cluster_file_f;
    clout.close();
}
