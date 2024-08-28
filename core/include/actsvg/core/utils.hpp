// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <tuple>

namespace actsvg {

namespace utils {

/** Helper method to run enumerate(...) with structured binding
 * over STL type containers.
 *
 * @param iterable is a std-like iterable container type
 *
 **/
template <typename container_type,
          typename container_type_iter =
              decltype(std::begin(std::declval<container_type>())),
          typename = decltype(std::end(std::declval<container_type>()))>
constexpr auto enumerate(container_type &&iterable) {
    struct iterator {
        size_t i;
        container_type_iter iter;

        bool operator!=(const iterator &rhs) const { return iter != rhs.iter; }

        /** Increase index and iterator at once */
        void operator++() {
            ++i;
            ++iter;
        }

        /** Tie them together for returning */
        auto operator*() const { return std::tie(i, *iter); }
    };
    struct iterable_wrapper {
        container_type iterable;
        auto begin() { return iterator{0, std::begin(iterable)}; }
        auto end() { return iterator{0, std::end(iterable)}; }
    };
    return iterable_wrapper{std::forward<container_type>(iterable)};
}

/** Helper for perp
 *
 * @param p_ the point
 *
 * @return a scalar of the norm
 **/
template <typename point2_type>
scalar perp(const point2_type &p_) {
    return std::sqrt(p_[0] * p_[0] + p_[1] * p_[1]);
}

/** Helper method for calculating the distance between two points
 *
 * @param x0_ the starting point
 * @param x1_ the ending point
 *
 * @return the distance
 **/
template <typename point2_type>
scalar distance(const point2_type &x0_, const point2_type &x1_) {
    return std::sqrt(std::pow(x1_[0] - x0_[0], 2) +
                     std::pow(x1_[1] - x0_[1], 2));
}

/** Helper method to intersect two lines in 2-d
 * @param x0_ the starting point of line 0
 * @param d0_ the direction of line 0
 * @param x1_ the starting point of line 1
 * @param d1_ the direction of line 1
 *
 * @return the intersection point
 **/
template <typename point2_type>
point2_type intersect(const point2_type &x0_, const point2_type &d0_,
                      const point2_type &x1_, const point2_type &d1_) {
    point2_type p_int;

    struct line {
        scalar a;
        scalar d;

        line(const point2_type &x_, const point2_type &d_) {
            a = d_[1u] / d_[0u];
            d = x_[1u] - a * x_[0u];
        }
    };

    if (d0_[0u] == 0) {

        line l1(x1_, d1_);
        p_int[0u] = x0_[0u];
        p_int[1u] = l1.a * p_int[0u] + l1.d;

    } else if (d1_[0u] == 0) {

        line l0(x0_, d0_);
        p_int[0u] = x1_[0u];
        p_int[1u] = l0.a * p_int[0u] + l0.d;

    } else {

        line l0(x0_, d0_);
        line l1(x1_, d1_);

        p_int[0u] = (l1.d - l0.d) / (l0.a - l1.a);
        p_int[1u] = l0.a * p_int[0u] + l0.d;
    }

    return p_int;
}

/** Helper from id to url
 * @param id_ the idnetification to be transformed
 **/
std::string id_to_url(const std::string &id_);

/** Helper to format point2
 *
 * @param s_ the scalar
 * @param pr_ the precision
 *
 * @return a string
 */
std::string to_string(const scalar &s_, size_t pr_ = 4);

/** Helper to format point2
 *
 * @param p_ the point
 * @param pr_ the precision
 *
 * @return a string
 */
template <typename point2_type>
std::string to_string(const point2_type &p_, size_t pr_ = 4) {
    std::stringstream sstream;
    sstream << std::setprecision(pr_) << "(" << p_[0] << "," << p_[1] << ")";
    return sstream.str();
}

/** Check if a parameter is within range
 *
 * @param p_ parameter to be checked
 * @param range_ the range
 *
 * @note borders are included in this check
 **/
template <typename range_type>
bool inside_range(scalar p_, const range_type &range_) {
    return (range_[0] <= p_ and p_ <= range_[1]);
}

/** Helper method to calculate the barycenter
 *
 * @param vs_ vertices
 *
 * @return a new barycenter
 **/
template <typename point2_type>
point2_type barycenter(const std::vector<point2_type> &vs_) {
    point2_type rv = {0., 0.};
    for (const auto &v : vs_) {
        rv[0] += v[0];
        rv[1] += v[1];
    }
    rv[0] /= vs_.size();
    rv[1] /= vs_.size();
    return rv;
}

/** Helper method to rotate a 2-d point
 * @param p_ the point to be rotated
 * @param a_ the angle alpha
 *
 * @return the rotated point
 **/
template <typename point2_type>
point2_type rotate(const point2_type &p_, scalar a_) {
    point2_type p_rot;
    p_rot[0] = std::cos(a_) * p_[0] - std::sin(a_) * p_[1];
    p_rot[1] = std::sin(a_) * p_[0] + std::cos(a_) * p_[1];
    return p_rot;
}

/** Helper method to scale object
 *
 * @param p0_ the point object and @param s_ the scale
 *
 * @return a new point_type object
 **/
template <typename point_type, std::size_t kDIM = 3u>
point_type scale(const point_type &p0_, scalar s_) {
    point_type scaled;
    scaled[0] = s_ * p0_[0];
    scaled[1] = s_ * p0_[1];
    if constexpr (kDIM == 3u) {
        scaled[2] = s_ * p0_[2];
    }
    return scaled;
}

/** Helper method for making a unit vector
 *
 * @param x0_ the starting point
 * @param x1_ the ending point
 *
 * @return the distance
 **/
template <typename point_type, std::size_t kDIM = 2u>
point_type unit_direction(const point_type &x0_, const point_type &x1_) {
    // 2 dimensional case
    if constexpr (kDIM == 2u) {
        scalar length = std::sqrt(std::pow(x1_[0] - x0_[0], 2) +
                                  std::pow(x1_[1] - x0_[1], 2));
        scalar norml = 1. / length;
        return point_type{norml * (x1_[0] - x0_[0]), norml * (x1_[1] - x0_[1])};
    }
    // 3 dimensional case
    if constexpr (kDIM == 3u) {
        scalar length = std::sqrt(std::pow(x1_[0] - x0_[0], 2) +
                                  std::pow(x1_[1] - x0_[1], 2) +
                                  std::pow(x1_[2] - x0_[2], 2));
        scalar norml = 1. / length;
        return point_type{norml * (x1_[0] - x0_[0]), norml * (x1_[1] - x0_[1]),
                          norml * (x1_[2] - x0_[2])};
    }
    return point_type{};
}

/** Helper method to add two point3 objects
 *
 * @param p0_ and @param p1_ the  points
 *
 * @return a new point_type object
 **/
template <typename point_type, std::size_t kDIM = 3u>
point_type add(const point_type &p0_, const point_type &p1_) {
    point_type added;
    added[0] = p0_[0] + p1_[0];
    added[1] = p0_[1] + p1_[1];
    if constexpr (kDIM == 3u) {
        added[2] = p0_[2] + p1_[2];
    }
    return added;
}

/** Helper method to rotate a 3D point @param p_ under
 * roatation @param rt_
 *
 * @return new rotated point3
 **/
template <typename point3_type>
point3_type rotate(const std::array<point3_type, 3> &rt_,
                   const point3_type &p_) {

    point3_type rotated;
    rotated[0] = rt_[0][0] * p_[0] + rt_[0][1] * p_[1] + rt_[0][2] * p_[2];
    rotated[1] = rt_[1][0] * p_[0] + rt_[1][1] * p_[1] + rt_[1][2] * p_[2];
    rotated[2] = rt_[2][0] * p_[0] + rt_[2][1] * p_[1] + rt_[2][2] * p_[2];
    return rotated;
}

/** Helper method to place a 3D point @param p_ under
 * @param tr_ transform and @param rt_ rotation
 *
 * @return new transform
 **/
template <typename point3_type>
point3_type place(const point3_type &p_, const point3_type &tr_,
                  const std::array<point3_type, 3> &rt_) {
    point3_type placed = rotate(rt_, p_);
    placed[0] = p_[0] + tr_[0];
    placed[1] = p_[1] + tr_[1];
    placed[2] = p_[2] + tr_[2];
    return placed;
}

/** Helper method to place a 3D points @param pc_ under
 * @param tr_ transform and @param rt_ rotation
 *
 * @return new transformed point3 collection
 **/
template <typename point3_container, typename point3_type>
point3_container place_vertices(const point3_container &pc_,
                                const point3_type &tr_,
                                const std::array<point3_type, 3> &rt_) {
    point3_container placed;
    for (const auto &p : pc_) {
        placed.push_back(place(p, tr_, rt_));
    }
    return placed;
}

}  // namespace utils

}  // namespace actsvg
