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

/** Helper from id to url
 * @param id_ the idnetification to be transformed
 **/
static inline std::string id_to_url(const std::string &id_) {
    return std::string("url(#") + id_ + ")";
}

/** Helper to format point2
 *
 * @param s_ the scalar
 * @param pr_ the precision
 *
 * @return a string
 */
static inline std::string to_string(const scalar &s_, size_t pr_ = 4) {
    std::stringstream sstream;
    sstream << std::setprecision(pr_) << s_;
    return sstream.str();
}

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

/** Helper method toscala a point3 object
 *
 * @param p0_ the point3 object and @param s_ the scale
 *
 * @return a new point3
 **/
template <typename point3_type>
point3_type scale(const point3_type &p0_, scalar s_) {
    point3_type scaled;
    scaled[0] = s_ * p0_[0];
    scaled[1] = s_ * p0_[1];
    scaled[2] = s_ * p0_[2];
    return scaled;
}

/** Helper method to add two point3 objects
 *
 * @param p0_ and @param p1_ the 3D points
 *
 * @return a new point3
 **/
template <typename point3_type>
point3_type add(const point3_type &p0_, const point3_type &p1_) {
    point3_type added;
    added[0] = p0_[0] + p1_[0];
    added[1] = p0_[1] + p1_[1];
    added[2] = p0_[2] + p1_[2];
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
template <typename point3_collection, typename point3_type>
point3_collection place_vertices(const point3_collection &pc_,
                                 const point3_type &tr_,
                                 const std::array<point3_type, 3> &rt_) {
    point3_collection placed;
    for (const auto &p : pc_) {
        placed.push_back(place(p, tr_, rt_));
    }
    return placed;
}

}  // namespace utils

}  // namespace actsvg