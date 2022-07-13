// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <iostream>
#include <strstream>

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
    sstream << std::setw(pr_) << s_;
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
    sstream << std::setw(pr_) << "(" << p_[0] << "," << p_[1] << ")";
    return sstream.str();
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
};

}  // namespace utils

}  // namespace actsvg