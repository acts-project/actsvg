// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

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

}  // namespace utils

}  // namespace actsvg