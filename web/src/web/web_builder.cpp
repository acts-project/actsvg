// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "actsvg/web/web_builder.hpp"

#include <filesystem>
#include <fstream>
#include <functional>
#include <numeric>
#include <string>

#include "actsvg/core/svg.hpp"
#include "actsvg/web/webpage_text.hpp"

namespace actsvg::web {

namespace compare {
bool alphanumeric(const actsvg::svg::object& svg1,
                  const actsvg::svg::object& svg2) {
    const std::string& str1 = svg1._id;
    const std::string& str2 = svg2._id;
    // Iterate through both strings simultaneously
    size_t i = 0, j = 0;
    while (i < str1.length() && j < str2.length()) {
        if (std::isdigit(str1[i]) && std::isdigit(str2[j])) {
            // If both characters are digits, compare them as integers
            int num1 = 0, num2 = 0;
            while (i < str1.length() && std::isdigit(str1[i])) {
                num1 = num1 * 10 + (str1[i] - '0');
                i++;
            }
            while (j < str2.length() && std::isdigit(str2[j])) {
                num2 = num2 * 10 + (str2[j] - '0');
                j++;
            }
            if (num1 != num2) {
                return num1 < num2;
            }
        } else {
            // Compare characters as alphabetic
            if (str1[i] != str2[j]) {
                return str1[i] < str2[j];
            }
            i++;
            j++;
        }
    }

    // Handle the case where one string is a prefix of the other
    return i == str1.length() && j != str2.length();
}
}  // namespace compare

void web_builder::create_template(
    const std::filesystem::path& output_directory) {
    create_directory_structure(output_directory);
    write_file(output_directory / "index.html", index_text);
    write_file(output_directory / "script.js", script_text);
    write_file(output_directory / "styles.css", css_text);
    write_file(output_directory / "rebuild.py", rebuild_text);
}

void web_builder::create_directory_structure(
    const std::filesystem::path& directory_path) {
    if (!std::filesystem::exists(directory_path)) {
        std::filesystem::create_directory(directory_path);
    }
    const auto svg_directory = directory_path / "svgs";
    if (!std::filesystem::exists(svg_directory)) {
        std::filesystem::create_directory(svg_directory);
    }
}  // namespace actsvg::web
}  // namespace actsvg::web
