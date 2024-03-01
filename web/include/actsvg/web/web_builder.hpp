// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "actsvg/core/svg.hpp"
#include "actsvg/web/webpage_text.hpp"
#include <string>
#include <fstream>
#include <filesystem>
#include <functional>
#include <numeric>


namespace actsvg::web {

namespace compare{
// Function to compare two strings alphanumerically.
bool alphanumeric(const actsvg::svg::object& svg1, const actsvg::svg::object& svg2){
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
} // namespace comparison_functions

/// @brief Class for generating a web page to view and merge svgs.
/// @note When used a debugging tool and rebuilding mutliple times,
/// if the webpage does not refresh as expected it is likely caused by browser caching.
class web_builder{
    public:

    using comparison_function = bool (*)(const actsvg::svg::object&, const actsvg::svg::object&);

    /// @brief Generates a webpage configured with the given svgs.
    /// @param output_directory the root of the web page.
    /// @param svgs the svgs avaible for selection on the web page.
    /// @param order_comparator a compartor function to determine 
    /// the display order of the svgs.
    /// @note When used a debugging tool and rebuilding mutliple times,
    /// if the webpage does not refresh as expected it is likely caused by browser caching.
    template <typename iterator_t>
    void build(const std::filesystem::path& output_directory, const iterator_t& svgs, comparison_function order_comparator)
    {
        create_template(output_directory);
        configure_svgs(output_directory, svgs, order_comparator);
    }

    /// @brief Generates a webpage configured with the given svgs.
    /// @param output_directory the root of the web page.
    /// @param svgs the svgs avaible for selection on the web page.
    /// @note When used a debugging tool and rebuilding mutliple times,
    /// if the webpage does not refresh as expected it is likely caused by browser caching.
    template <typename iterator_t>
    void build(const std::filesystem::path& output_directory, const iterator_t& svgs)
    {
        create_template(output_directory);
        configure_svgs(output_directory, svgs);
    }

    private:

    void create_template(const std::filesystem::path& output_directory)
    {
        create_directory_structure(output_directory);
        write_file(output_directory / "index.html", index_text);
        write_file(output_directory / "script.js", script_text);
        write_file(output_directory / "styles.css", css_text);
        write_file(output_directory / "rebuild.py", rebuild_text);
    }

    /// @brief Writes the svgs to the svgs folder in the output directory and generates the config.json file.
    template <typename iterator_t>
    void configure_svgs(const std::filesystem::path& output_directory, const iterator_t& it_svgs){
        std::vector<std::string> file_names;
        // Write the svgs.
        for (const auto& svg : it_svgs)
        {
            svg::file output_file;
            output_file.add_object(svg);

            const std::string file_name = svg._id + ".svg";

            write_file(output_directory / "svgs" / file_name, output_file);

            file_names.push_back(file_name);
        }

        // Write the config.
        write_file(output_directory / "config.json", get_config(file_names));
    }

    /// @brief Writes the svgs to the svgs folder in the output directory and generates the config.json file.
    template <typename iterator_t>
    void configure_svgs(const std::filesystem::path& output_directory, const iterator_t& it_svgs, comparison_function order_comparator){
        auto svgs = std::vector(it_svgs.cbegin(), it_svgs.cend());
        std::sort(svgs.begin(), svgs.end(), order_comparator);
        configure_svgs(output_directory, svgs);
    }

    /// @brief Writes a file.
    template <typename content_t>
    void write_file(const std::filesystem::path& path, const content_t& content)
    {
        std::ofstream rstream;
        rstream.open(path);
        rstream << content;
        rstream.close();
    }

    /// @returns the json of the config file.
    template <typename iterator_t>
    std::string get_config(const iterator_t& file_names_it)
    {
        std::string res;
        for (const auto& file_name : file_names_it){
            if (res.size() == 0){
                res += "\"" + file_name + "\"";
            }
            else{
                res += ", \"" + file_name + "\"";
            }
        }
        return "[" + res + "]";
    }

    /// @brief Creates the folder structure for the webpage.
    void create_directory_structure(const std::filesystem::path& directory_path){
        if (!std::filesystem::exists(directory_path)) {
            std::filesystem::create_directory(directory_path);
        }
        const auto svg_directory = directory_path / "svgs";
        if (!std::filesystem::exists(svg_directory)) {
            std::filesystem::create_directory(svg_directory);
        }
    }
};

}  // namespace actsvg::web
