# This file is part of the actsvg packge.
#
# Copyright (C) 2022 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import sys
import datetime


current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
js_path = os.path.join(current_directory, "template", "script.js")
html_path = os.path.join(current_directory, "template", "index.html")
css_path = os.path.join(current_directory, "template", "styles.css")
rebuild_path = os.path.join(current_directory, "template", "rebuild.py")

namespace = "actsvg::web"
file_text_variables = [("index_text", html_path), ("script_text", js_path), ("css_text", css_path), ("rebuild_text", rebuild_path)]

current_year = str(datetime.datetime.now().year)
license = r"""// This file is part of the actsvg package.
//
// Copyright (C) """ + current_year + r""" CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""

file_content = license + "#pragma once\n\n#include <string>\n\n" + "namespace " + namespace + " {\n\n"
for variable_name, path in file_text_variables:
    with open(path, "r") as file:
        file_content += "const std::string " + variable_name + " = R\"(" + file.read() + ")\";\n\n"
file_content += "} //  namespace " + namespace

outputFile = os.path.join(current_directory,"webpage_text.hpp")
with open(outputFile, "w") as file:
    file.write(file_content)
print('Header File Created')
