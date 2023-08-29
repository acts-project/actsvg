// This file is part of the actsvg packge.
//
// Copyright (C) 2022 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

namespace actsvg::web {

const std::string index_text = R"(<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>SVG Viewer</title>
    <link rel="stylesheet" type="text/css" href="styles.css">
</head>
<body>
    <div class="app">
        <h1>Select SVGs</h1>
        <div id="formContainer" class="svg-form">
            <!-- The checkboxes will be added here -->
        </div>

        <div id="result" class="result">
            <!-- The resulting SVG will be added here -->
        </div>
    </div>

    <script src="script.js"></script>
</body>
</html>

)";

const std::string script_text = R"(// Initialize the page by loading a form with a button for each SVG available.
// This is done by reading the paths in the config.json file.
document.addEventListener("DOMContentLoaded", function () {
    fetch('./config.json').then(response => response.json()).then(json =>{
        load_form(json);
    });
});

// Loads a form containing a button for each SVG available.
function load_form(group)
{
    let form = document.createElement("form");
    form.id = "checkboxForm";

    group.forEach(item =>{

        itemDiv = document.createElement("div");
        itemDiv.classList.add("form-item");

        let label = document.createElement("label");
        itemDiv.appendChild(label);

        let checkbox = document.createElement("input");
        checkbox.type = "checkbox";
        checkbox.name = "checkbox";
        checkbox.value = "./svgs/" + item;
        label.append(checkbox);

        let span = document.createElement("span");
        span.textContent = get_name(item);
        label.appendChild(span);

        form.appendChild(itemDiv);

    });

    form.append(document.createElement("br"));

    // Create a button to apply the changes
    let apply_button = document.createElement("input");
    apply_button.type = "button";
    apply_button.value = "Apply Selection";
    apply_button.onclick = getSelectedValues;
    apply_button.classList.add("apply-button");

    // Append the apply button to the form
    form.append(apply_button);

    // Append the form to the formContainer div
    const formContainer = document.getElementById("formContainer");
    formContainer.appendChild(form);
}

// For removing the file extension.
function get_name(path){
    return path.replace(".svg", "");
}

// Updates the displayed SVG by combining the selected SVGs.
async function getSelectedValues() {
    var checkboxes = document.getElementsByName("checkbox");
    var selectedValues = [];

    for (var i = 0; i < checkboxes.length; i++) {
        if (checkboxes[i].checked) {
            selectedValues.push(checkboxes[i].value);
        }
    }
    
    selectedValues = selectedValues.reverse();
    let svg = await combineSVGS(selectedValues);
    var resultDiv = document.getElementById('result');
    resultDiv.innerHTML = svg;
}

// Removes the <svg> and </svg> tag to obtain its content.
function removeSVGTag(data){
    let startTag = /<svg[^>]*>/i;
    let endTag = /<\/svg>/i;
    data = data.replace(startTag, '');
    data = data.replace(endTag, '');
    return data;
}

// Given a collection of paths, returns a collection of containing the respective file text.
async function readFiles(paths){
    let result = []
    for (const path of paths) {
        try{
            let data = await (await fetch(path)).text();
            result.push(data);
        }
        catch(err) 
        {
            console.error(err);
        }
    }
    return result;
}

// Creates an SVG from a list of svg objects by joining and appending svg start and end tags.
function createSVG(contents){
    let result = []
    result.push('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="-300 -300 600 600">\n');
    contents.forEach(e => {
        result.push(e);
        result.push('\n');
    });
    result.push('</svg>');
    return result.join('');
}

// Given a list of paths to svgs, returns an svg of their combined content.
async function combineSVGS(paths){
    let contents = await readFiles(paths);
    contents = contents.map(removeSVGTag);
    return createSVG(contents);
}
)";

const std::string css_text = R"(.app{
    font-family: monospace;
    text-align: center;
}

.app h1{
    font-size: 40px;
}

.svg-form{
    text-align: center;
    align-items: center;
}

.form-item{
    display: inline-block;
    border-color: white;
    background-color: transparent;
    font-size: 16px;
    font-weight: 1000;
    margin: 8px;
    max-width: 1000px;
    height: 20px;
}

.form-item label {
    border-radius: 10px;
    user-select: none;
    cursor: pointer;
}

.form-item label span {
    text-align: center;
    padding: 3px 3px;
    display: block;
    border-radius: 10px;
    border: grey;
    border: solid;
    border-width: 2px;
}

.form-item label input {
    display: none;
}

.form-item input:checked + span{
    background: linear-gradient(
        45deg,
        rgba(24, 67, 90, 0.8),
        rgba(88, 12, 31, 0.8) 100%
      );
  color: white;
}

.apply-button{
    border: 10px;
    margin: 20px;
    padding: 10px;
    font-family: monospace;
}

.result{
    box-shadow: 0 4px 17px rgba(0, 0, 0, 0.35);
    width: 60%;
    height: 60%;
    display: inline-block;
}
)";

} //  namespace actsvg::web