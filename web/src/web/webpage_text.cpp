// This file is part of the actsvg package.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string>

namespace actsvg::web {

std::string index_text = R"(<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta http-equiv="Cache-Control" content="no-cache, no-store, must-revalidate" />
    <meta http-equiv="Pragma" content="no-cache" />
    <meta http-equiv="Expires" content="0" />
    <title>SVG Viewer</title>
    <link rel="stylesheet" type="text/css" href="styles.css">
</head>
<body>
    <div class="app">
        <h1>Select SVGs</h1>
        <div id="formContainer" class="svg-form">
            <!-- The checkboxes will be added here -->
        </div>

        <div id="result-div" class="result-div" onmouseover="document.body.style.overflow='hidden';"  onmouseout="document.body.style.overflow='auto';">
            <svg id="result-svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="-500 -500 1000 1000">
            
            </svg>
        </div>
        <br>
        <div class="tooltip">&#9432;
            <span>Click and drag to move around. Use the mouse wheel to zoom.</span>
        </div>
        
    </div>
  


    <script src="script.js"></script>
</body>
</html>

)";

std::string script_text =
    R"(let SVGResult = document.getElementById('result-svg');
let formContainer = document.getElementById('formContainer');
let SVGContainer = document.getElementById("result-div");

// Given a collection of paths, returns a collection of containing the respective file text.
async function readFiles(paths){
    let result = []
    for (const path of paths) {
        let data = await (await fetch(path)).text();
        result.push(data);
    }
    return result;
}

// Removes the <svg> and </svg> tag to obtain its content.
function removeSVGTag(data){
    let startTag = /<svg[^>]*>/i;
    let endTag = /<\/svg>/i;
    data = data.replace(startTag, '');
    data = data.replace(endTag, '');
    return data;
}

async function SVGContentMerge(paths){
    let contents = await readFiles(paths);
    return contents.map(removeSVGTag).join('\n');
}

// Loads a form containing a button for each SVG available.
function load_form(group)
{
    let form = document.createElement('form');
    form.id = 'checkboxForm';

    group.forEach(item =>{

        itemDiv = document.createElement('div');
        itemDiv.classList.add('form-item');

        let label = document.createElement("label");
        itemDiv.appendChild(label);

        let checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.name = 'checkbox';
        checkbox.value = './svgs/' + item;
        label.append(checkbox);

        let span = document.createElement('span');
        span.textContent = item.replace('.svg', '');
        label.appendChild(span);

        form.appendChild(itemDiv);

    });

    form.append(document.createElement('br'));

    // Create a button to apply the changes
    let apply_button = document.createElement('input');
    apply_button.type = 'button';
    apply_button.value = 'Apply Selection';
    apply_button.onclick = applyChanges;
    apply_button.classList.add('apply-button');

    // Append the apply button to the form
    form.append(apply_button);

    // Append the form to the formContainer div
    formContainer.appendChild(form);
}

// Initialize the page by loading a form with a button for each SVG available.
// This is done by reading the paths in the config.json file.
document.addEventListener('DOMContentLoaded', function () {
    fetch('./config.json').then(response => response.json()).then(json =>{
        load_form(json);
    });
});

// Updates the displayed SVG by combining the selected SVGs.
async function applyChanges() {
    var checkboxes = document.getElementsByName('checkbox');
    var selectedValues = [];

    for (var i = 0; i < checkboxes.length; i++) {
        if (checkboxes[i].checked) {
            selectedValues.push(checkboxes[i].value);
        }
    }
    
    paths = selectedValues.reverse();
    SVGResult.innerHTML = await SVGContentMerge(paths);
}


// For navigation on the svg:

//--- Adjustable paramters:
// Speed of zoom on scroll.
const zoomFactor = 1.1;
// To calculate the maximum width and height of the viewbox.
const maxHalfLengths = { x: 3000, y: 3000 };
//----

let pivot = { x: 0, y: 0 };
let position = { x: 0, y: 0 };
let halfSize = { x: 500, y: 500 };
const clamp = (num, min, max) => Math.min(Math.max(num, min), max);

setViewBox();

SVGResult.addEventListener("mousemove", onMousemove);
SVGResult.addEventListener("mousedown", onMousedown);
SVGResult.addEventListener("wheel", onWheel);

const mouse = {
    position: { x: 0, y: 0 },
    isDown: false,
};

function viewBoxDim(){
    return {x: position.x - halfSize.x, y: position.y - halfSize.y, w: 2*halfSize.x, h: 2*halfSize.y};
}

function onMousedown(e) {
    mouse.position = screenToViewBoxPosition(e.pageX, e.pageY);
    window.addEventListener("mouseup", onMouseup);
    mouse.isDown = true;
}

function setViewBox() {
    let dim = viewBoxDim();
    SVGResult.setAttribute("viewBox", `${dim.x} ${dim.y} ${dim.w} ${dim.h}`);
}

function screenToViewBoxPosition(screenX, screenY){
    return {
        x: screenX * 2 * halfSize.x/SVGResult.clientWidth,
        y: screenY * 2 * halfSize.y/SVGResult.clientHeight
    }
}

function onMousemove(e) {
    if (mouse.isDown) {
        let pos = screenToViewBoxPosition(e.pageX, e.pageY);
        let dx = (pos.x - mouse.position.x);
        let dy = (pos.y - mouse.position.y);
        position = {x: clamp(position.x - dx, -maxHalfLengths.x, maxHalfLengths.x), y: clamp(position.y - dy, -maxHalfLengths.y, maxHalfLengths.y)};
        mouse.position = pos;
        setViewBox();
    }
}

function onMouseup(e) {
    window.removeEventListener("mouseup", onMouseup);
    mouse.isDown = false;
}

function onWheel(e) {
    const scale = e.deltaY > 0 ? zoomFactor : 1/zoomFactor;
    halfSize = {x: clamp(halfSize.x * scale, 1, maxHalfLengths.x), y: clamp(halfSize.y * scale, 1, maxHalfLengths.y)};
    setViewBox();
}
)";

std::string css_text = R"(.app{
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

.result-div{
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.35);
    height: 40%;
    width: 40%;
    display: inline-block;
    cursor: grab;
}

.result-div svg g{
    cursor: pointer;
}

.tooltip {
  font-size: 16px;
  position: relative;
  display: inline-block;
  border-bottom: 1px dotted black;
}

.tooltip span {
  width: 300px;
  visibility: hidden;
  background-color: black;
  color: #fff;
  text-align: center;
  border-radius: 2px;
  padding: 8px;
  position: absolute;
  z-index: 1;
}

.tooltip:hover span {
  visibility: visible;
}
)";

std::string rebuild_text = R"(import os
import sys
import re

# Sets up the config.json.
# This allows the website script to know the available SVGs.

svg_relative_path = 'svgs'
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
svg_directory = os.path.join(current_directory, svg_relative_path)

if os.path.exists(svg_directory) and os.path.isdir(svg_directory):
    # Get the file names.
    file_names = os.listdir(svg_directory)
    print("Found SVG directory with " + str(len(file_names)) + " files")
    
    # Sort the file names alphanumerically.
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    file_names.sort(key=alphanum_key)
    
    # Write the json file used by the webpage.
    outputFile = os.path.join(current_directory,"config.json")
    with open(outputFile, "w") as file:
        content = "[" + ", ".join(["\"" + f + "\"" for f in file_names]) + "]"
        file.write(content)
    print('Build Complete')
else:
    print("Build Error: Directory Not Found")
)";

}  //  namespace actsvg::web
