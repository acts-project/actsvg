let SVGResult = document.getElementById('result-svg');
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
