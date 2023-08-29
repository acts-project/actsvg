// Initialize the page by loading a form with a button for each SVG available.
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
