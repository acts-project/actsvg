import os
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
