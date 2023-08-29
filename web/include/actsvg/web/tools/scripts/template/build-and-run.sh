#!/bin/bash

# Change to the directory of the Bash script
cd "$(dirname "$0")"

# Run the build.py script
python3 build.py

# Start the Python 3 HTTP server
python3 -m http.server
