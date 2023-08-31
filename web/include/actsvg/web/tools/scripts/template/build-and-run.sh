#!/bin/bash

# Change to the directory of the Bash script
cd "$(dirname "$0")"

# Run the build script
python3 rebuild.py

# Start the Python 3 HTTP server
python3 -m http.server
