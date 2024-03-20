#!/usr/bin/env python

import argparse
import actsvg
import json

from actsvg import io, display, proto, style

if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Script to read surface mapts from a json files"
    )
    p.add_argument("--input", type=str, default="", help="JSON file")

    args = p.parse_args()

    with open(args.input, "r") as f:
        surface_material_maps = json.load(f)
        maps = io.json.read_surface_material_maps(surface_material_maps)
        print("** pyactsvg **: read", len(maps), "material maps from", args.input)

        for i, m in enumerate(maps):
            # Lets make a gradient
            sm_displayed = display.surface_material("surface_material_" + str(i), m)
            sm_file = actsvg.io.file()
            sm_file.add_object(sm_displayed)
            sm_file.write("surface_material_" + str(i) + ".svg")
