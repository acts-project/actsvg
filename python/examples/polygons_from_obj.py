#!/usr/bin/env python

import argparse
import actsvg
from actsvg import io, core, display


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Script to read polygon surface objects from an obj files"
    )
    p.add_argument("--obj", type=str, default="", help="OBJ file")
    p.add_argument("--mtl", type=str, default="", help="MTL file")
    p.add_argument("--view", type=str, default="xy", help="View type (xy, rz)")
    p.add_argument("--range", type=float, nargs="+", help="Range restriction")

    args = p.parse_args()

    polygons = io.obj.read_polygons([args.obj], [args.mtl])
    print("** pyactsvg **: read", len(polygons), "polygons from", args.obj)

    if args.range and args.view == "xy":
        polygons = display.select_view_range(polygons, "z", args.range)

    surfaces_xy = display.surfaces(polygons, args.view)

    s_xy_file = actsvg.io.file()
    s_xy_file.add_objects(surfaces_xy)
    s_xy_file.write("surfaces_xy.svg")
