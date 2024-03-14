import json as jsn
import actsvg
from actsvg import proto

""" read a grid axis from a json file"""
@staticmethod
def read_grid_axis( json_grid_axis ) :
    edges = []
    
    return edges

""" read a grid from a json file"""
@staticmethod
def read_grid( json_grid ) :
    grid = proto.grid()
    try:
        grid_axes = json_grid["axes"]
        
    except KeyError:
        print("** pyactsvg **: `read_grid` problem in key access ")
    return grid

""" read surface mateiral from a json file"""
@staticmethod
def read_surface_material( json_surface_material) :
    psm = proto.surface_material()
    try:
        surface_material = json_surface_material["value"]["material"]
        surface_material_accessor = surface_material["accessor"]
        # Convert the grid
        accessor_grid = surface_material_accessor["grid"]
        psm.grid = read_grid(accessor_grid)
    except KeyError:
        print("** pyactsvg **: `read_surface_material` problem in key access")
        return psm
    return psm
    
""" read surface material maps from a json file """
@staticmethod
def read_surface_material_maps( json_surface_material_maps ) :
    proto_maps = []
    try:
        for smap in json_surface_material_maps["Surfaces"]["entries"]  :
            psm = read_surface_material(smap)
            proto_maps.append(psm)
    except KeyError :
            print('** pyactsvg **: `read_surface_material_maps` problem in key access')
    return proto_maps

""" read a surface from a json file"""
@staticmethod
def read_surface( json_surface ) :
    ps = proto.surface()
    return ps

    
