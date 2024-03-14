import json as jsn

from .pyactsvg import *
from . import pyactsvg

from actsvg import io, proto

class json :
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
            psm.grid = json.read_grid(accessor_grid)
        except KeyError:
            print("** pyactsvg **: `read_surface_material` problem in key access")
            return psm
        return psm
        
    """ read surface material maps from a json file """
    @staticmethod
    def read_surface_material_maps( filename ) :
        proto_maps = []
        # read the file from disk
        with open(filename, 'r') as f:
            smaps = jsn.load(f)
            try :
                for smap in smaps["Surfaces"]["entries"]  :
                    psm = json.read_surface_material(smap)
                    proto_maps.append(psm)
            except KeyError :
                print('** pyactsvg **: no "Surfaces/entries" in', filename)
        return proto_maps
    
setattr(io, 'json', json)