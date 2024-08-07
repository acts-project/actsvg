import json as jsn
import actsvg
from actsvg import proto
from actsvg import style

""" read a grid axis from a json file"""


def read_grid_axis(json_grid_axis):
    edges = []
    if json_grid_axis["type"] == 0:
        n_bins = json_grid_axis["bins"]
        a_range = json_grid_axis["range"]
        for i in range(n_bins + 1):
            edges.append(a_range[0] + i * (a_range[1] - a_range[0]) / n_bins)
    return edges


""" read a grid from a json file"""


def read_grid(json_grid):
    grid = proto.grid()
    try:
        grid_axes = json_grid["axes"]
        grid.edges_0 = read_grid_axis(grid_axes[0])
        if len(grid_axes) > 1:
            grid.edges_1 = read_grid_axis(grid_axes[1])
        index_data = json_grid["data"]

    except KeyError:
        print("** pyactsvg **: `read_grid` problem in key access ")
    return grid, index_data


""" read surface mateiral from a json file"""


def read_surface_material(json_surface_material):
    psm = proto.surface_material()
    try:
        surface_material = json_surface_material["value"]["material"]
        surface_material_accessor = surface_material["accessor"]
        # Convert the grid
        accessor_grid = surface_material_accessor["grid"]
        psm.grid, index_data = read_grid(accessor_grid)

        # indicate the grid type
        global_casts = surface_material["global_casts"]
        if global_casts == ["binZ", "binPhi"]:
            psm.grid.type = proto.grid.grid_type.z_phi
            psm.grid.reference_r = 100
        elif global_casts == ["binR", "binPhi"]:
            psm.grid.type = proto.grid.grid_type.r_phi
        elif global_casts == ["binX", "binY"]:
            psm.grid.type = proto.grid.grid_type.x_y

        # Create the matrix material  matrix
        nbins0 = len(psm.grid.edges_0) - 1
        nbins1 = 1
        if len(psm.grid.edges_1) > 0:
            nbins1 = len(psm.grid.edges_1) - 1
        psm.material_matrix = [
            [proto.material_slab() for i in range(nbins0)] for j in range(nbins1)
        ]

        surface_material_vector_json = surface_material_accessor["material"]
        surface_material_vector = []
        for smj in surface_material_vector_json:
            if smj["material"] == None:
                sm = proto.material_slab()
            else:
                sm = proto.material_slab(smj["material"], smj["thickness"])
            surface_material_vector.append(sm)
        psm = proto.surface_material.fill_indexed_material(
            psm, surface_material_vector, index_data
        )
        psm.evaluate_material_ranges()

    except KeyError:
        print("** pyactsvg **: `read_surface_material` problem in key access")
        return psm

    return psm


""" read surface material maps from a json file """


def read_surface_material_maps(json_surface_material_maps):
    proto_maps = []
    try:
        for smap in json_surface_material_maps["Surfaces"]["entries"]:
            psm = read_surface_material(smap)
            proto_maps.append(psm)
    except KeyError:
        print(
            "** pyactsvg **: `read_surface_material_maps(...)` raised key access error"
        )
    return proto_maps


""" read a surface from a json file"""


def read_surface(json_surface, apply_transform=True):

    # Construct the name
    surface_name = (
        "module_vol"
        + str(json_surface["volume"])
        + "_lay"
        + str(json_surface["layer"])
        + "_sen"
        + str(json_surface["sensitive"])
    )

    # Read the surface description

    surface_description = json_surface["value"]
    surface_translation = (0.0, 0.0, 0.0)
    surface_rotation = (
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
    )

    if apply_transform:
        surface_transform = surface_description["transform"]
        surface_translation = surface_transform["translation"]
        if surface_transform["rotation"] is not None:
            surface_rotation = surface_transform["rotation"]
            surface_rotation = (
                (surface_rotation[0], surface_rotation[1], surface_rotation[2]),
                (surface_rotation[3], surface_rotation[4], surface_rotation[5]),
                (surface_rotation[6], surface_rotation[7], surface_rotation[8]),
            )

    surface_vertices = []
    bounds_type = surface_description["bounds"]["type"]
    bounds_values = surface_description["bounds"]["values"]
    if bounds_type == "TrapezoidBounds":
        hx_miny = bounds_values[0]
        hx_maxy = bounds_values[1]
        hy = bounds_values[2]
        surface_vertices = [
            (-hx_miny, -hy, 0.0),
            (hx_miny, -hy, 0.0),
            (hx_maxy, hy, 0.0),
            (-hx_maxy, hy, 0.0),
        ]
        # Trapezoid translates into a polygon
        ps = proto.surface.polygon_from_vertices_and_transform(
            surface_name,
            surface_vertices,
            surface_translation,
            surface_rotation,
            style.defaults.sensitive_fill(),
            style.defaults.sensitive_stroke(),
        )
    elif bounds_type == "RectangleBounds":
        hx = bounds_values[0]
        hy = bounds_values[1]
        surface_vertices = [
            (-hx, -hy, 0.0),
            (hx, -hy, 0.0),
            (hx, hy, 0.0),
            (-hx, hy, 0.0),
        ]
        # Rectangle translates into a polygon
        ps = proto.surface.polygon_from_vertices_and_transform(
            surface_name,
            surface_vertices,
            surface_translation,
            surface_rotation,
            style.defaults.sensitive_fill(),
            style.defaults.sensitive_stroke(),
        )
    elif bounds_type == "AnnulusBounds":
        ps = proto.surface.annulus_from_bounds_and_transform(
            surface_name,
            bounds_values,
            surface_translation,
            surface_rotation,
            style.defaults.sensitive_fill(),
            style.defaults.sensitive_stroke(),
        )
    else:
        print("** pyactsvg **: `json.read_surface` bounds type not (yet) supported")
        raise ValueError

    return ps
