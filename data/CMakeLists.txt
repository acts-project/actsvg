message(STATUS "Building 'actsvg::data' component")

# Set up the library.
actsvg_add_library( actsvg_data data
   "include/actsvg/data/odd_pixel_barrel.hpp"
   "include/actsvg/data/odd_pixel_ec.hpp")
