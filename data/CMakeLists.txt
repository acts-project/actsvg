message(STATUS "Building 'actsvg::data' component")

# Set up the library.
actsvg_add_library( actsvg_data data
   "include/actsvg/data/odd_pixel_barrel.hpp"
   "include/actsvg/data/odd_pixel_endcap.hpp"

   "src/data/odd_pixel_barrel.cpp"
   "src/data/odd_pixel_endcap.cpp"
)

target_link_libraries( actsvg_data PUBLIC actsvg::core)
