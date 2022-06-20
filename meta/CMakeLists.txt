message(STATUS "Building 'actsvg::meta' component")

# Set up the library.
actsvg_add_library( actsvg_meta meta
   "include/actsvg/proto/detector.hpp"
   "include/actsvg/proto/volume.hpp"
   "include/actsvg/proto/surface.hpp"
   "include/actsvg/proto/grid.hpp"
   "include/actsvg/display/datamodel.hpp"
   "include/actsvg/display/helpers.hpp"
   "include/actsvg/display/sheets.hpp"
   "include/actsvg/styles/defaults.hpp")

target_link_libraries( actsvg_meta INTERFACE actsvg::core)

