message(STATUS "Building 'actsvg::meta' component")

# Set up the library.
actsvg_add_library( actsvg_meta meta
   "include/actsvg/display/datamodel.hpp"
   "include/actsvg/display/geometry.hpp"
   "include/actsvg/display/grids.hpp"
   "include/actsvg/display/helpers.hpp"
   "include/actsvg/display/materials.hpp"
   "include/actsvg/display/sheets.hpp"
   "include/actsvg/proto/cluster.hpp"
   "include/actsvg/proto/detector.hpp"
   "include/actsvg/proto/grid.hpp"
   "include/actsvg/proto/material.hpp"
   "include/actsvg/proto/portal.hpp"
   "include/actsvg/proto/surface.hpp"
   "include/actsvg/proto/track.hpp"
   "include/actsvg/proto/volume.hpp"
   "include/actsvg/styles/defaults.hpp"
   "include/actsvg/meta.hpp"

   "src/display/geometry.cpp"
   "src/display/grids.cpp"
   "src/display/helpers.cpp"
   "src/display/materials.cpp"
   "src/proto/material.cpp"
   "src/styles/defaults.cpp"
)

target_link_libraries( actsvg_meta PUBLIC actsvg::core)
