message(STATUS "Building 'actsvg::core' component")

# Set up the library.
actsvg_add_library( actsvg_core core
   "include/actsvg/actsvg.hpp"
   "include/actsvg/core/connectors.hpp"
   "include/actsvg/core/defs.hpp"
   "include/actsvg/core/draw.hpp"
   "include/actsvg/core/generators.hpp"
   "include/actsvg/core/style.hpp"
   "include/actsvg/core/svg.hpp"
   "include/actsvg/core/utils.hpp"
   "include/actsvg/core/views.hpp")

