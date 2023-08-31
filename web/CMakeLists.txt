message(STATUS "Building 'actsvg::web' component")

# Set up the library.
actsvg_add_library( actsvg_web web
   "include/actsvg/web/web_viwer.hpp")
target_link_libraries( actsvg_web INTERFACE actsvg::core)
