message(STATUS "Building 'actsvg::web' component")

# Set up the library.
actsvg_add_library( actsvg_web web
   "include/actsvg/web/web_builder.hpp"
   "include/actsvg/web/webpage_text.hpp"

   "src/web/web_builder.cpp"
   "src/web/webpage_text.cpp"
)

target_link_libraries( actsvg_web PUBLIC actsvg::core)
