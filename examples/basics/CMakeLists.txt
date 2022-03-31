message(STATUS "Building 'actsvg::examples'")

actsvg_add_executable(basic_shapes shapes.cpp
                               LINK_LIBRARIES actsvg::core Boost::program_options)
