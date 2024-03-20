message(STATUS "Building 'actsvg::python' bindings")


set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/actsvg)

set(_python_install_dir "python/actsvg")

set(modules_src 
  src/core.cpp
  src/display.cpp
  src/io.cpp
  src/modules.cpp  
  src/proto.cpp
  src/style.cpp
)

if (ACTSVG_BUILD_EXAMPLES)
  list(APPEND modules_src
    src/examples.cpp)
else()
  list(APPEND modules_src
    src/dummies/examples.cpp)
endif()
  

pybind11_add_module(pyactsvg 
  ${modules_src}
)


target_link_libraries(pyactsvg PUBLIC actsvg::core actsvg::meta)
if (ACTSVG_BUILD_EXAMPLES)
  target_link_libraries(pyactsvg PUBLIC actsvg::data)
endif()

install(TARGETS pyactsvg DESTINATION ${_python_install_dir})

configure_file(setup.sh.in ${_python_dir}/setup.sh @ONLY)
install(FILES ${_python_dir}/setup.sh DESTINATION "python")

set_target_properties(pyactsvg PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}")
set_target_properties(pyactsvg PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/actsvg)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/actsvg/ DESTINATION ${_python_install_dir})
