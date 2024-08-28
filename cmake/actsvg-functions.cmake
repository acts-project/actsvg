# CMake include(s).
include( CMakeParseArguments )

# Helper function for setting up the actsvg libraries.
#
# Usage: actsvg_add_library( actsvg_core core "header1.cpp"... )
#
function( actsvg_add_library fullname basename )
  cmake_parse_arguments( ARG "" "" "" ${ARGN} )

  # Create the library.
  add_library( ${fullname} SHARED ${ARG_UNPARSED_ARGUMENTS} )

  # Set up how clients should find its headers.
  target_include_directories( ${fullname} PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> )

  set_target_properties(${fullname}
    PROPERTIES
    COMPILE_FLAGS ${ACTSVG_CXX_FLAGS})

  # Make sure that the library is available as "actsvg::${basename}" in every
  # situation.
  set_target_properties( ${fullname} PROPERTIES EXPORT_NAME ${basename} )
  add_library( actsvg::${basename} ALIAS ${fullname} )


  # Set up the installation of the library and its headers.
  install( TARGETS ${fullname}
      EXPORT actsvg-exports
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" )
  install( DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" OPTIONAL )

endfunction( actsvg_add_library )



# Helper function for setting up the actsvg executables.
#
# The actsvg executables are *not* installed with the project, as they are only
# used for testing / benchmarking the code. Clients of actsvg do not need them.
#
# Usage: actsvg_add_executable( foo bar.cpp
#                               LINK_LIBRARIES actsvg::core )
#
function( actsvg_add_executable name )

  # Parse the function's options.
  cmake_parse_arguments( ARG "" "" "LINK_LIBRARIES" ${ARGN} )

  # Create the executable.
  set( exe_name "actsvg_${name}" )
  add_executable( ${exe_name} ${ARG_UNPARSED_ARGUMENTS} )
  if( ARG_LINK_LIBRARIES )
    target_link_libraries( ${exe_name} PRIVATE ${ARG_LINK_LIBRARIES} )
  endif()

endfunction( actsvg_add_executable )


# Helper function for setting up the actsvg tests.
#
# Usage: actsvg_add_test( core source1.cpp source2.cpp
#                         LINK_LIBRARIES actsvg::core )
#
function( actsvg_add_test name )

  # Parse the function's options.
  cmake_parse_arguments( ARG "" "" "LINK_LIBRARIES" ${ARGN} )

  # Create the test executable.
  set( test_exe_name "actsvg_test_${name}" )
  add_executable( ${test_exe_name} ${ARG_UNPARSED_ARGUMENTS} )
  if( ARG_LINK_LIBRARIES )
    target_link_libraries( ${test_exe_name} PRIVATE ${ARG_LINK_LIBRARIES} )
  endif()

  # Run the executable as the test.
  add_test( NAME ${test_exe_name}
      COMMAND ${test_exe_name} )

  # Set all properties for the test.
  set_tests_properties( ${test_exe_name} PROPERTIES
      ENVIRONMENT actsvg_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/data/ )

endfunction( actsvg_add_test )
