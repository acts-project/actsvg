cmake_minimum_required(VERSION 3.11)

project (actsvg VERSION 0.1 LANGUAGES CXX )

# Set up the used C++ standard(s).
set( CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "Disable C++ extensions" )

# CMake include(s).
include( CMakeDependentOption )
include( GNUInstallDirs )

# Explicitly set the output directory for the binaries. Such that if this
# project is included by another project, the main project's configuration would
# win out.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}" CACHE PATH
   "Directory for the built binaries" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
   "Directory for the built libraries" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
   "Directory for the built static libraries" )

# Include the Detray CMake code.
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" )
include( actsvg-functions )

# Core library component
add_subdirectory(core)

option(ACTSVG_BUILD_META "Build the meta level interface of ACTSVG" TRUE )
if ( ACTSVG_BUILD_META )
   add_subdirectory(meta)
endif()

option(ACTSVG_BUILD_EXAMPLES "Build the example applications of ACTSVG" TRUE )
if ( ACTSVG_BUILD_EXAMPLES )
   # Setup boost for the examples 
   set(_actsvg_boost_version 1.71.0)
   option(ACTSVG_SETUP_BOOST "Explicitly set up Boost for the project" ON)
   option(ACTSVG_USE_SYSTEM_BOOST "Use a system-provided boost" ON)
   if (ACTSVG_SETUP_BOOST)
      if (ACTSVG_USE_SYSTEM_BOOST)
        find_package(Boost ${_actsvg_boost_version} REQUIRED COMPONENTS program_options)
      else()
           add_subdirectory(extern/boost)
      endif()
 endif()
 
   # Add the examples
   add_subdirectory(examples)
endif()

option(ACTSVG_BUILD_TESTING "Build the (unit) tests of ACTSVG" TRUE )

if (ACTSVG_BUILD_TESTING OR ACTSVG_BUILD_EXAMPLES)
   add_subdirectory(data)
endif()

# Set up the test(s).
include( CTest )
if( ACTSVG_BUILD_TESTING )
  add_subdirectory( tests )
endif()

# Set up GoogleTest.
option( ACTSVG_SETUP_GOOGLETEST
   "Set up the GoogleTest target(s) explicitly" TRUE )
option( ACTSVG_USE_SYSTEM_GOOGLETEST
   "Pick up an existing installation of GoogleTest from the build environment"
   FALSE )
if( ACTSVG_SETUP_GOOGLETEST )
   if( ACTSVG_USE_SYSTEM_GOOGLETEST )
      find_package( GTest REQUIRED )
   else()
      add_subdirectory( extern/googletest )
   endif()
endif()

# Set up the packaging of the project.
include( actsvg-packaging )

