cmake_minimum_required(VERSION 3.11)

project (actsvg VERSION 0.4.41 LANGUAGES CXX )

# Set up the used C++ standard(s).
set( CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "Disable C++ extensions" )

set(ACTSVG_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wshadow -Wno-unused-local-typedefs")
# This adds some useful conversion checks like float-to-bool, float-to-int, etc.
# However, at the moment this is only added to clang builds, since GCC's -Wfloat-conversion 
# is much more aggressive and also triggers on e.g., double-to-float
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang") 
 set(ACTSVG_CXX_FLAGS "${ACTSVG_CXX_FLAGS} -Wfloat-conversion")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ACTSVG_CXX_FLAGS}")

# CMake include(s).
include( CMakeDependentOption )
include( GNUInstallDirs )

# Flags controlling the meta-build system.
option( ACTSVG_USE_SYSTEM_LIBS "Use system libraries be default" FALSE )

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

option(ACTSVG_BUILD_META "Build the meta level interface of ACTSVG" ON )
if ( ACTSVG_BUILD_META )
   add_subdirectory(meta)
endif()

option(ACTSVG_BUILD_WEB "Build the webpage builder interface of ACTSVG" ON )
if ( ACTSVG_BUILD_WEB )
   add_subdirectory(web)
endif()


option(ACTSVG_BUILD_TESTING "Build the (unit) tests of ACTSVG" OFF )

if (ACTSVG_BUILD_TESTING OR ACTSVG_BUILD_EXAMPLES)
   add_subdirectory(data)
endif()

# Set up the test(s).
include( CTest )
if( ACTSVG_BUILD_TESTING )
  add_subdirectory( tests )
endif()

if (ACTSVG_BUILD_TESTING OR ACTSVG_BUILD_EXAMPLES)
   # Set up GoogleTest.
   option( ACTSVG_SETUP_GOOGLETEST
      "Set up the GoogleTest target(s) explicitly" TRUE )
   option( ACTSVG_USE_SYSTEM_GOOGLETEST
      "Pick up an existing installation of GoogleTest from the build environment"
      ${ACTSVG_USE_SYSTEM_LIBS} )
   if( ACTSVG_SETUP_GOOGLETEST )
      if( ACTSVG_USE_SYSTEM_GOOGLETEST )
         find_package( GTest REQUIRED )
      else()
         add_subdirectory( extern/googletest )
      endif()
   endif()
endif()

option(ACTSVG_BUILD_PYTHON_BINDINGS "Build the python bindings of ACTSVG" OFF)
if(ACTSVG_BUILD_PYTHON_BINDINGS)
  find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development)
  option(ACTSVG_USE_SYSTEM_PYBIND11 "Pick up an existing installation of pybind11")
  if(ACTSVG_USE_SYSTEM_PYBIND11) 
    find_package(pybind11 CONFIG REQUIRED)
  else()
    add_subdirectory(extern/pybind11)
  endif()
  # Add the python sub directory
  add_subdirectory(python)
endif()


# Set up the packaging of the project.
include( actsvg-packaging )
