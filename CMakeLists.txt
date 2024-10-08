cmake_minimum_required(VERSION 3.9)

# Set the project name. The chosen project name must be in CamelCase. Your CMake
# tree must then provide a properly namespaced target with the same name as
# your project.
project(Readuct
  VERSION 6.0.0
  DESCRIPTION "This is the SCINE module ReaDuct."
)

# Set the module path for universal cmake files inclusion
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake)

# Component setup
include(ComponentSetup)
scine_setup_component()

# Don't build tests of dependencies
set(_build_tests ${SCINE_BUILD_TESTS})
set(SCINE_BUILD_TESTS OFF)
option(BUILD_SPARROW "Will download and build Sparrow (the SCINE semi-empirical module)." OFF)
if(BUILD_SPARROW)
  include(ImportSparrow)
  import_sparrow()
endif()

option(BUILD_XTB "Will download and build SCINE Xtb (the SCINE wrapper around the xtb API by the Grimme group)." OFF)
if(BUILD_XTB)
  include(ImportXtb)
  import_xtb()
endif()

option(BUILD_SWOOSE "Will download and build Swoose (the SCINE QM/MM module)." OFF)
if(BUILD_SWOOSE)
  include(ImportSwoose)
  import_swoose()
endif()
set(SCINE_BUILD_TESTS ${_build_tests})

# Enable testing
if(SCINE_BUILD_TESTS)
  enable_testing()
endif()

option(SCINE_USE_MKL "Use the optimized MKL library for linear algebra operations of Eigen" OFF)

# Subdirectories
add_subdirectory(src)

# Documentation
scine_component_documentation(UtilsOSDocumentation)
