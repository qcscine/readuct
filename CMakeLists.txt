cmake_minimum_required(VERSION 3.9)

# Set the project name. The chosen project name must be in CamelCase. Your CMake
# tree must then provide a properly namespaced target with the same name as
# your project.
project(Readuct
  VERSION 1.0.0
  DESCRIPTION "This is the SCINE module Readuct."
)

# Set the module path for universal cmake files inclusion
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Component setup
include(ComponentSetup)
scine_setup_component()

# Enable testing
if(SCINE_BUILD_TESTS)
  enable_testing()
endif()

# Enable testing
option(BUILD_SPARROW "Will download and build Sparrow (the Scine semi-empirice module)." OFF)
if(BUILD_SPARROW)
  include(ImportSparrow)
  import_sparrow()
endif()

option(SCINE_USE_MKL "Use the optimized MKL library for linear algebra operations of Eigen" OFF)

# Subdirectories
add_subdirectory(src)
