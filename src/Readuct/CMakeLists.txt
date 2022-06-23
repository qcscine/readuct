# Load lists of source files
include(Files.cmake)

# Find/Import dependecies
include(ImportUtilsOS)
include(ImportCore)
import_utils_os()
import_core()
include(ImportYamlCpp)
import_yamlcpp()

if(NOT TARGET Boost::filesystem OR NOT TARGET Boost::program_options)
  find_package(Boost REQUIRED COMPONENTS filesystem system program_options)
endif()

function(target_set_rpath)
  set(options "")
  set(oneValueArgs "")
  set(multiValueArgs "TARGETS;UNIX_RPATH")
  cmake_parse_arguments(SCINE_RPATH_TARGETS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(APPLE)
    string(REPLACE $ORIGIN @loader_path osx_rpath "${SCINE_RPATH_UNIX_RPATH}")
    set_target_properties(${SCINE_RPATH_TARGETS} PROPERTIES
      BUILD_WITH_INSTALL_RPATH OFF
      INSTALL_RPATH "${osx_rpath}"
    )
  elseif(UNIX)
    set_target_properties(${SCINE_RPATH_TARGETS} PROPERTIES
      BUILD_WITH_INSTALL_RPATH OFF
      INSTALL_RPATH "${SCINE_RPATH_UNIX_RPATH}"
    )
  endif()
endfunction()

# Executable
add_executable(ReaductApp ${READUCT_APP_FILES})
add_executable(Scine::ReaductApp ALIAS ReaductApp)
set_target_properties(ReaductApp PROPERTIES OUTPUT_NAME readuct)

# Shared Readuct library / Readuct module
add_library(Readuct ${READUCT_ELEMENTARY_STEP_OPTIMIZATION_FILES})
if(NOT BUILD_SHARED_LIBS)
  set_target_properties(Readuct PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
target_include_directories(Readuct
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
  PRIVATE
    ${YAML_CPP_INCLUDE_DIR}
)
set_target_properties(Readuct PROPERTIES OUTPUT_NAME readuct)
target_link_libraries(Readuct PUBLIC Scine::UtilsOS PRIVATE yaml-cpp)
target_compile_options(Readuct PUBLIC
  $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
)
add_library(Scine::Readuct ALIAS Readuct)
install(
  TARGETS Readuct
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Tests
if(SCINE_BUILD_TESTS)
  include(ImportGTest)
  import_gtest()
  add_executable(Readuct_tests ${READUCT_TEST_FILES})
  target_link_libraries(Readuct_tests PRIVATE
    GTest::Main
    GMock::GMock
    Scine::Readuct
    yaml-cpp
  )
  target_include_directories(Readuct_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/App/
    ${YAML_CPP_INCLUDE_DIR}
  )
  add_test(
    NAME Readuct
    COMMAND Readuct_tests --gtest_output=xml:${CMAKE_CURRENT_BINARY_DIR}/readuct_binary_report.xml
  )
  target_compile_options(Readuct_tests PUBLIC
    $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
  )
endif()

target_include_directories(ReaductApp PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/App/
  ${YAML_CPP_INCLUDE_DIR}
)
target_link_libraries(ReaductApp PRIVATE
  Scine::UtilsOS
  Scine::Core
  Scine::Readuct
  Boost::boost
  Boost::filesystem
  Boost::program_options
  yaml-cpp
)
target_compile_options(ReaductApp PUBLIC
  $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
)

target_set_rpath(
  TARGETS Readuct ReaductApp
  UNIX_RPATH "\$ORIGIN;\$ORIGIN/../lib"
)

# Targets
install(TARGETS ReaductApp EXPORT reaductTargets RUNTIME DESTINATION bin)

scine_install_component_cmake_files(
  COMPONENT ScineReaduct
  EXPORT_NAME reaductTargets
)

# Python bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  include(ImportPybind11)
  import_pybind11()
  set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

  # Python module
  pybind11_add_module(scine_readuct ${READUCT_PYTHON_FILES})
  set_target_properties(scine_readuct PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scine_readuct
  )
  target_include_directories(scine_readuct PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/App
    ${YAML_CPP_INCLUDE_DIR}
  )
  target_link_libraries(scine_readuct PRIVATE
    Scine::UtilsOS
    Scine::Core
    Scine::Readuct
    Boost::boost
    Boost::filesystem
    yaml-cpp
  )
  target_set_rpath(TARGETS scine_readuct UNIX_RPATH "\$ORIGIN")
  install(CODE
    "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade --no-deps ${CMAKE_CURRENT_BINARY_DIR}
                     RESULT_VARIABLE retcode)
     if(NOT \$\{retcode\} EQUAL 0)
       message(FATAL_ERROR \"Fatal error when installing Python module using PIP.\")
     endif()"
  )

  # Copy the setup file and README.rst into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python/pkginit.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_readuct
  )
  file(
    RENAME ${CMAKE_CURRENT_BINARY_DIR}/scine_readuct/pkginit.py
    ${CMAKE_CURRENT_BINARY_DIR}/scine_readuct/__init__.py
  )

  # Figure out which targets we have to copy
  include(TargetLibName)
  set(_py_targets_to_copy Scine::Core)  # Core is always shared
  target_lib_type(Scine::UtilsOS _utils_libtype)
  if(_utils_libtype STREQUAL "SHARED_LIBRARY")
    list(APPEND _py_targets_to_copy Scine::UtilsOS)
  endif()
  unset(_utils_libtype)
  if(BUILD_SHARED_LIBS)
    list(APPEND _py_targets_to_copy Readuct)
  endif()

  # Generate generator expressions for each target and figure out filenames
  # for the python setup file
  set(readuct_PY_DEPS "")
  foreach(target ${_py_targets_to_copy})
    list(APPEND _py_target_gen_exprs "\$<TARGET_FILE:${target}>")
    target_lib_filename(${target} _target_filename)
    string(APPEND readuct_PY_DEPS ", \"${_target_filename}\"")
  endforeach()
  message(STATUS "Targets to copy for python bindings: ${_py_targets_to_copy}")
  unset(_py_targets_to_copy)

  add_custom_command(TARGET Readuct POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${_py_target_gen_exprs} ${CMAKE_CURRENT_BINARY_DIR}/scine_readuct
    COMMENT "Copying required shared libraries into python package directory"
  )
  unset(_py_target_gen_exprs)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
    ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    @ONLY
  )
  unset(readuct_PY_DEPS)

  file(
    COPY ${PROJECT_SOURCE_DIR}/README.rst
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )

  include(SphinxDocumentation)
  scine_sphinx_documentation(
    TARGET scine_readuct
    CONFIGURATION ${CMAKE_CURRENT_SOURCE_DIR}/Python/sphinx/conf.py.in
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Python/sphinx
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doc-py
    LINK UtilsPythonDoc
    DOCTEST
    DOCTEST_REQUIRES scine_utilities
  )
endif()

# Setup App tests if possible
if((SCINE_BUILD_TESTS) AND (TARGET Scine::Sparrow))
  string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
  if(build_type MATCHES debug)
    add_test(
      NAME Readuct_Application
      COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/Tests/AppTests/test_readuct.py -k TestReaductFast --junitxml=${CMAKE_CURRENT_BINARY_DIR}/readuct_application_report.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  else()
    add_test(
      NAME Readuct_Application
      COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/Tests/AppTests/test_readuct.py --junitxml=${CMAKE_CURRENT_BINARY_DIR}/readuct_application_report.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()

  include(TargetLibName)
  target_lib_directory(Scine::Sparrow _sparrow_module_path)
  message(STATUS "Sparrow module found at ${_sparrow_module_path}")
  target_lib_directory(Scine::UtilsOS _utils_module_path)
  message(STATUS "Utils module found at ${_utils_module_path}")
  set_tests_properties(Readuct_Application PROPERTIES
    ENVIRONMENT "PATH=$ENV{PATH}:${CMAKE_CURRENT_BINARY_DIR};SCINE_MODULE_PATH=${_sparrow_module_path}:${_utils_module_path}"
  )
  unset(_sparrow_module_path)
  unset(_utils_module_path)
endif()

if(SCINE_BUILD_PYTHON_BINDINGS)
  file(COPY ${PROJECT_SOURCE_DIR}/LICENSE.txt  ${EXTERNAL_LICENSE_FILES}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_readuct
  )
  file(COPY ${EXTERNAL_LICENSE_FILES}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/external_licenses
  )
endif()
