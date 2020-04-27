# Load lists of source files
include(Files.cmake)

# Find/Import dependecies
include(ImportUtilsOS)
include(ImportCore)
import_utils_os()
import_core()
find_package(Eigen3 3.3.1 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS filesystem system program_options)
set(Boost_USE_STATIC_LIBS OFF)

# Set the RPATH to be used when installing. $ORIGIN refers to the installation
# directory of the binary.
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

# Executable
set(Boost_USE_STATIC_LIBS OFF)

add_executable(ReaductApp ${READUCT_APP_FILES})
add_executable(Scine::ReaductApp ALIAS ReaductApp)
set_target_properties(ReaductApp PROPERTIES
  OUTPUT_NAME readuct
)
target_include_directories(ReaductApp
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/App/
)

# Shared Readuct library / Readuct module
add_library(Readuct SHARED ${READUCT_ELEMENTARY_STEP_OPTIMIZATION_FILES})
set_property(TARGET Readuct PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(Readuct PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
set_target_properties(Readuct PROPERTIES OUTPUT_NAME Readuct)
target_link_libraries(Readuct PUBLIC Scine::Core Scine::UtilsOS)
target_compile_definitions(Readuct PUBLIC BOOST_ALL_DYN_LINK)
target_compile_options(Readuct
  PUBLIC
  $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
)
add_library(Scine::Readuct ALIAS Readuct)
install(TARGETS Readuct LIBRARY DESTINATION lib)

# Tests
if(SCINE_BUILD_TESTS)
  add_executable(Readuct_tests ${READUCT_TEST_FILES})
  target_link_libraries(Readuct_tests
    PRIVATE
    gtest_main
    gmock
    Scine::Readuct
  )
  add_test(NAME Readuct COMMAND Readuct_tests)
  target_compile_options(Readuct_tests
    PUBLIC
    $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
  )
endif()

target_link_libraries(ReaductApp
 PRIVATE
  Scine::UtilsOS
  Scine::Core
  Scine::Readuct
  Boost::boost
  Boost::filesystem
  Boost::program_options
)
target_compile_options(ReaductApp
  PUBLIC
    $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
)

target_compile_definitions(ReaductApp
                           PUBLIC
                           BOOST_ALL_DYN_LINK
                           )

# Targets
install(TARGETS ReaductApp EXPORT reaductTargets RUNTIME DESTINATION bin)

# Python bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  include(ImportPybind11)
  import_pybind11()
  set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

  # Python module
  pybind11_add_module(scine_readuct ${READUCT_PYTHON_FILES})
  set_target_properties(scine_readuct PROPERTIES SUFFIX ".so")
  set_target_properties(scine_readuct PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scine_readuct)
  target_include_directories(scine_readuct PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/App)
  target_link_libraries(scine_readuct
    PRIVATE
      Scine::UtilsOS
      Scine::Core
      Scine::Readuct
      Boost::boost
      Boost::filesystem
  )
  target_compile_options(scine_readuct
    PUBLIC
      $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
  )
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade ${CMAKE_CURRENT_BINARY_DIR})")

  # Python tests
  # TODO: Uncomment once tests are available
  #if(SCINE_BUILD_TESTS)
  #  add_test(
  #    NAME ReaDuct_Python_Bindings
  #    COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/Python/Tests --junitxml=${CMAKE_CURRENT_BINARY_DIR}/pytest_report.xml
  #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  #  )
  #  set_tests_properties(ReaDuct_Python_Bindings PROPERTIES
  #    ENVIRONMENT SCINE_MODULE_PATH=${CMAKE_CURRENT_BINARY_DIR}
  #    ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}
  #  )
  #endif()

  # Copy the setup file and README.md into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python/__init__.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_readuct
  )
  file(
    COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
      ${PROJECT_SOURCE_DIR}/README.md
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
