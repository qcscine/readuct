# Load lists of source files
include(Files.cmake)

# Find/Import dependecies
include(ImportUtilsOS)
include(ImportCore)
import_utils_os()
import_core()
find_package(Eigen3 3.3.1 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS filesystem system)

# Set the RPATH to be used when installing. $ORIGIN refers to the installation directory of the binary.
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

# Executable
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS program_options REQUIRED)

add_executable(ReaductApp ${READUCT_APP_FILES})
add_executable(Scine::ReaductApp ALIAS ReaductApp)
set_target_properties(ReaductApp PROPERTIES
  OUTPUT_NAME readuct
)
target_include_directories(ReaductApp 
  PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/App/
)
target_link_libraries(ReaductApp
 PRIVATE
  Scine::UtilsOS
  Scine::Core
  Boost::boost
  Boost::filesystem
)

# Targets
install(TARGETS ReaductApp 
        EXPORT reaductTargets
        RUNTIME DESTINATION bin
)

# Python bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  include(ImportPybind11)
  import_pybind11()
  set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})
  
  # Python module
  pybind11_add_module(scine_readuct ${READUCT_PYTHON_FILES})
  set_target_properties(scine_readuct PROPERTIES SUFFIX ".so")
  target_include_directories(scine_readuct 
    PRIVATE 
      ${CMAKE_CURRENT_SOURCE_DIR}/App
  )
  target_link_libraries(scine_readuct
    PRIVATE
      Scine::UtilsOS
      Scine::Core
      Boost::boost
      Boost::filesystem
  )
  install(TARGETS scine_readuct LIBRARY DESTINATION lib)
  
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
    COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
      ${PROJECT_SOURCE_DIR}/README.md
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
endif(SCINE_BUILD_PYTHON_BINDINGS)
