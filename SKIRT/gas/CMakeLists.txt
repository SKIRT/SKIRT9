# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# All the code which directly accesses the gas module is built as a
# separate library. This allows users who are not interested in gas,
# to build the code without having to download build the gas module.
# ------------------------------------------------------------------

# set the target name
set(TARGET gas)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SKIRT library dependencies
target_link_libraries(${TARGET} utils mpi)
include_directories(../utils ../mpi)

# add SMILE library dependencies
target_link_libraries(${TARGET} fundamentals)
include_directories(../../SMILE/fundamentals)

# define a user-configurable option to build with gas support, which
# requires the compiled gas module to be present on the system
option(BUILD_WITH_GAS "build with gas support - requires gas module installed at GASMODULE_INSTALL_DIR")

# if requested, find gas module installation and configure accordingly
if (BUILD_WITH_GAS)
  # find and configure gas module.
  if (NOT DEFINED GASMODULE_INSTALL_DIR)
    message("Please set GASMODULE_INSTALL_DIR - trying /usr/local/ by default")
    set(GASMODULE_INSTALL_DIR "/usr/local")
  endif()
  message("Looking for gas module in ${GASMODULE_INSTALL_DIR}")
  find_library(GASMODULE_LIBRARIES
    NAMES libgasmodule.a
    PATHS "${GASMODULE_INSTALL_DIR}/lib")
  find_path(GASMODULE_INCLUDE_PATH
    NAMES "GasInterface.hpp"
    PATHS "${GASMODULE_INSTALL_DIR}/include")
  target_link_libraries(${TARGET} ${GASMODULE_LIBRARIES})
  include_directories(${GASMODULE_INCLUDE_PATH})
  # add the C++ preprocessor symbol that triggers compilation with gas
  add_definitions(-DBUILD_WITH_GAS)

  # add GSL dependencies
  find_package(GSL REQUIRED)
  target_link_libraries(${TARGET} GSL::gsl GSL::gslcblas)
endif()

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
