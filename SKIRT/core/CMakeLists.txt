# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for the SKIRT core layer
# ------------------------------------------------------------------

# set the target name
set(TARGET skirtcore)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} schema fundamentals)
include_directories(../../SMILE/schema ../../SMILE/fundamentals)

# add SKIRT library dependencies
target_link_libraries(${TARGET} fitsio voro mpi utils gas)
include_directories(../fitsio ../voro ../mpi ../utils ../gas)

# add gasmodule dependency
set(GASMODULE_INSTALL_DIR "$ENV{HOME}/Repositories/gasmodule/cmake_release")
include_directories("${GASMODULE_INSTALL_DIR}/include")
find_library(GASMODULE_LIB
  NAMES libgasmodule.a
  PATHS "${GASMODULE_INSTALL_DIR}/lib"
  NO_DEFAULT_PATH)
target_link_libraries(${TARGET} ${GASMODULE_LIB})

# add GSL dependencies
find_package(GSL REQUIRED)
target_link_libraries(${TARGET} GSL::gsl GSL::gslcblas)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
