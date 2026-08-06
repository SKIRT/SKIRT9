# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with support for MPI if requested and available,
# or a stub representing a single process if not. Concentrating all
# MPI calls in this library allows all other code to be compiled
# without the MPI extra's.
# ------------------------------------------------------------------

# set the target name
set(TARGET mpi)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} fundamentals)
include_directories(../../SMILE/fundamentals)

# define a user-configurable option to build with MPI support,
# which requires some MPI implementation to be installed on the system
option(BUILD_WITH_MPI "build with MPI support - requires MPI installation")

# if requested, find MPI installation and configure accordingly
if (BUILD_WITH_MPI)
    # find and configure MPI
    find_package(MPI REQUIRED)
    target_link_libraries(${TARGET} ${MPI_CXX_LIBRARIES})
    include_directories(${MPI_CXX_INCLUDE_PATH})
    # add the C++ preprocessor symbol that triggers compilation with MPI
    add_definitions(-DBUILD_WITH_MPI)
    # cause the relevant variables to be shown in non-advanced mode
    mark_as_advanced(CLEAR MPI_CXX_COMPILER MPI_CXX_LIBRARIES MPI_CXX_INCLUDE_PATH)
    mark_as_advanced(FORCE MPI_LIBRARY MPI_EXTRA_LIBRARY)
endif()

# adjust C++ compiler flags to our needs
set(NO_EXTRA_WARNINGS true)  # to avoid warnings in the MPI headers
include("../../SMILE/build/CompilerFlags.cmake")
