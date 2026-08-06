# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with tool classes that depend on SimulationItem:
# SimulationItem itself, configuration, logging, random numbers,
# unit systems, cosmology, parallel execution, input output
# (file paths, stored-table/columns, text, FITS).
# ------------------------------------------------------------------

# set the target name
set(TARGET tools)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} schema fundamentals)
include_directories(../../SMILE/schema ../../SMILE/fundamentals)

# add direct SKIRT library dependencies
target_link_libraries(${TARGET} fitsio mpi utils)
include_directories(../fitsio ../mpi ../utils)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
