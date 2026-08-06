# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with 3D density distributions (Geometry and its
# decorators), SPH smoothing kernels, and the imported N-body/hydro
# snapshot readers (Snapshot hierarchy) that some geometries build on
# ------------------------------------------------------------------

# set the target name
set(TARGET geometry)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} schema fundamentals)
include_directories(../../SMILE/schema ../../SMILE/fundamentals)

# add direct SKIRT library dependencies
target_link_libraries(${TARGET} mpi tools utils voro)
include_directories(../mpi ../tools ../utils ../voro)

# add indirect SKIRT library dependencies (needed transitively via the direct ones above)
target_link_libraries(${TARGET} fitsio)
include_directories(../fitsio)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
