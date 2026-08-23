# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with MediumSystem and the Medium component hierarchy,
# the material normalization classes, the spatial cell library classes,
# the dynamic state recipe classes, the SpatialGrid hierarchy with
# its tree-node and tree-policy helpers, the 1D Mesh hierarchy, the
# SecondarySourceSystem, and the secondary source classes.
# ------------------------------------------------------------------

# set the target name
set(TARGET medium)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} schema fundamentals)
include_directories(../../SMILE/schema ../../SMILE/fundamentals)

# add direct SKIRT library dependencies
target_link_libraries(${TARGET} field geometry material mpi tetgen tools utils wavelength)
include_directories(../field ../geometry ../material ../mpi ../tools ../utils ../wavelength)
include_directories(SYSTEM ../tetgen)  # suppress warnings in tetgen header

# add indirect SKIRT library dependencies (needed transitively via the direct ones above)
target_link_libraries(${TARGET} fitsio grain voro)
include_directories(../fitsio ../grain ../voro)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
