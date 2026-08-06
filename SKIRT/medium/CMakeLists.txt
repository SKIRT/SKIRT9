# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with MediumSystem, the Medium component hierarchy,
# the spatial cell library classes, the SpatialGrid hierarchy (with
# its tree-node and tree-policy helpers), the 1D Mesh hierarchy, the
# secondary source classes and SecondarySourceSystem, and the
# material normalization and dynamic state recipe classes used only
# by these
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

# add SKIRT library dependencies
target_link_libraries(${TARGET} fitsio voro mpi utils tetgen tools geometry wavelength field grain material)
include_directories(../fitsio ../voro ../mpi ../utils ../tools ../geometry ../wavelength ../field ../grain ../material)
include_directories(SYSTEM ../tetgen)  # suppress warnings in tetgen header

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
