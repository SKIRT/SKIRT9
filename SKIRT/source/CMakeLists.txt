# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with primary emitting sources (Source), their
# angular emission patterns (AngularDistribution) and polarization
# profiles, and the luminosity normalizations used to scale them.
# Secondary (re-emission) sources are deliberately excluded.
# ------------------------------------------------------------------

# set the target name
set(TARGET source)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} schema fundamentals)
include_directories(../../SMILE/schema ../../SMILE/fundamentals)

# add direct SKIRT library dependencies
target_link_libraries(${TARGET} field geometry mpi sed tools utils wavelength)
include_directories(../field ../geometry ../mpi ../sed ../tools ../utils ../wavelength)

# add indirect SKIRT library dependencies (needed transitively via the direct ones above)
target_link_libraries(${TARGET} fitsio voro)
include_directories(../fitsio ../voro)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
