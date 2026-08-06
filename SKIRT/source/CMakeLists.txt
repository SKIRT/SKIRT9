# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with primary emitting sources (Source), their
# angular emission patterns (AngularDistribution) and polarization
# profiles, and the luminosity normalizations used to scale them.
# Secondary (dust/line re-emission) sources are deliberately excluded
# for now; they belong elsewhere later.
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

# add SKIRT library dependencies
target_link_libraries(${TARGET} utils mpi tools geometry sed wavelength field)
include_directories(../utils ../mpi ../tools ../geometry ../sed ../wavelength ../field)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
