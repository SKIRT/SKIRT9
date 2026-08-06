# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with spectral energy distributions (SED) and
# parameterized SED libraries (SEDFamily) from which an SED can be
# obtained by interpolation
# ------------------------------------------------------------------

# set the target name
set(TARGET sed)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} schema fundamentals)
include_directories(../../SMILE/schema ../../SMILE/fundamentals)

# add direct SKIRT library dependencies
target_link_libraries(${TARGET} tools utils)
include_directories(../tools ../utils)

# add indirect SKIRT library dependencies (needed transitively via the direct ones above)
target_link_libraries(${TARGET} fitsio mpi)
include_directories(../fitsio ../mpi)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
