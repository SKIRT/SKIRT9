# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with InstrumentSystem, the Instrument hierarchy,
# the all-sky projection helpers used by AllSkyInstrument and
# HEALPixSkyInstrument, the TimeGrid hierarchy used by TimeInstrument
# and FluxRecorder, and FluxRecorder itself
# ------------------------------------------------------------------

# set the target name
set(TARGET instrument)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} schema fundamentals)
include_directories(../../SMILE/schema ../../SMILE/fundamentals)

# add direct SKIRT library dependencies
target_link_libraries(${TARGET} medium mpi tools utils wavelength)
include_directories(../medium ../mpi ../tools ../utils ../wavelength)

# add indirect SKIRT library dependencies (needed transitively via the direct ones above)
target_link_libraries(${TARGET} field fitsio geometry grain material tetgen voro)
include_directories(../field ../fitsio ../geometry ../grain ../material ../voro)
include_directories(SYSTEM ../tetgen)  # suppress warnings in tetgen header

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
