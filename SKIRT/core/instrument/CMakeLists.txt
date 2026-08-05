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
include_directories(../../../SMILE/schema ../../../SMILE/fundamentals)

# add SKIRT library dependencies
target_link_libraries(${TARGET} fitsio voro mpi utils tetgen tools geometry wavelength grain material medium)
include_directories(../../fitsio ../../voro ../mpi ../utils ../tools ../geometry ../wavelength ../grain ../material ../medium)
include_directories(SYSTEM ../../tetgen)  # suppress warnings in tetgen header

# adjust C++ compiler flags to our needs
include("../../../SMILE/build/CompilerFlags.cmake")
