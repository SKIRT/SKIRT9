# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds the SKIRT executable: the command-line entry point
# (SkirtMain, SkirtCommandLineHandler) together with the top-level
# simulation orchestration classes (Simulation, MonteCarloSimulation,
# ConfigurationSetup, SimulationItemRegistry) that only this
# executable ever needs
# ------------------------------------------------------------------

# set the target name
set(TARGET simulation)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the executable target; the produced binary is named "skirt" and, for
# backward compatibility with procedures that hardcode this path, is placed
# under SKIRT/main in the build tree rather than under this source directory
add_executable(${TARGET} ${SOURCES} ${HEADERS})
set_target_properties(${TARGET} PROPERTIES
    OUTPUT_NAME skirt
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/SKIRT/main
)

# enable multi-threading
find_package(Threads REQUIRED)
target_link_libraries(${TARGET} Threads::Threads)

# add SMILE library dependencies
target_link_libraries(${TARGET} serialize schema fundamentals build)
include_directories(../../SMILE/serialize ../../SMILE/schema ../../SMILE/fundamentals ../../SMILE/build)

# add direct SKIRT library dependencies
target_link_libraries(${TARGET} field geometry grain instrument material medium mpi probe sed source tools utils wavelength)
include_directories(../field ../geometry ../grain ../instrument ../material ../medium ../mpi ../probe ../sed ../source ../tools ../utils ../wavelength)

# add indirect SKIRT library dependencies (needed transitively via the direct ones above)
target_link_libraries(${TARGET} fitsio tetgen voro)
include_directories(../fitsio ../voro)
include_directories(SYSTEM ../tetgen)  # suppress warnings in tetgen header

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
