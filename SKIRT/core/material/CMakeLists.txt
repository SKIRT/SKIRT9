# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with the MaterialMix hierarchy (dust, electron and
# gas mixes), MaterialMixFamily, and the shared medium-state types
# (MaterialState, StateVariable, MediumState) and calculation helpers
# used only by these classes
# ------------------------------------------------------------------

# set the target name
set(TARGET material)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} schema fundamentals)
include_directories(../../../SMILE/schema ../../../SMILE/fundamentals)

# add SKIRT library dependencies
target_link_libraries(${TARGET} mpi utils tools wavelength grain)
include_directories(../mpi ../utils ../tools ../wavelength ../grain)

# adjust C++ compiler flags to our needs
include("../../../SMILE/build/CompilerFlags.cmake")
