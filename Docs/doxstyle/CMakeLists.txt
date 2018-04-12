# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds doxstyle, the doxygen documentation block streamliner
# ------------------------------------------------------------------

# set the target name
set(TARGET doxstyle)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the executable target
add_executable(${TARGET} ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} fundamentals)
include_directories(../../SMILE/fundamentals)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
