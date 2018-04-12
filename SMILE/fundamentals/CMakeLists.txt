# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for the SMILE fundamentals layer
# ------------------------------------------------------------------

# set the target name
set(TARGET fundamentals)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# on Windows, add some system libraries:
#  - advapi32: for getting user name
#  - dbghelp: for constructing a stack trace
#  - ws2_32: for getting host name
#  - Pathcch: for getting canonical path
if (CMAKE_HOST_WIN32)
    target_link_libraries(${TARGET} advapi32 dbghelp ws2_32 Pathcch)
endif()

# on Mac, add the core foundation framework (for getting the executable path)
if (CMAKE_HOST_APPLE)
    find_library(CORE_FOUNDATION CoreFoundation)
    target_link_libraries(${TARGET} ${CORE_FOUNDATION})
endif()

# adjust C++ compiler flags to our needs
include("../build/CompilerFlags.cmake")
