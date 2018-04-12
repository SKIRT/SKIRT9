# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for the SMILE GUI wizard layer (requires Qt5)
# ------------------------------------------------------------------

# set the target name
set(TARGET wizard)

# enable Qt functionality
find_package(Qt5Widgets)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library and Qt5 dependencies
target_link_libraries(${TARGET} serialize schema fundamentals Qt5::Widgets)
include_directories(../serialize ../schema ../fundamentals)

# adjust C++ compiler flags to our needs
include("../build/CompilerFlags.cmake")
