# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for the SMILE GUI wizard layer (requires Qt5)
# ------------------------------------------------------------------

# set the target name
set(TARGET wizard)

# enable Qt functionality, supporting Qt5 and Qt6
find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core Widgets)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library and Qt5 dependencies
target_link_libraries(${TARGET} serialize schema fundamentals Qt${QT_VERSION_MAJOR}::Core Qt${QT_VERSION_MAJOR}::Widgets)
include_directories(../serialize ../schema ../fundamentals)

# adjust C++ compiler flags to our needs
include("../build/CompilerFlags.cmake")
