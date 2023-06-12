# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds the MakeUp desktop GUI wizard for editing SMILE data sets
# ------------------------------------------------------------------

# set the target name
set(TARGET MakeUp)

# enable Qt functionality, supporting Qt5 and Qt6
find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core Widgets)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# add the application icon to the source files
if (CMAKE_HOST_APPLE)
    set(SOURCES ${SOURCES} mirror.icns)
    set_source_files_properties(mirror.icns PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
endif()
if (CMAKE_HOST_WIN32)
    set(SOURCES ${SOURCES} mirror.ico MakeUp.rc)
endif()

# if SKIRT is available, have it generate a SMILE schema file
set(SMILE_SCHEMA_DIR ${CMAKE_CURRENT_BINARY_DIR}/../schemas)
file(MAKE_DIRECTORY ${SMILE_SCHEMA_DIR})
if (BUILD_SKIRT)
    add_custom_command(OUTPUT ${SMILE_SCHEMA_DIR}/skirt.smile
                       COMMAND skirt -x
                       DEPENDS skirtcore
                       WORKING_DIRECTORY ${SMILE_SCHEMA_DIR})
    set(SOURCES ${SOURCES} ${SMILE_SCHEMA_DIR}/skirt.smile)
endif()

# create the executable target
add_executable(${TARGET} ${SOURCES} ${HEADERS})
target_link_libraries (${TARGET} Qt${QT_VERSION_MAJOR}::Core Qt${QT_VERSION_MAJOR}::Widgets)

# add SMILE library dependencies
target_link_libraries(${TARGET} wizard fundamentals build)
include_directories(../../SMILE/wizard ../../SMILE/fundamentals ../../SMILE/build)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")

# get the project-wide version string
include ("../../SMILE/build/Version.cmake")

# configure the executable as a GUI application
if (CMAKE_HOST_APPLE)
    set_target_properties(${TARGET} PROPERTIES
                                        MACOSX_BUNDLE true
                                        MACOSX_BUNDLE_INFO_STRING "MakeUp -- wizard for editing SMILE data sets"
                                        MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}"
                                        MACOSX_BUNDLE_COPYRIGHT "© Astronomical Observatory, Ghent University"
                                        MACOSX_BUNDLE_GUI_IDENTIFIER "be.ugent.observatory.make-up"
                                        MACOSX_BUNDLE_ICON_FILE "mirror.icns")
endif()
if (CMAKE_HOST_WIN32)
    set_target_properties(${TARGET} PROPERTIES WIN32_EXECUTABLE true)
endif()
