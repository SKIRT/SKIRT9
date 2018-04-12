# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Provide version and build information
# ------------------------------------------------------------------

# set the target name
set(TARGET build)

# create the library target
add_library(${TARGET} STATIC BuildInfo.cpp BuildInfo.hpp)

# add include path for Basics.hpp
include_directories(../fundamentals)

# add include path for headers created by CMake in the build directory
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# adjust C++ compiler flags to our needs
include("CompilerFlags.cmake")

# ------------------------------------------------------------------

# get the project-wide version string
include("Version.cmake")

# create the version.h header containing the version string
configure_file("version_in.h" "${CMAKE_CURRENT_BINARY_DIR}/version.h")

# ------------------------------------------------------------------

# write a cmake script that writes the timestamp.h header containing a build timestamp and the git commit hash
set(TS_PATH "${PROJECT_BINARY_DIR}/timestamp.cmake")
file(WRITE ${TS_PATH}  "string(TIMESTAMP BUILD_DATE \"%d/%m/%Y\")\n")
file(APPEND ${TS_PATH} "string(TIMESTAMP BUILD_TIME \"%H:%M:%S\")\n")
find_program(GIT_EXECUTABLE git)
if (GIT_EXECUTABLE)
  file(APPEND ${TS_PATH} "execute_process(COMMAND \"${GIT_EXECUTABLE}\" describe --dirty --always\n")
  file(APPEND ${TS_PATH} "    WORKING_DIRECTORY \"${CMAKE_SOURCE_DIR}\" OUTPUT_VARIABLE GIT_COMMIT_HASH\n")
  file(APPEND ${TS_PATH} "    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)\n")
else()
  file(APPEND ${TS_PATH} "set(GIT_COMMIT_HASH 0)\n")
endif()
file(APPEND ${TS_PATH} "file(WRITE timestamp.h \"#ifndef TIMESTAMP_H\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#define TIMESTAMP_H\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#define BUILD_DATE \\\"\${BUILD_DATE}\\\"\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#define BUILD_TIME \\\"\${BUILD_TIME}\\\"\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#define COMMIT_HASH \\\"\${GIT_COMMIT_HASH}\\\"\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#endif\\n\")\n")

# provide a custom target that will execute the above cmake script at build time
add_custom_target(timestamp
    COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BINARY_DIR}/timestamp.cmake
    ADD_DEPENDENCIES ${PROJECT_BINARY_DIR}/timestamp.cmake)

# add a dependency to ensure that the above custom target is executed for each build
add_dependencies(${TARGET} timestamp)

# ------------------------------------------------------------------
