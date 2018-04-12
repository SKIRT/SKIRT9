# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds all targets defined in the SMILE subproject
# ------------------------------------------------------------------

# define user-configurable options to build the smiletool and /or shapes executables;
# the wizard library is built when the calling project's BUILD_MAKE_UP option is enabled
option(BUILD_SMILE_TOOL "build smiletool, the SMILE console Q&A tool")
option(BUILD_SMILE_SHAPES "build shapes, a basic example for using SMILE")

# add all relevant subdirectories; each subdirectory defines a single target
add_subdirectory(build)
add_subdirectory(fundamentals)
add_subdirectory(schema)
add_subdirectory(serialize)
if (BUILD_MAKE_UP)
    add_subdirectory(wizard)
endif()
if (BUILD_SMILE_TOOL)
    add_subdirectory(smiletool)
endif()
if (BUILD_SMILE_SHAPES)
    add_subdirectory(shapes)
endif()
