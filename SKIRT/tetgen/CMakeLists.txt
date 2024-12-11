# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library that supports building Tetra meshes.
# The source code is taken from the TetGen library written by
# Hang Si (Weierstrass Institute for Applied Analysis and Stochastics (WIAS))
# at https://github.com/libigl/tetgen (git commit 4f3bfba).
#
# Changes:
#  - changed TetSizeFunc from function pointer to std::function so lambda with captures can be used
#  - adjusted the calling of the tetunsuitable to remove the last two unused parameters
#---------------------------------------------------------------------

# set the target name
set(TARGET tetgen)

# list the source files in this directory
file(GLOB SOURCES "*.cxx")
file(GLOB HEADERS "*.h")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})
target_compile_definitions(tetgen PRIVATE -DTETLIBRARY) # to compile TetGen as a library

# adjust C++ compiler flags to our needs
set(NO_EXTRA_WARNINGS true)  # to avoid warnings in the tetgen code (warnings still present)
include("../../SMILE/build/CompilerFlags.cmake")


