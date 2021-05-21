# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library that supports building Voronoi meshes.
# The source code is taken from the Voro++ library written by
# Chris H. Rycroft (Harvard University/Lawrence Berkeley Laboratory)
# at https://github.com/chr1shr/voro (git commit 122531f).
#
# Changes:
#  - rename v_base_wl.cc to v_base_wl.hh
#  - fix warnings in cell.hh and rad_option.hh
#  - remove wall.hh/cc, pre_container.hh/cc and all non-source files
#---------------------------------------------------------------------

# set the target name
set(TARGET voro)

# list the source files in this directory
file(GLOB SOURCES "*.cc")
file(GLOB HEADERS "*.hh")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# adjust C++ compiler flags to our needs
set(NO_EXTRA_WARNINGS true)  # to avoid warnings in the Voro++ code
include("../../SMILE/build/CompilerFlags.cmake")
