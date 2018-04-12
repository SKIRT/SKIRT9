# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library that supports building Voronoi meshes.
# The source code is adapted from the Voro++ library written by
# Chris H. Rycroft (LBL / UC Berkeley)
# which can be downloaded at http://math.lbl.gov/voro++/about.html
#---------------------------------------------------------------------

# set the target name
set(TARGET voro)

# list the source files in this directory
file(GLOB SOURCES "*.cc")
file(GLOB HEADERS "*.hh")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
