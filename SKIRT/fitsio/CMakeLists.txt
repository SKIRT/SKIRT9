# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library providing input/output capabilities for FITS files.
# It contains a subset of the source files offered in tarball
# cfitsio-3.49.tar.gz downloaded from NASA's HEASARC software site
# at http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html .
# No changes were made to the source files except for removing the
# "#line" compiler directives (to avoid problems on Windows).
#---------------------------------------------------------------------

# set the target name
set(TARGET fitsio)

# list the source files in this directory
file(GLOB SOURCES "*.c")
file(GLOB HEADERS "*.h")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# avoid non-ANSI include
add_definitions(-DFF_NO_UNISTD_H=1)

# suppress all compiler warnings
target_compile_options(${TARGET} PRIVATE -w)
