# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Calculates the Faddeeva function.
# The package is taken from http://ab-initio.mit.edu/faddeeva//
# This is specifically used to evaluate the Voigt profile by setting H(x, a) = Real(w(x+ia)).
#---------------------------------------------------------------------

# set the target name
set(TARGET faddeeva)

# list the source files in this directory
file(GLOB SOURCES "*.cc")
file(GLOB HEADERS "*.hh")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# adjust C++ compiler flags to our needs
set(NO_EXTRA_WARNINGS true)  # to avoid warnings in the faddeeva code (warnings still present)
include("../../SMILE/build/CompilerFlags.cmake")


