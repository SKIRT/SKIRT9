# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds all targets defined in the SKIRT subproject
# ------------------------------------------------------------------

# add all relevant subdirectories; each subdirectory defines a single target
add_subdirectory(fitsio)
add_subdirectory(voro)
add_subdirectory(tetgen)
add_subdirectory(utils)
add_subdirectory(mpi)
add_subdirectory(tools)
add_subdirectory(geometry)
add_subdirectory(sed)
add_subdirectory(wavelength)
add_subdirectory(field)
add_subdirectory(source)
add_subdirectory(grain)
add_subdirectory(material)
add_subdirectory(medium)
add_subdirectory(instrument)
add_subdirectory(probe)
add_subdirectory(simulation)
