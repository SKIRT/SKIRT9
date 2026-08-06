# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds all targets defined in the SKIRT subproject
# ------------------------------------------------------------------

# add all relevant subdirectories; each subdirectory defines a single target
add_subdirectory(field)
add_subdirectory(fitsio)
add_subdirectory(geometry)
add_subdirectory(grain)
add_subdirectory(instrument)
add_subdirectory(material)
add_subdirectory(medium)
add_subdirectory(mpi)
add_subdirectory(probe)
add_subdirectory(sed)
add_subdirectory(simulation)
add_subdirectory(source)
add_subdirectory(tetgen)
add_subdirectory(tools)
add_subdirectory(utils)
add_subdirectory(voro)
add_subdirectory(wavelength)
