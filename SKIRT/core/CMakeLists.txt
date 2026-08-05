# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Groups the targets that make up the SKIRT core layer; this
# directory itself does not define a target
# ------------------------------------------------------------------

# add all relevant subdirectories; each subdirectory defines a single target
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
