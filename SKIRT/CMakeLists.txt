# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds all targets defined in the SKIRT subproject
# ------------------------------------------------------------------

# add all relevant subdirectories; each subdirectory defines a single target,
# except for core, which just groups the utils/mpi/framework targets
add_subdirectory(fitsio)
add_subdirectory(voro)
add_subdirectory(tetgen)
add_subdirectory(core)
add_subdirectory(main)
