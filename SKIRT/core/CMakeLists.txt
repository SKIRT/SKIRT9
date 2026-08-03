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
add_subdirectory(framework)
