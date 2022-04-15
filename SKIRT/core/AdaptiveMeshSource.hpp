/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHSOURCE_HPP
#define ADAPTIVEMESHSOURCE_HPP

#include "MeshSource.hpp"

////////////////////////////////////////////////////////////////////

/** An AdaptiveMeshSource instance represents a primary radiation source with a spatial and
    spectral luminosity distribution described by an Adaptive Mesh Refinement (AMR) grid
    partitioning a cuboidal domain. The data is usually extracted from a cosmological simulation
    snapshot, and it must be provided in a column text file formatted as described below.

    Refer to the description of the AdaptiveMeshSnapshot class for information on the structure of
    an adaptive mesh and on how to represent it in text column file format. Refer to the
    description of the TextInFile class for information on overall formatting and on how to include
    header lines specifying the units for each column in the input file. In case the input file has
    no unit specifications, the default units mentioned below are used instead. The number of
    columns expected in the input file depends on the options configured by the user for this
    AdaptiveMeshSource instance, including the selected %SEDFamily.

    \f[ [ v_x\,(\mathrm{km/s}) \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) \quad [
    \sigma_v\,(\mathrm{km/s}) ] ] \quad [M_\mathrm{curr}\,(\mathrm{M}_\odot)] \quad \dots
    \text{SED family parameters}\dots \f]

    If the \em importVelocity option is enabled, the first three columns specify the \f$v_x\f$,
    \f$v_y\f$, \f$v_z\f$ bulk velocity components of the source population represented by the cell.
    If additionally the \em importVelocityDispersion option is enabled, the next column specifies
    the velocity dispersion \f$\sigma_v\f$, adjusting the velocity for each photon packet launch
    with a random offset sampled from a spherically symmetric Gaussian distribution.

    The remaining columns specify the parameters required by the configured %SED family to select
    and scale the appropriate %SED. For example for the Bruzual-Charlot %SED family, the remaining
    columns provide the initial mass, the metallicity, and the age of the stellar population
    represented by the cell corresponding to the site. Refer to the documentation of the configured
    type of SEDFamily for information about the expected parameters and their default units. */
class AdaptiveMeshSource : public MeshSource
{
    ITEM_CONCRETE(AdaptiveMeshSource, MeshSource,
                  "a primary source imported from data represented on an adaptive mesh (AMR grid)")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new AdaptiveMeshSnapshot object, calls its open() function,
        passes it the domain extent configured by the user (using properties offered by the
        MeshSource base class), and returns a pointer to the object. Ownership of the Snapshot
        object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif
