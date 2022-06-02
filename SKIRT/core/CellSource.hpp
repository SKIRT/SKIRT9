/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CELLSOURCE_HPP
#define CELLSOURCE_HPP

#include "ImportedSource.hpp"

////////////////////////////////////////////////////////////////////

/** A CellSource instance represents a primary radiation source with a spatial and spectral
    luminosity distribution described by a list of cuboidal cells lined up with the coordinate
    axes; refer to the CellSnapshot class for more information. The cell data is usually extracted
    from a cosmological simulation snapshot, and it must be provided in a column text file
    formatted as described below.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns expected in the input file depends on the options configured by the user for
    this CellSource instance, including the selected %SEDFamily.

    \f[ x_\mathrm{min}\,(\mathrm{pc}) \quad y_\mathrm{min}\,(\mathrm{pc}) \quad
    z_\mathrm{min}\,(\mathrm{pc}) \quad x_\mathrm{max}\,(\mathrm{pc}) \quad
    y_\mathrm{max}\,(\mathrm{pc}) \quad z_\mathrm{max}\,(\mathrm{pc}) \quad [ v_x\,(\mathrm{km/s})
    \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) \quad [ \sigma_v\,(\mathrm{km/s}) ] ]
    \quad [M_\mathrm{curr}\,(\mathrm{M}_\odot)] \quad \dots \text{SED family parameters}\dots \f]

    The first six columns specify the coordinates of the lower-left and upper-right corners of the
    cell. If the \em importVelocity option is enabled, the next three columns specify the
    \f$v_x\f$, \f$v_y\f$, \f$v_z\f$ bulk velocity components of the source population represented
    by the cell. If additionally the \em importVelocityDispersion option is enabled, the next
    column specifies the velocity dispersion \f$\sigma_v\f$, adjusting the velocity for each photon
    packet launch with a random offset sampled from a spherically symmetric Gaussian distribution.

    The remaining columns specify the parameters required by the configured %SED family to select
    and scale the appropriate %SED. For example for the Bruzual-Charlot %SED family, the remaining
    columns provide the initial mass, the metallicity, and the age of the stellar population
    represented by the cell. Refer to the documentation of the configured type of SEDFamily for
    information about the expected parameters and their default units. */
class CellSource : public ImportedSource
{
    ITEM_CONCRETE(CellSource, ImportedSource, "a primary source imported from cuboidal cell data")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new CellSnapshot object, calls its open() function, and returns
        a pointer to the object. Ownership of the Snapshot object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif
