/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDRICALCELLSOURCE_HPP
#define CYLINDRICALCELLSOURCE_HPP

#include "ImportedSource.hpp"

////////////////////////////////////////////////////////////////////

/** A CylindricalCellSource instance represents a primary radiation source with a spatial and
    spectral luminosity distribution described by a list of cylindrical cells lined up with the
    cylindrical coordinate axes; refer to the CylindricalCellSnapshot class for more information.
    The cell data is usually extracted from a cosmological simulation snapshot, and it must be
    provided in a column text file formatted as described below.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns expected in the input file depends on the options configured by the user for
    this CylindricalCellSource instance, including the selected %SEDFamily:

    \f[ R_\mathrm{min}\,(\mathrm{pc}) \quad \varphi_\mathrm{min}\,(\mathrm{deg}) \quad
    z_\mathrm{min}\,(\mathrm{pc}) \quad R_\mathrm{max}\,(\mathrm{pc}) \quad
    \varphi_\mathrm{max}\,(\mathrm{deg}) \quad z_\mathrm{max}\,(\mathrm{pc}) \quad [
    v_x\,(\mathrm{km/s}) \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) \quad [
    \sigma_v\,(\mathrm{km/s}) ] \quad [M_\mathrm{curr}\,(\mathrm{M}_\odot)] \quad b\,(1) \quad
    \dots \text{SED family parameters} \dots \f]

    The first six columns specify the coordinates of the bordering planes and cylinders of the
    cell. The \em autoRevolve property controls a feature to automatically revolve 2D data to 3D.
    See the CylindricalCellSnapshot class for more information.

    If the \em importVelocity option is enabled, the next three columns specify the \f$v_R\f$,
    \f$v_\varphi\f$, \f$v_z\f$ components of the bulk velocity, in cylindrical coordinates, for the
    source population represented by the cell. If additionally the \em importVelocityDispersion
    option is enabled, the next column specifies the velocity dispersion \f$\sigma_v\f$, adjusting
    the velocity for each photon packet launch with a random offset sampled from a spherically
    symmetric Gaussian distribution. If the \em importCurrentMass option is enabled, the next
    column provides the current mass of the cell, \f$M_\mathrm{curr}\f$. This mass is currently
    only used for probing the input model. If the \em importBias option is enabled, the next column
    specifies the bias parameter, \f$b\f$, which is used to bias the photon sampling for each cell
    (see the documentation of the ImportedSource class).

    The remaining columns specify the parameters required by the configured %SED family to select
    and scale the appropriate %SED. For example for the Bruzual-Charlot %SED family, the remaining
    columns provide the initial mass, the metallicity, and the age of the stellar population
    represented by the cell. Refer to the documentation of the configured type of SEDFamily for
    information about the expected parameters and their default units. */
class CylindricalCellSource : public ImportedSource
{
    ITEM_CONCRETE(CylindricalCellSource, ImportedSource, "a primary source imported from cylindrical cell data")
        ATTRIBUTE_TYPE_DISPLAYED_IF(CylindricalCellSource, "Level2")

        PROPERTY_BOOL(autoRevolve, "automatically revolve 2D data to a 3D model")
        ATTRIBUTE_DEFAULT_VALUE(autoRevolve, "false")

        PROPERTY_INT(numAutoRevolveBins, "the number of azimuth bins for auto-revolving 2D data")
        ATTRIBUTE_RELEVANT_IF(numAutoRevolveBins, "autoRevolve")
        ATTRIBUTE_MIN_VALUE(numAutoRevolveBins, "2")
        ATTRIBUTE_MAX_VALUE(numAutoRevolveBins, "1024")
        ATTRIBUTE_DEFAULT_VALUE(numAutoRevolveBins, "16")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new CylindricalCellSource object, calls its open() function, and
        returns a pointer to the object. Ownership of the Snapshot object is transferred to the
        caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif
