/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERICALCELLSOURCE_HPP
#define SPHERICALCELLSOURCE_HPP

#include "ImportedSource.hpp"

////////////////////////////////////////////////////////////////////

/** A SphericalCellSource instance represents a primary radiation source with a spatial and
    spectral luminosity distribution described by a list of spherical cells lined up with the
    spherical coordinate axes; refer to the SphericalCellSnapshot class for more information. The
    cell data is usually extracted from a cosmological simulation snapshot, and it must be provided
    in a column text file formatted as described below.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns expected in the input file depends on the options configured by the user for
    this SphericalCellSource instance, including the selected %SEDFamily:

    \f[ r_\mathrm{min}\,(\mathrm{pc}) \quad \theta_\mathrm{min}\,(\mathrm{deg}) \quad
    \varphi_\mathrm{min}\,(\mathrm{deg}) \quad r_\mathrm{max}\,(\mathrm{pc}) \quad
    \theta_\mathrm{max}\,(\mathrm{deg}) \quad \varphi_\mathrm{max}\,(\mathrm{deg}) \quad [
    v_x\,(\mathrm{km/s}) \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) \quad [
    \sigma_v\,(\mathrm{km/s}) ] \quad [M_\mathrm{curr}\,(\mathrm{M}_\odot)] \quad b\,(1) \quad
    \dots \text{SED family parameters} \dots \f]

    The first six columns specify the coordinates of the bordering surfaces of the cell. The \em
    autoRevolve property controls a feature to automatically revolve 2D data to 3D. See the
    SphericalCellSnapshot class for more information.

    If the \em importVelocity option is enabled, the next three columns specify the \f$v_r\f$,
    \f$v_\theta\f$, \f$v_\varphi\f$ components of the bulk velocity, in spherical coordinates, for
    the source population represented by the cell. If additionally the \em importVelocityDispersion
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
class SphericalCellSource : public ImportedSource
{
    /** The enumeration type for selecting the auto-revolve option. */
    ENUM_DEF(AutoRevolveType, None, Inclination, Azimuth, InclinationAndAzimuth)
        ENUM_VAL(AutoRevolveType, None, "no auto-revolve")
        ENUM_VAL(AutoRevolveType, Inclination, "auto-revolve inclination")
        ENUM_VAL(AutoRevolveType, Azimuth, "auto-revolve azimuth")
        ENUM_VAL(AutoRevolveType, InclinationAndAzimuth, "auto-revolve both inclination and azimuth")
    ENUM_END()

    ITEM_CONCRETE(SphericalCellSource, ImportedSource, "a primary source imported from spherical cell data")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SphericalCellSource, "Level2")

        PROPERTY_ENUM(autoRevolve, AutoRevolveType, "automatically revolve 1D or 2D data to a 3D model")
        ATTRIBUTE_DEFAULT_VALUE(autoRevolve, "None")

        PROPERTY_INT(numInclinationRevolveBins, "the number of inclination bins for auto-revolving 1D or 2D data")
        ATTRIBUTE_RELEVANT_IF(numInclinationRevolveBins, "autoRevolveInclination|autoRevolveInclinationAndAzimuth")
        ATTRIBUTE_MIN_VALUE(numInclinationRevolveBins, "2")
        ATTRIBUTE_MAX_VALUE(numInclinationRevolveBins, "1024")
        ATTRIBUTE_DEFAULT_VALUE(numInclinationRevolveBins, "16")

        PROPERTY_INT(numAzimuthRevolveBins, "the number of azimuth bins for auto-revolving 1D or 2D data")
        ATTRIBUTE_RELEVANT_IF(numAzimuthRevolveBins, "autoRevolveAzimuth|autoRevolveInclinationAndAzimuth")
        ATTRIBUTE_MIN_VALUE(numAzimuthRevolveBins, "2")
        ATTRIBUTE_MAX_VALUE(numAzimuthRevolveBins, "1024")
        ATTRIBUTE_DEFAULT_VALUE(numAzimuthRevolveBins, "16")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new SphericalCellSource object, calls its open() function, and
        returns a pointer to the object. Ownership of the Snapshot object is transferred to the
        caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif
