/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERICALCELLMEDIUM_HPP
#define SPHERICALCELLMEDIUM_HPP

#include "ImportedMedium.hpp"

////////////////////////////////////////////////////////////////////

/** A SphericalCellMedium instance represents a transfer medium with a spatial density distribution
    (and, optionally, other properties) described by a list of spherical cells lined up with the
    spherical coordinate axes; refer to the SphericalCellSnapshot class for more information. The
    cell data is usually extracted from a cosmological simulation snapshot, and it must be provided
    in a column text file formatted as described below.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns expected in the input file depends on the options configured by the user for
    this SphericalCellMedium instance:

    \f[ r_\mathrm{min}\,(\mathrm{pc}) \quad \theta_\mathrm{min}\,(\mathrm{deg}) \quad
    \varphi_\mathrm{min}\,(\mathrm{deg}) \quad r_\mathrm{max}\,(\mathrm{pc}) \quad
    \theta_\mathrm{max}\,(\mathrm{deg}) \quad \varphi_\mathrm{max}\,(\mathrm{deg}) \quad \{\,
    \rho\,(\text{M}_\odot\,\text{pc}^{-3}) \;\;|\;\; M\,(\text{M}_\odot) \;\;|\;\;
    n\,(\text{cm}^{-3}) \;\;|\;\; N\,(1) \,\} \quad [Z\,(1)] \quad [T\,(\mathrm{K})] \f] \f[ [
    v_x\,(\mathrm{km/s}) \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) ] \quad [
    B_x\,(\mu\mathrm{G}) \quad B_y\,(\mu\mathrm{G}) \quad B_z\,(\mu\mathrm{G}) ] \quad [
    \dots\text{mix family params}\dots ] \f]

    The first six columns specify the coordinates of the bordering surfaces of the cell. The \em
    autoRevolve property controls a feature to automatically revolve 1D or 2D data to 3D. See the
    SphericalCellSnapshot class for more information.

    Depending on the value of the \em massType option, the seventh column lists the average mass
    density \f$\rho\f$, the integrated mass \f$M\f$, the average number density \f$n\f$, or the
    integrated number density \f$N\f$ for the cell.

    If the \em importMetallicity option is enabled, the next column specifies a "metallicity"
    fraction, which is multiplied with the mass/density column to obtain the actual mass/density of
    the cell. If the \em importTemperature option is enabled, the next column specifies a
    temperature. If this temperature is higher than the maximum configured temperature, the mass
    and density of the cell are set to zero, regardless of the mass or density specified in the
    seventh column. If the \em importTemperature option is disabled, or the maximum temperature
    value is set to zero, such a cutoff is not applied.

    If both the \em importMetallicity and \em importTemperature options are enabled, this leads to
    the following expression for the mass of an imported particle:

    \f[ M_\mathrm{imported} = \begin{cases} f_\mathrm{mass}\,Z\,M & \mathrm{if}\; T<T_\mathrm{max}
    \;\mathrm{or}\; T_\mathrm{max}=0 \\ 0 & \mathrm{otherwise} \end{cases} \f]

    If the \em importVelocity option is enabled, the subsequent three columns specify the
    \f$v_R\f$, \f$v_\varphi\f$, \f$v_z\f$ components of the bulk velocity, in spherical
    coordinates, for the material represented by the cell.

    If the \em importMagneticField option is enabled, the subsequent three columns specify the
    \f$B_R\f$, \f$B_\varphi\f$, \f$B_z\f$ magnetic field vector components for the cell, in
    spherical coordinates.

    Finally, if the \em importVariableMixParams option is enabled, the remaining columns specify
    the parameters used by the configured material mix family to select a particular material mix
    for the cell. */
class SphericalCellMedium : public ImportedMedium
{
    /** The enumeration type indicating the type of mass quantity to be imported. */
    ENUM_DEF(MassType, MassDensity, Mass, NumberDensity, Number)
        ENUM_VAL(MassType, MassDensity, "mass density")
        ENUM_VAL(MassType, Mass, "mass (volume-integrated)")
        ENUM_VAL(MassType, NumberDensity, "number density")
        ENUM_VAL(MassType, Number, "number (volume-integrated)")
    ENUM_END()

    /** The enumeration type for selecting the auto-revolve option. */
    ENUM_DEF(AutoRevolveType, None, Inclination, Azimuth, InclinationAndAzimuth)
        ENUM_VAL(AutoRevolveType, None, "no auto-revolve")
        ENUM_VAL(AutoRevolveType, Inclination, "auto-revolve inclination")
        ENUM_VAL(AutoRevolveType, Azimuth, "auto-revolve azimuth")
        ENUM_VAL(AutoRevolveType, InclinationAndAzimuth, "auto-revolve both inclination and azimuth")
    ENUM_END()

    ITEM_CONCRETE(SphericalCellMedium, ImportedMedium, "a transfer medium imported from spherical cell data")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SphericalCellMedium, "Level2")

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

        PROPERTY_ENUM(massType, MassType, "the type of mass quantity to be imported")
        ATTRIBUTE_DEFAULT_VALUE(massType, "MassDensity")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new SphericalCellSnapshot object, calls its open() function,
        configures it to import a mass or density column, and finally returns a pointer to the
        object. Ownership of the Snapshot object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif
