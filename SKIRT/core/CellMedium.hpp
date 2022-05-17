/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CELLMEDIUM_HPP
#define CELLMEDIUM_HPP

#include "ImportedMedium.hpp"

////////////////////////////////////////////////////////////////////

/** A CellMedium instance represents a transfer medium with a spatial density distribution (and,
    optionally, other properties) described by a list of cuboidal cells lined up with the
    coordinate axes; refer to the CellSnapshot class for more information. The cell data is usually
    extracted from a cosmological simulation snapshot, and it must be provided in a column text
    file formatted as described below.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns expected in the input file depends on the options configured by the user for
    this CellMedium instance:

    \f[ x_\mathrm{min}\,(\mathrm{pc}) \quad y_\mathrm{min}\,(\mathrm{pc}) \quad
    z_\mathrm{min}\,(\mathrm{pc}) \quad x_\mathrm{max}\,(\mathrm{pc}) \quad
    y_\mathrm{max}\,(\mathrm{pc}) \quad z_\mathrm{max}\,(\mathrm{pc}) \quad \{\,
    \rho\,(\text{M}_\odot\,\text{pc}^{-3}) \;\;|\;\; M\,(\text{M}_\odot) \;\;|\;\;
    n\,(\text{cm}^{-3}) \;\;|\;\; N\,(1) \,\} \quad [Z\,(1)] \quad [T\,(\mathrm{K})] \f] \f[ [
    v_x\,(\mathrm{km/s}) \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) ] \quad [
    B_x\,(\mu\mathrm{G}) \quad B_y\,(\mu\mathrm{G}) \quad B_z\,(\mu\mathrm{G}) ] \quad [
    \dots\text{mix family params}\dots ] \f]

    The first six columns specify the coordinates of the lower-left and upper-right corners of the
    cell. Depending on the value of the \em massType option, the seventh column lists the average
    mass density \f$\rho\f$, the integrated mass \f$M\f$, the average number density \f$n\f$, or
    the integrated number density \f$N\f$ for the cell.

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
    \f$v_x\f$, \f$v_y\f$, \f$v_z\f$ velocity components of the particle (considered as the bulk
    velocity for the mass represented by the particle).

    If the \em importMagneticField option is enabled, the subsequent three columns specify the
    \f$B_x\f$, \f$B_y\f$, \f$B_z\f$ magnetic field vector components for the cell.

    Finally, if the \em importVariableMixParams option is enabled, the remaining columns specify
    the parameters used by the configured material mix family to select a particular material mix
    for the cell. */
class CellMedium : public ImportedMedium
{
    /** The enumeration type indicating the type of mass quantity to be imported. */
    ENUM_DEF(MassType, MassDensity, Mass, NumberDensity, Number)
        ENUM_VAL(MassType, MassDensity, "mass density")
        ENUM_VAL(MassType, Mass, "mass (volume-integrated)")
        ENUM_VAL(MassType, NumberDensity, "number density")
        ENUM_VAL(MassType, Number, "number (volume-integrated)")
    ENUM_END()

    ITEM_CONCRETE(CellMedium, ImportedMedium, "a transfer medium imported from cuboidal cell data")

        PROPERTY_ENUM(massType, MassType, "the type of mass quantity to be imported")
        ATTRIBUTE_DEFAULT_VALUE(massType, "MassDensity")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new CellSnapshot object, calls its open() function, configures
        it to import a mass or density column, it, and finally returns a pointer to the object.
        Ownership of the Snapshot object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif
