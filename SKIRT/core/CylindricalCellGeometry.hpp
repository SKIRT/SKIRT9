/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDRICALCELLGEOMETRY_HPP
#define CYLINDRICALCELLGEOMETRY_HPP

#include "ImportedGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** A CylindricalCellGeometry instance represents a 3D geometry with a spatial density distribution
    described by a list of cylindrical cells lined up with the cylindrical coordinate axes; refer
    to the CylindricalCellSnapshot class for more information. The cell data is usually extracted
    from a cosmological simulation snapshot, and it must be provided in a column text file
    formatted as described below. The total mass in the geometry is normalized to unity after
    importing the data.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns in the input file depends on the options configured by the user for this
    CylindricalCellGeometry instance:

    \f[ R_\mathrm{min}\,(\mathrm{pc}) \quad \varphi_\mathrm{min}\,(\mathrm{deg}) \quad
    z_\mathrm{min}\,(\mathrm{pc}) \quad R_\mathrm{max}\,(\mathrm{pc}) \quad
    \varphi_\mathrm{max}\,(\mathrm{deg}) \quad z_\mathrm{max}\,(\mathrm{pc}) \quad \{\,
    \rho\,(\text{M}_\odot\,\text{pc}^{-3}) \;\;|\;\; M\,(\text{M}_\odot) \;\;|\;\;
    n\,(\text{cm}^{-3}) \;\;|\;\; N\,(1) \,\} \quad [Z\,(1)] \quad [T\,(\mathrm{K})] \f]

    The first six columns specify the coordinates of the bordering planes and cylinders of the
    cell. The \em autoRevolve property controls a feature to automatically revolve 2D data to 3D.
    See the CylindricalCellSnapshot class for more information.

    Depending on the value of the \em massType option, the seventh column lists the average mass
    density \f$\rho\f$, the integrated mass \f$M\f$, the average number density \f$n\f$, or the
    integrated number density \f$N\f$ for the cell. The precise units for this field are irrelevant
    because the total mass in the geometry will be normalized to unity after importing the data.
    However, the import procedure still insists on knowing the units.

    If the \em importMetallicity option is enabled, the next column specifies a "metallicity"
    fraction, which in this context is simply multiplied with the mass/density column to obtain the
    actual mass/density of the cell. If the \em importTemperature option is enabled, the next
    column specifies a temperature. If this temperature is higher than the maximum configured
    temperature, the mass and density of the cell are set to zero, regardless of the mass or
    density specified in the seventh column. If the \em importTemperature option is disabled, or
    the maximum temperature value is set to zero, such a cutoff is not applied. */
class CylindricalCellGeometry : public ImportedGeometry
{
    /** The enumeration type indicating the type of mass quantity to be imported. */
    ENUM_DEF(MassType, MassDensity, Mass, NumberDensity, Number)
        ENUM_VAL(MassType, MassDensity, "mass density")
        ENUM_VAL(MassType, Mass, "mass (volume-integrated mass density)")
        ENUM_VAL(MassType, NumberDensity, "number density")
        ENUM_VAL(MassType, Number, "number (volume-integrated number density)")
    ENUM_END()

    ITEM_CONCRETE(CylindricalCellGeometry, ImportedGeometry, "a geometry imported from cylindrical cell data")
        ATTRIBUTE_TYPE_DISPLAYED_IF(CylindricalCellMedium, "Level2")

        PROPERTY_BOOL(autoRevolve, "automatically revolve 2D data to a 3D model")
        ATTRIBUTE_DEFAULT_VALUE(autoRevolve, "false")

        PROPERTY_INT(numAutoRevolveBins, "the number of azimuth bins for auto-revolving 2D data")
        ATTRIBUTE_RELEVANT_IF(numAutoRevolveBins, "autoRevolve")
        ATTRIBUTE_MIN_VALUE(numAutoRevolveBins, "2")
        ATTRIBUTE_MAX_VALUE(numAutoRevolveBins, "1024")
        ATTRIBUTE_DEFAULT_VALUE(numAutoRevolveBins, "16")

        PROPERTY_ENUM(massType, MassType, "the type of mass quantity to be imported")
        ATTRIBUTE_DEFAULT_VALUE(massType, "MassDensity")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new CylindricalCellGeometry object, calls its open() function,
        configures it to import a mass or density column, and finally returns a pointer to the
        object. Ownership of the Snapshot object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif
