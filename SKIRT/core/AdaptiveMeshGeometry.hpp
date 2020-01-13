/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHGEOMETRY_HPP
#define ADAPTIVEMESHGEOMETRY_HPP

#include "AdaptiveMeshInterface.hpp"
#include "MeshGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** An AdaptiveMeshGeometry instance represents a 3D geometry with a spatial density distribution
    described by an Adaptive Mesh Refinement (AMR) grid partitioning a cuboidal domain. The data is
    usually extracted from a cosmological simulation snapshot, and it must be provided in a column
    text file formatted as described below. The total mass in the geometry is normalized to unity
    after importing the data.

    Refer to the description of the AdaptiveMeshSnapshot class for information on the structure of
    an adaptive mesh and on how to represent it in text column file format. Refer to the
    description of the TextInFile class for information on overall formatting and on how to include
    header lines specifying the units for each column in the input file. In case the input file has
    no unit specifications, the default units mentioned below are used instead. The input file
    should contain 1, 2, or 3 columns, depending on the options configured by the user for this
    AdaptiveMeshGeometry instance:

    \f[ \{\, \rho\,(\text{M}_\odot\,\text{pc}^{-3}) \;\;|\;\; M\,(\text{M}_\odot) \;\;|\;\;
    n\,(\text{cm}^{-3}) \;\;|\;\; N\,(1) \,\} \quad [Z\,(1)] \quad [T\,(\mathrm{K})] \f]

    Depending on the value of the \em massType option, the first column lists the average mass
    density \f$\rho\f$, the integrated mass \f$M\f$, the average number density \f$n\f$, or the
    integrated number density \f$N\f$ for the cell corresponding to the site. The precise units for
    this field are irrelevant because the total mass in the geometry will be normalized to unity
    after importing the data. However, the import procedure still insists on knowing the units.

    If the \em importMetallicity option is enabled, the next column specifies a "metallicity"
    fraction, which in this context is simply multiplied with the mass/density column to obtain the
    actual mass/density of the cell. If the \em importTemperature option is enabled, the next
    column specifies a temperature. If this temperature is higher than the maximum configured
    temperature, the mass and density for the site are set to zero, regardless of the mass or
    density specified in the fourth column. If the \em importTemperature option is disabled, or the
    maximum temperature value is set to zero, such a cutoff is not applied. */
class AdaptiveMeshGeometry : public MeshGeometry, public AdaptiveMeshInterface
{
    ITEM_CONCRETE(AdaptiveMeshGeometry, MeshGeometry,
                  "a geometry imported from data represented on an adaptive mesh (AMR grid)")
        ATTRIBUTE_TYPE_INSERT(AdaptiveMeshGeometry, "AdaptiveMeshInterface")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new AdaptiveMeshSnapshot object, calls its open() function,
        passes it the domain extent configured by the user, configures it to import a mass or a
        density column, and finally returns a pointer to the object. Ownership of the Snapshot
        object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;

    //=================== Other functions ==================

protected:
    /** This function implements the AdaptiveMeshInterface interface. It returns a pointer to the
        adaptive mesh snapshot maintained by this geometry. */
    AdaptiveMeshSnapshot* adaptiveMesh() const override;

    //===================== Data members ====================

private:
    // an extra pointer to our snapshot used to implement AdaptiveMeshInterface (ownership is passed to base class)
    AdaptiveMeshSnapshot* _adaptiveMeshSnapshot{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
