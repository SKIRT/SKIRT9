/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIMESHGEOMETRY_HPP
#define VORONOIMESHGEOMETRY_HPP

#include "MeshGeometry.hpp"
#include "VoronoiMeshInterface.hpp"

////////////////////////////////////////////////////////////////////

/** A VoronoiMeshGeometry instance represents a 3D geometry with a spatial density distribution
    described by a list of sites generating a Voronoi tesselation of a cubodail domain. The data is
    usually extracted from a cosmological simulation snapshot, and it must be provided in a column
    text file formatted as described below. The total mass in the geometry is normalized to unity
    after importing the data.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    input file should contain 4, 5, or 6 columns, depending on the options configured by the user
    for this VoronoiMeshGeometry instance:

    \f[ x\,(\mathrm{pc}) \quad y\,(\mathrm{pc}) \quad z\,(\mathrm{pc}) \quad \{\,
    \rho\,(\text{M}_\odot\,\text{pc}^{-3}) \;\;|\;\; M\,(\text{M}_\odot) \;\;|\;\;
    n\,(\text{cm}^{-3}) \;\;|\;\; N\,(1) \,\} \quad [Z\,(1)] \quad [T\,(\mathrm{K})] \f]

    The first three columns are the \f$x\f$, \f$y\f$ and \f$z\f$ coordinates of the Voronoi site
    (i.e. the location defining a particular Voronoi cell).

    Depending on the value of the \em massType option, the fourth column lists the average mass
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
class VoronoiMeshGeometry : public MeshGeometry, public VoronoiMeshInterface
{
    ITEM_CONCRETE(VoronoiMeshGeometry, MeshGeometry, "a geometry imported from data represented on a Voronoi mesh")
        ATTRIBUTE_TYPE_INSERT(VoronoiMeshGeometry, "VoronoiMeshInterface")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new VoronoiMeshSnapshot object, calls its open() function,
        passes it the domain extent configured by the user, configures it to import a mass or a
        density column, and finally returns a pointer to the object. Ownership of the Snapshot
        object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;

    //=================== Other functions ==================

protected:
    /** This function implements the VoronoiMeshInterface interface. It returns a pointer to the
        Voronoi mesh snapshot maintained by this geometry. */
    VoronoiMeshSnapshot* voronoiMesh() const override;

    //===================== Data members ====================

private:
    // an extra pointer to our snapshot used to implement VoronoiMeshInterface (ownership is passed to base class)
    VoronoiMeshSnapshot* _voronoiMeshSnapshot{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
