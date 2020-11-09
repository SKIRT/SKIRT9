/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef READFITS3DGEOMETRY_HPP
#define READFITS3DGEOMETRY_HPP

#include "Array.hpp"
#include "GenGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The ReadFits3DGeometry class is a subclass of the GenGeometry class, and describes an arbitrary
    3D geometry for a single component read from a FITS file datacube. The model geometry is
    defined by two parameters: the input filename and the pixel scale (i.e. the physical length per
    pixel). The geometry is automatically centered on the origin.

    The input geometry should be provided as a 3D datacube, i.e. an ndarray with shape (nx, ny,
    nz), stored in an HDU of the FITS input file. The order of axes is x, y, z and the header
    corresponding to this data array resembles:

        NAXIS   =                    3 / number of array dimensions
        NAXIS1  =                   nx
        NAXIS2  =                   ny
        NAXIS3  =                   nz

    with nx, ny and nz the number of cells in each spatial dimension. The input filename should
    include the ".fits" filename extension. By default, the first HDU with nonempty data is used.
    Using the data from another HDU in the provided file is possible by specifying this HDU name
    between square brackets after the input filename, for example: "geometry.fits[USETHISHDU]". */
class ReadFits3DGeometry : public GenGeometry
{
    ITEM_CONCRETE(ReadFits3DGeometry, GenGeometry, "a geometry read from a 3D FITS file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ReadFits3DGeometry, "Level2")

        PROPERTY_STRING(filename, "the filename of the datacube")

        PROPERTY_DOUBLE(pixelScale, "the pixel scale for the datacube (i.e. the physical length per pixel)")
        ATTRIBUTE_QUANTITY(pixelScale, "length")
        ATTRIBUTE_MIN_VALUE(pixelScale, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads in the .fits datacube. A vector with the normalized cumulative
        distribution in each spaxel is computed, satisfying the condition that the total mass
        equals 1. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(x,y,z)\f$ at the position (x,y,z). */
    double density(Position bfr) const override;

    /** This function generates a random position \f$(x,y,z)\f$ from the geometry, by drawing a
        random point from the appropriate probability density distribution function. The
        \f$(x,y,z)\f$ coordinates are derived from the normalized cumulative distribution vector.
        */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the density along
        the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0) \, \text{d} x \f] */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density along
        the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0) \, \text{d} y \f] */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z) \, \text{d} z \f] */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // the input datacube and the cumulative distribution
    Array _datacube;
    int _nx{0}, _ny{0}, _nz{0};
    Array _Xv;

    // other data members initialized during setup
    double _xmin{0.}, _ymin{0.}, _zmin{0.};
};

////////////////////////////////////////////////////////////////////

#endif
