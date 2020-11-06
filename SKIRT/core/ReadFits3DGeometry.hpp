/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef READFITS3DGEOMETRY_HPP
#define READFITS3DGEOMETRY_HPP

#include "Array.hpp"
#include "GenGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The ReadFits3DGeometry class is a subclass of the GenGeometry class, and describes an arbitary
    3D geometry for a single component, read from a 3D (.fits file) datacube. The model geometry is
    set by two parameters: the input filename and the pixel scale. */
class ReadFits3DGeometry : public GenGeometry
{
    ITEM_CONCRETE(ReadFits3DGeometry, GenGeometry, "a 3D geometry read from a FITS file")
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
    double _xmin{0.}, _xmax{0.}, _ymin{0.}, _ymax{0.}, _zmin{0.}, _zmax{0.};
};

////////////////////////////////////////////////////////////////////

#endif
