/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GEOMETRICSOURCE_HPP
#define GEOMETRICSOURCE_HPP

#include "NormalizedSource.hpp"
#include "Geometry.hpp"

//////////////////////////////////////////////////////////////////////

/** GeometricSource represents a primary radiation source for which the spatial luminosity
    distribution is characterized by a Geometry object. The spectral distribution is identical in
    all locations and is characterized by an SED object. The bolometric output is characterized by
    a LuminosityNormalization object. The emitted radiation is isotropic and unpolarized. The
    source can have a single bulk velocity, i.e. the bulk velocity is identical in all locations.
    */
class GeometricSource : public NormalizedSource
{
    ITEM_CONCRETE(GeometricSource, NormalizedSource, "a primary source with a built-in geometry")

    PROPERTY_ITEM(geometry, Geometry, "the geometry of the spatial luminosity distribution for the source")
        ATTRIBUTE_DEFAULT_VALUE(geometry, "PlummerGeometry")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry configured for the source. */
    int geometryDimension() const override;

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index, wavelength, weighted luminosity contribution, and redshift interface.
        The position of the emission is determined randomly from the geometry configured for the
        source. The emission is unpolarized and isotropic; the emission direction is simply sampled
        from a uniform distribution on the unit sphere. */
    void launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                          VelocityInterface* bvi) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
