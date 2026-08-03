/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GEOMETRICSOURCE_HPP
#define GEOMETRICSOURCE_HPP

#include "Geometry.hpp"
#include "NormalizedSource.hpp"
#include "VectorField.hpp"

//////////////////////////////////////////////////////////////////////

/** GeometricSource represents a primary radiation source for which the spatial luminosity
    distribution is characterized by a Geometry object. The spectral distribution is identical in
    all locations and is characterized by an SED object. The bolometric output is characterized by
    a LuminosityNormalization object. The emitted radiation is isotropic and unpolarized. The
    source can also be assigned a vector field determining the bulk velocity in each location. */
class GeometricSource : public NormalizedSource
{
    ITEM_CONCRETE(GeometricSource, NormalizedSource, "a primary source with a built-in geometry")

        PROPERTY_ITEM(geometry, Geometry, "the geometry of the spatial luminosity distribution for the source")
        ATTRIBUTE_DEFAULT_VALUE(geometry, "PlummerGeometry")

        PROPERTY_ITEM(velocityDistribution, VectorField, "the spatial distribution of the velocity, if any")
        ATTRIBUTE_REQUIRED_IF(velocityDistribution, "false")
        ATTRIBUTE_RELEVANT_IF(velocityDistribution, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityDistribution, "Level2")

        PROPERTY_DOUBLE(velocityMagnitude, "the magnitude of the velocity (multiplier)")
        ATTRIBUTE_QUANTITY(velocityMagnitude, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityMagnitude, "[-100000 km/s")
        ATTRIBUTE_MAX_VALUE(velocityMagnitude, "100000 km/s]")
        ATTRIBUTE_RELEVANT_IF(velocityMagnitude, "Panchromatic&velocityDistribution")
        ATTRIBUTE_DISPLAYED_IF(velocityMagnitude, "Level2")
        ATTRIBUTE_INSERT(velocityMagnitude, "velocityDistribution&velocityMagnitude:SourceVelocity")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches a flag indicating whether there is a nonzero velocity field. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the source, which is the same as the dimension of
        its spatial distribution combined with the dimension of the velocity field, if there is
        one. */
    int dimension() const override;

    /** This function returns true if the source has a nonzero velocity field. */
    bool hasVelocity() const override;

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index, wavelength, weighted luminosity contribution, and redshift interface.
        The position of the emission is determined randomly from the geometry configured for the
        source, and the bulk velocity is derived from the assigned velocity field, if there is one.
        The emission is unpolarized and isotropic; the emission direction is simply sampled from a
        uniform distribution on the unit sphere. */
    void launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw) const override;

    //======================== Data Members ========================

private:
    bool _hasVelocity{false};
};

//////////////////////////////////////////////////////////////////////

#endif
