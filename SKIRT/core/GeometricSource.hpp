/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GEOMETRICSOURCE_HPP
#define GEOMETRICSOURCE_HPP

#include "Geometry.hpp"
#include "NormalizedSource.hpp"
#include "VelocityInterface.hpp"

//////////////////////////////////////////////////////////////////////

/** GeometricSource represents a primary radiation source for which the spatial luminosity
    distribution is characterized by a Geometry object. The spectral distribution is identical in
    all locations and is characterized by an SED object. The bolometric output is characterized by
    a LuminosityNormalization object. The emitted radiation is isotropic and unpolarized. The
    source can have a single bulk velocity, i.e. the bulk velocity is identical in all locations.
    */
class GeometricSource : public NormalizedSource, public VelocityInterface
{
    ITEM_CONCRETE(GeometricSource, NormalizedSource, "a primary source with a built-in geometry")

        PROPERTY_ITEM(geometry, Geometry, "the geometry of the spatial luminosity distribution for the source")
        ATTRIBUTE_DEFAULT_VALUE(geometry, "PlummerGeometry")

        PROPERTY_DOUBLE(velocityX, "the bulk velocity of the source, x component")
        ATTRIBUTE_QUANTITY(velocityX, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityX, "[0")
        ATTRIBUTE_MAX_VALUE(velocityX, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityX, "0")
        ATTRIBUTE_RELEVANT_IF(velocityX, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityX, "Level2")
        ATTRIBUTE_INSERT(velocityX, "Panchromatic&velocityX:Dimension3")

        PROPERTY_DOUBLE(velocityY, "the bulk velocity of the source, y component")
        ATTRIBUTE_QUANTITY(velocityY, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityY, "[0")
        ATTRIBUTE_MAX_VALUE(velocityY, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityY, "0")
        ATTRIBUTE_RELEVANT_IF(velocityY, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityY, "Level2")
        ATTRIBUTE_INSERT(velocityY, "Panchromatic&velocityY:Dimension3")

        PROPERTY_DOUBLE(velocityZ, "the bulk velocity of the source, z component")
        ATTRIBUTE_QUANTITY(velocityZ, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityZ, "[0")
        ATTRIBUTE_MAX_VALUE(velocityZ, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityZ, "0")
        ATTRIBUTE_RELEVANT_IF(velocityZ, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(velocityZ, "Level2")
        ATTRIBUTE_INSERT(velocityZ, "Panchromatic&velocityZ:Dimension2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function creates a private object offering the velocity interface if the bulk velocity
        is nonzero. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the source, which is the same as the dimension of
        its spatial distribution (to be provided by the subclass), except if there is a nonzero
        bulk velocity. */
    int dimension() const override;

    /** This function returns true if the velocity of the source is nonzero. */
    bool hasVelocity() const override;

    /** This function implements the VelocityInterface interface. It returns the bulk velocity
         of this source, as configured by the user. */
    Vec velocity() const override;

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index, wavelength, weighted luminosity contribution, and redshift interface.
        The position of the emission is determined randomly from the geometry configured for the
        source. The emission is unpolarized and isotropic; the emission direction is simply sampled
        from a uniform distribution on the unit sphere. */
    void launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw) const override;

    //======================== Data Members ========================

private:
    // pointer to an object offering the redshift interface; either "this" or null pointer if the bulk velocity is zero
    VelocityInterface* _bvi{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif
