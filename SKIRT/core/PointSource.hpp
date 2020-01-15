/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POINTSOURCE_HPP
#define POINTSOURCE_HPP

#include "AngularDistribution.hpp"
#include "PolarizationProfile.hpp"
#include "SpecialtySource.hpp"

//////////////////////////////////////////////////////////////////////

/** PointSource represents a primary radiation source that is limited to a single point in space.
    The spectral distribution is characterized by an SED object, and the bolometric output is
    characterized by a LuminosityNormalization object. The emitted radiation can be anisotropic
    (configured through an AngularDistribution object) and/or be polarized (configured through a
    PolarizationState object). The point source can also have a velocity. */
class PointSource : public SpecialtySource
{
    ITEM_CONCRETE(PointSource, SpecialtySource, "a primary point source")

        PROPERTY_DOUBLE(positionX, "the position of the point source, x component")
        ATTRIBUTE_QUANTITY(positionX, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionX, "0")
        ATTRIBUTE_INSERT(positionX, "positionX:Dimension3")

        PROPERTY_DOUBLE(positionY, "the position of the point source, y component")
        ATTRIBUTE_QUANTITY(positionY, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionY, "0")
        ATTRIBUTE_INSERT(positionY, "positionY:Dimension3")

        PROPERTY_DOUBLE(positionZ, "the position of the point source, z component")
        ATTRIBUTE_QUANTITY(positionZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionZ, "0")
        ATTRIBUTE_INSERT(positionZ, "positionZ:Dimension2")

        PROPERTY_ITEM(angularDistribution, AngularDistribution, "the angular luminosity distribution of the emission")
        ATTRIBUTE_DEFAULT_VALUE(angularDistribution, "IsotropicAngularDistribution")
        ATTRIBUTE_REQUIRED_IF(angularDistribution, "false")
        ATTRIBUTE_DISPLAYED_IF(angularDistribution, "Level2")

        PROPERTY_ITEM(polarizationProfile, PolarizationProfile, "the polarization profile of the emission")
        ATTRIBUTE_DEFAULT_VALUE(polarizationProfile, "NoPolarizationProfile")
        ATTRIBUTE_REQUIRED_IF(polarizationProfile, "false")
        ATTRIBUTE_DISPLAYED_IF(polarizationProfile, "Level2")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the point source, taking into account anisotropic
        emission or polarization, if present, but ignoring the velocity (because this is
        handled in the base class). */
    int geometryDimension() const override;

    /** This function causes the photon packet \em pp to be launched from the source using the
        given history index, wavelength, weighted luminosity contribution, and redshift interface.
        The position of the emission is determined randomly from the geometry configured for the
        source. The emission is unpolarized and isotropic; the emission direction is simply sampled
        from a uniform distribution on the unit sphere. */
    void launchSpecialty(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                         VelocityInterface* bvi) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
