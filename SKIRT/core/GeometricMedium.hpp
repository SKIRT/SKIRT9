/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GEOMETRICMEDIUM_HPP
#define GEOMETRICMEDIUM_HPP

#include "Medium.hpp"
#include "Geometry.hpp"
#include "MaterialNormalization.hpp"
#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** GeometricMedium represents a transfer medium for which the spatial material density
    distribution is characterized by a Geometry object. The material properties are characterized
    by a single MaterialMix object and are thus identical in all locations. The medium can have a
    single bulk velocity, i.e. the bulk velocity is identical in all locations. */
class GeometricMedium : public Medium
{
    ITEM_CONCRETE(GeometricMedium, Medium, "a transfer medium with a built-in geometry")

    PROPERTY_ITEM(geometry, Geometry, "the geometry of the spatial density distribution for the medium")
        ATTRIBUTE_DEFAULT_VALUE(geometry, "PlummerGeometry")

    PROPERTY_ITEM(materialMix, MaterialMix, "the material type and properties throughout the medium")
        ATTRIBUTE_DEFAULT_VALUE(materialMix, "MeanTrustBenchmarkDustMix")

    PROPERTY_ITEM(normalization, MaterialNormalization, "the type of normalization for the amount of material")
        ATTRIBUTE_DEFAULT_VALUE(normalization, "OpticalDepthMaterialNormalization")

    PROPERTY_DOUBLE(velocityX, "the bulk velocity of the medium, x component")
        ATTRIBUTE_QUANTITY(velocityX, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityX, "[0")
        ATTRIBUTE_MAX_VALUE(velocityX, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityX, "0")

    PROPERTY_DOUBLE(velocityY, "the bulk velocity of the medium, y component")
        ATTRIBUTE_QUANTITY(velocityY, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityY, "[0")
        ATTRIBUTE_MAX_VALUE(velocityY, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityY, "0")

    PROPERTY_DOUBLE(velocityZ, "the bulk velocity of the medium, z component")
        ATTRIBUTE_QUANTITY(velocityZ, "velocity")
        ATTRIBUTE_MIN_VALUE(velocityZ, "[0")
        ATTRIBUTE_MAX_VALUE(velocityZ, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(velocityZ, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function determines the normalization for number and mass density in this medium.
        Given that a geometry has a total mass of one by definition, the normalization factors
        correspond to the total number and the total mass of the medium. */
    void setupSelfAfter() override;

    //======================== Functions implemented by subclasses =======================

    /** This function returns the dimension of the medium, which is the same as the dimension of
        its spatial distribution, except if there is a nonzero bulk velocity. */
    int dimension() const override;

    /** This function returns the MaterialMix object defining the material properties for the
        medium. The same object is returned regardless of position. */
    const MaterialMix* mix(Position bfr) const override;

    /** This function returns the number density of the medium at the specified position. */
    double numberDensity(Position bfr) const override;

    /** This function returns the total number of material entities in the medium. */
    double number() const override;

    /** This function returns the mass density of the medium at the specified position. */
    double massDensity(Position bfr) const override;

    /** This function returns the total mass in the medium. */
    double mass() const override;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$
        along the full X axis of the model coordinate system. */
    double opticalDepthX(double lambda) const override;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$
        along the full Y axis of the model coordinate system. */
    double opticalDepthY(double lambda) const override;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$
        along the full Z axis of the model coordinate system. */
    double opticalDepthZ(double lambda) const override;

    /** This function generates a random position sampled from the medium's spatial density
        distribution. Because the conversion from number to mass is the same throughout the
        medium's spatial domain, there is no difference between sampling from the number density
        or the mass density. */
    Position generatePosition() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    double _number{0.};     // the total number of entities in the medium
    double _mass{0.};       // the total mass in the medium
};

////////////////////////////////////////////////////////////////////

#endif
