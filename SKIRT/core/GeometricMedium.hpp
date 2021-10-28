/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GEOMETRICMEDIUM_HPP
#define GEOMETRICMEDIUM_HPP

#include "Geometry.hpp"
#include "MaterialMix.hpp"
#include "MaterialNormalization.hpp"
#include "Medium.hpp"
#include "VectorField.hpp"

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
        ATTRIBUTE_DEFAULT_VALUE(materialMix, "MeanInterstellarDustMix")

        PROPERTY_ITEM(normalization, MaterialNormalization, "the type of normalization for the amount of material")
        ATTRIBUTE_DEFAULT_VALUE(normalization, "OpticalDepthMaterialNormalization")

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
        ATTRIBUTE_INSERT(velocityMagnitude, "velocityDistribution&velocityMagnitude:MediumVelocity")

        PROPERTY_ITEM(magneticFieldDistribution, VectorField, "the spatial distribution of the magnetic field, if any")
        ATTRIBUTE_REQUIRED_IF(magneticFieldDistribution, "false")
        ATTRIBUTE_DISPLAYED_IF(magneticFieldDistribution, "Level3")

        PROPERTY_DOUBLE(magneticFieldStrength, "the strength of the magnetic field  (multiplier)")
        ATTRIBUTE_QUANTITY(magneticFieldStrength, "magneticfield")
        ATTRIBUTE_MIN_VALUE(magneticFieldStrength, "[-1 T")
        ATTRIBUTE_MAX_VALUE(magneticFieldStrength, "1 T]")
        ATTRIBUTE_RELEVANT_IF(magneticFieldStrength, "magneticFieldDistribution")
        ATTRIBUTE_DISPLAYED_IF(magneticFieldStrength, "Level3")
        ATTRIBUTE_INSERT(magneticFieldStrength, "magneticFieldDistribution&magneticFieldStrength:MagneticField")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function determines the normalization for number and mass density in this medium.
        Given that a geometry has a total mass of one by definition, the normalization factors
        correspond to the total number and the total mass of the medium. */
    void setupSelfAfter() override;

    //======================== Functions implemented by subclasses =======================

    /** This function returns the dimension of the medium, which is determined by the symmetries of
        the spatial distribution of its density, bulk velocity (if any) and magnetic field (if
        any). */
    int dimension() const override;

    /** This function returns the MaterialMix object defining the material properties for the
        medium. The same object is returned regardless of position. */
    const MaterialMix* mix(Position bfr) const override;

    /** This function returns the MaterialMix object defining the material properties for the
        medium. */
    const MaterialMix* mix() const override;

    /** This function always returns false because a geometric medium has the same dust mix
        throughout its spatial domain. */
    bool hasVariableMix() const override;

    /** This function returns true if the bulk velocity of the medium is nonzero. */
    bool hasVelocity() const override;

    /** This function returns the bulk velocity of the medium. The same velocity is returned
        regardless of position. */
    Vec bulkVelocity(Position bfr) const override;

    /** This function returns true if the medium has a nonempty \em magneticFieldDistribution and a
        nonzero \em magneticFieldStrength. */
    bool hasMagneticField() const override;

    /** This function returns the magnetic field vector of the medium at the specified position. If
        hasMagneticField() returns true, this function returns the product of the configured
        magnetic field strength and the vector retrieved from the configured normalized vector
        field for the given position; otherwise it returns a zero magnetic field. */
    Vec magneticField(Position bfr) const override;

    /** This function returns false because geometric media do not define a metallicity. */
    bool hasMetallicity() const override;

    /** This function always returns zero because geometric media do not define a metallicity */
    double metallicity(Position bfr) const override;

    /** This function returns false because geometric media do not define a temperature. */
    bool hasTemperature() const override;

    /** This function always returns zero because geometric media do not define a temperature. */
    double temperature(Position bfr) const override;

    /** This function returns false because geometric media do not define custom parameters. */
    bool hasParameters() const override;

    /** This function resizes the specified array to zero length because a geometric medium does
        not provide custom input model parameters. */
    void parameters(Position bfr, Array& params) const override;

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
    double _number{0.};  // the total number of entities in the medium
    double _mass{0.};    // the total mass in the medium
};

////////////////////////////////////////////////////////////////////

#endif
