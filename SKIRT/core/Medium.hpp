/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUM_HPP
#define MEDIUM_HPP

#include "Array.hpp"
#include "Position.hpp"
#include "SimulationItem.hpp"
class MaterialMix;

////////////////////////////////////////////////////////////////////

/** Medium is an abstract class representing a radiation transfer medium in the simulation. The
    configuration for a simulation may include multiple, spatially superimposed, transfer media.
    Each medium (an instance of a Medium subclass) must define the following information at each
    point in the spatial domain of the medium (which may differ from the size of the spatial grid
    in the simulation):
      - The material type and the optical material properties at a given wavelength.
      - Depending on the material type, the properties required to calculate medium state and
        secondary emission (for dust, calorimetric properties; for gas, kinetic temperature).
      - The number density and mass density of the material (implying a total number and mass).
      - The bulk velocity of the medium relative to the model coordinate frame.
      - If the material consists of non-spherical (spheroidal) particles, the magnetic field
        vector, which will be used to determine the direction and degree of alignment.

    The material properties at a particular spatial position in a medium are defined by an instance
    of a MaterialMix subclass. While the specific properties (such as, e.g., the scattering albedo)
    may vary with location, the material type and the level of support for various physical
    processes (such as, e.g., polarization) must be the same throughout the spatial domain (for a
    particular medium). Often, a Medium subclass employs just a single MaterialMix object to define
    material properties everywhere. In more complex situations, a family of MaterialMix objects of
    the same type may provided, for example, parameterized on a quantity imported from a
    hydrodynamical snapshot.

    Because the material must be of the same type throughout the spatial domain, a particular
    Medium instance has a well-defined material type (i.e. dust, electrons, or gas). It is not
    possible for a single Medium object to carry more than one material type. This limitation can
    be overcome by configuring multiple media, each with its own type, even if their spatial
    distribution happens to be defined by the same synthetic geometry or is imported from the same
    hydrodynamical snapshot.

    There is a complication for situations where the simulation is requested to self-consistently
    calculate the densities for multiple material types from the radiation field. For example, one
    could determine the amount of dust in a given location from the local gas state, including
    parameters such as the kinetic gas temperature or the degree of ionization. Because changing
    the dust density will influence the radiation field, this is clearly an iterative process, and
    it requires interaction with multiple Media objects. In this case, the density distributions in
    the media configured by the user merely serve as an initial state, which will be updated by the
    iterative state calculations.

    The spatial distribution of a medium may be represented internally using number density or mass
    density, at the choice of the implementation. Refer to the MaterialMix class header for more
    information on the conversion between number of entities and masses for each type of material.
*/
class Medium : public SimulationItem
{
    ITEM_ABSTRACT(Medium, SimulationItem, "a transfer medium")
    ITEM_END()

    //======================== Functions implemented by subclasses =======================

public:
    /** This function returns the dimension of the medium, which depends on the (lack of) symmetry
        of its geometry. A value of 1 means spherical symmetry, 2 means axial symmetry and 3 means
        none of these symmetries. */
    virtual int dimension() const = 0;

    /** This function returns (a pointer to) a MaterialMix object defining the material properties
        for the medium at the specified position. While the specific object being returned may vary
        with location, the material type and the level of support for various physical processes
        (such as, e.g., polarization) must be the same for all positions. An object of the
        appropriate type must be returned even if the density of the material happens to be zero at
        the specified position. */
    virtual const MaterialMix* mix(Position bfr) const = 0;

    /** This function returns (a pointer to) a default MaterialMix object representative of the
        material properties of this medium. In other words, it returns an arbitrary material mix
        with the same material type and level of support for various physical processes as any of
        the material mixes that may be returned by the mix(bfr) function. */
    virtual const MaterialMix* mix() const = 0;

    /** This function returns true if the mix(bfr) function for this medium may return a different
        MaterialMix object depending on the specified position, or false when the same object is
        always returned. */
    virtual bool hasVariableMix() const = 0;

    /** This function returns true if this medium may have a nonzero bulk velocity for some
        positions. It may be called before setup of the receiving medium has completed. */
    virtual bool hasVelocity() const = 0;

    /** This function returns the bulk velocity of the medium at the specified position. */
    virtual Vec bulkVelocity(Position bfr) const = 0;

    /** This function returns true if the magneticField() function for this medium may return a
        nonzero vector for some positions. */
    virtual bool hasMagneticField() const = 0;

    /** This function returns the magnetic field vector for the medium at the specified position.
        */
    virtual Vec magneticField(Position bfr) const = 0;

    /** This function returns true if the metallicity() function for this medium may return a
        nonzero value for some positions. */
    virtual bool hasMetallicity() const = 0;

    /** This function returns the metallicity of the medium at the specified position as defined in
        the input model, or zero if the input model does not define a metallicity for this medium
        (at all, or at the given position). */
    virtual double metallicity(Position bfr) const = 0;

    /** This function returns true if the temperature() function for this medium may return a
        nonzero value for some positions. */
    virtual bool hasTemperature() const = 0;

    /** This function returns the temperature of the medium at the specified position as defined in
        the input model, or zero if the input model does not define a temperature for this medium
        (at all, or at the given position). */
    virtual double temperature(Position bfr) const = 0;

    /** This function returns true if the parameters() function for this medium returns a nonempty
        array. */
    virtual bool hasParameters() const = 0;

    /** If custom input model parameters are available for this medium, this function stores the
        parameter values at the specified position into the given array. If the position is outside
        the domain, the parameter values default to zero. If no custom input model parameters are
        available for this medium, the array is resized to zero length. */
    virtual void parameters(Position bfr, Array& params) const = 0;

    /** This function returns the number density of the medium at the specified position. */
    virtual double numberDensity(Position bfr) const = 0;

    /** This function returns the total number of material entities in the medium. */
    virtual double number() const = 0;

    /** This function returns the mass density of the medium at the specified position. */
    virtual double massDensity(Position bfr) const = 0;

    /** This function returns the total mass in the medium. */
    virtual double mass() const = 0;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$
        along the full X axis of the model coordinate system. */
    virtual double opticalDepthX(double lambda) const = 0;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$
        along the full Y axis of the model coordinate system. */
    virtual double opticalDepthY(double lambda) const = 0;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$
        along the full Z axis of the model coordinate system. */
    virtual double opticalDepthZ(double lambda) const = 0;

    /** This function generates a random position sampled from the medium's spatial density
        distribution. It is undefined whether the function uses the number density or the mass
        density. If the conversion from number to mass is the same throughout the medium's spatial
        domain, there is no difference. In cases where it matters, the function will often use the
        more logical choice, depending on the way the density distribution was imported or
        normalized, but there is no guarantee. */
    virtual Position generatePosition() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
