/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHDISTRIBUTION_HPP
#define WAVELENGTHDISTRIBUTION_HPP

#include "SimulationItem.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** An instance of a WavelengthDistribution subclass represents a plain wavelength probability
    distribution. It can used, for example, to define a weighting scheme for the wavelengths
    assigned to photon packets launched from a source. This abstract base class just defines an
    interface that must be implemented by each subclass. There are two key operations: drawing a
    random wavelength from the probability distribution, and returning the probability for a given
    wavelength.

    The wavelength probability distribution is automatically normalized to unity over the
    wavelength range of the associated source (obtained through the SourceWavelengthRangeInterface),
    intersected with the intrinsic wavelength range of the distribution. Consequently, the random
    wavelengths returned by the generateWavelength() function will always fall inside the
    intersected range, and the value returned by the probability() function can be nonzero only
    within that same range. */
class WavelengthDistribution : public SimulationItem
{
    ITEM_ABSTRACT(WavelengthDistribution, SimulationItem, "a wavelength probability distribution")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the probability of the distribution at the given wavelength, or zero
        if the wavelength is out of range (the intersection of the intrinsic and external
        wavelength ranges). */
    virtual double probability(double wavelength) const = 0;

    /** This function draws a random wavelength from the wavelength distribution. The returned
        value will always fall inside the intersection of the intrinsic and external wavelength
        ranges. */
    virtual double generateWavelength() const = 0;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data member initialized during setup
    Random* _random{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
