/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OLIGOWAVELENGTHDISTRIBUTION_HPP
#define OLIGOWAVELENGTHDISTRIBUTION_HPP

#include "Array.hpp"
#include "WavelengthDistribution.hpp"
class OligoWavelengthGrid;

//////////////////////////////////////////////////////////////////////

/** OligoWavelengthDistribution is a subclass of the WavelengthDistribution class representing the
    wavelength probability distribution of primary-source photon packets emitted in an
    oligochromatic simulation. In principle, the class implements a probability distribution that
    has a constant nonzero value inside a set of narrow ranges placed around the oligochromatic
    wavelengths in the simulation, corresponding to the bins defined by the oligo wavelength grid.
    These bins all have the same half bin width given by 1/1000 of the shortest wavelength in the
    list (see the OligoWavelengthGrid class). As a result, after normalization, the constant value
    of the probability distribution inside each of the bins is given by \f$1 / (N_\text{bins}
    \,w_\text{bin})\f$.

    In practice, in the interest of performance, rather than distributing the generated wavelengths
    across each bin, the class always generates exactly one of the discrete oligochromatic
    wavelengths. Also, the class always returns the normalized constant probability value
    regardless of the specified wavelength. This is acceptable because in an oligochromatic
    simulation this value is queried only for one of the discrete wavelengths anyway.

    The oligo wavelength distribution in a simulation is created programmatically rather than
    directly configured by the user. The list of wavelengths and the bin width are specified by
    passing an OligoWavelengthGrid instance to the constructor, which used the relevant
    information from this object. */
class OligoWavelengthDistribution : public WavelengthDistribution
{

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically to create an oligochromatic wavelength
        distribution. Before the constructor returns, the newly created object is hooked up as a
        child to the specified oligochromatic wavelength grid (so it will automatically be
        deleted), and its setup() function has been called. */
    explicit OligoWavelengthDistribution(OligoWavelengthGrid* wavelengthGrid);

    //======================== Other Functions =======================

public:
    /** This function returns the probability of the distribution at the given wavelength. For this
        class, the function always returns the normalized constant probability value \f$1 /
        (N_\text{bins} \,w_\text{bin})\f$ regardless of the specified wavelength. This is
        acceptable because in an oligochromatic simulation this function is called only for one of
        the discrete wavelengths anyway. */
    double probability(double wavelength) const override;

    /** This function draws a random wavelength from the wavelength distribution. For this class,
        it returns one of the discrete oligochromatic wavelengths with equal probability. */
    double generateWavelength() const override;

    //====================== Data members =====================

private:
    // data members initialized during construction
    const Array& _wavelengths;
    double _probability{0.};
};

//////////////////////////////////////////////////////////////////////

#endif
