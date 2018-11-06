/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEDISCRETEWAVELENGTHDISTRIBUTION_HPP
#define FILEDISCRETEWAVELENGTHDISTRIBUTION_HPP

#include "DiscreteWavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** A FileDiscreteWavelengthDistribution object represents a discrete wavelength probability
    distribution that is loaded from an input file. The floating point numbers in the first column
    of the text file specify the wavelengths to be included in the discrete distribution. The
    wavelengths are by default given in micron (the units can be overridden by column header info
    in the file) and they can be listed in arbitrary order (they are sorted automatically). Any
    wavelengths outside of the wavelength range of the associated source are automatically removed
    from the list before the distribution is constructed. If no wavelengths are inside the source
    range, a fatal error is issued.

    For more information on the implemented "discrete" distribution, refer to the documentation of
    the DiscreteWavelengthDistribution class. */
class FileDiscreteWavelengthDistribution : public DiscreteWavelengthDistribution
{
    ITEM_CONCRETE(FileDiscreteWavelengthDistribution, DiscreteWavelengthDistribution,
                  "a discrete wavelength probability distribution loaded from a text file")

    PROPERTY_STRING(filename, "the name of the file with the wavelengths")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads and returns the tabulated wavelengths from the input file. */
    Array getWavelengths() const override;
};

////////////////////////////////////////////////////////////////////

#endif
