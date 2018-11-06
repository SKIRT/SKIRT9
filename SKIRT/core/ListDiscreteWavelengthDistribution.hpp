/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTDISCRETEWAVELENGTHDISTRIBUTION_HPP
#define LISTDISCRETEWAVELENGTHDISTRIBUTION_HPP

#include "DiscreteWavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** A ListDiscreteWavelengthDistribution object represents a discrete wavelength probability
    distribution that is fully specified inside the configuration file (i.e. without referring to
    an input file) by a list of wavelengths. The wavelengths can be listed in arbitrary order (they
    are sorted automatically), and any wavelengths outside of the wavelength range of the
    associated source are automatically removed from the list before the distribution is
    constructed. If no wavelengths are inside the source range, a fatal error is issued.

    For more information on the implemented "discrete" distribution, refer to the documentation of
    the DiscreteWavelengthDistribution class. */
class ListDiscreteWavelengthDistribution : public DiscreteWavelengthDistribution
{
    ITEM_CONCRETE(ListDiscreteWavelengthDistribution, DiscreteWavelengthDistribution,
                  "a discrete wavelength probability distribution specified inside the configuration file")

    PROPERTY_DOUBLE_LIST(wavelengths, "the wavelengths")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns the configured wavelengths. */
    Array getWavelengths() const override;
};

////////////////////////////////////////////////////////////////////

#endif
