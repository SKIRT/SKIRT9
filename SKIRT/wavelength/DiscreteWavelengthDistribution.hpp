/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DISCRETEWAVELENGTHDISTRIBUTION_HPP
#define DISCRETEWAVELENGTHDISTRIBUTION_HPP

#include "GridWavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** DiscreteWavelengthDistribution is a specialty class for representing discrete wavelength
    probability distributions derived from a wavelength grid that can be configured by the user.
    This class is useful for mimicking codes that emit photon packets at discrete wavelengths
    rather than across a continuous range, perhaps to compare results in the context of a
    benchmark. One has to take great care when using this class in a simulation configuration, as
    indicated in more detail below.

    When a DiscreteWavelengthDistribution instance is used as the wavelength bias distribution for
    a source with a composite bias factor of one, the source will emit photon packets only at the
    characteristic wavelengths of the configured grid, and the photon packets will be distributed
    with equal probability among those wavelengths. More precisely, the
    DiscreteWavelengthDistribution instance uses only the characteristic wavelengths of the
    configured grid that fall inside the wavelength range of the associated source (obtained
    through the SourceWavelengthRangeInterface). If none of the characteristic wavelengths fall
    inside the source range, a fatal error is issued.

    When sending photon packets at a particular grid of discrete wavelengths, it makes little sense
    to configure wavelength grids with different bins for detecting photon packets in other areas
    of the simulation. Doing so might cause unexpected results. For example, some bins might
    receive no photon packets at all and thus incorrectly report zero influx. It is therefore best
    to configure the same wavelength grid throughout the simulation.

    To configure a simulation that uses discrete wavelengths, follow these guidelines:

    - Choose a wavelength grid with adjacent bins (i.e. representing a continuous range), such as a
    LogWavelengthGrid or a FileWavelengthGrid with \em relativeHalfWidth set to zero.

    - Configure all primary and secondary sources (including, e.g. dust emission) with a wavelength
    bias factor of one, and a wavelength bias distribution of type DiscreteWavelengthDistribution
    using the chosen wavelength grid.

    - Explicitly configure all wavelength grids in the simulation to be the same chosen wavelength
    grid. This includes the radiation field wavelength grid, the dust emission wavelength grid, and
    the wavelength grids for all instruments and probes.

    - Avoid using kinematics in the simulation. Small Doppler shifts from the emitted discrete
    wavelengths will stay within the same wavelength bin anyway. Larger shifts will cause photon
    packets to move between wavelength bins in a discontinuous way.

    See also BinnedWavelengthDistribution, which offers similar functionality but spreads the
    generated wavelengths across the full width of each bin instead of concentrating them at the
    bin's characteristic wavelength. */
class DiscreteWavelengthDistribution : public GridWavelengthDistribution
{
    ITEM_CONCRETE(DiscreteWavelengthDistribution, GridWavelengthDistribution,
                  "a discrete wavelength probability distribution derived from a wavelength grid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function returns the characteristic wavelength of the given bin, so that the
        distribution implements a Dirac-delta function within each of the (in-range) bins of the
        configured wavelength grid. */
    double wavelengthInBin(int ell) const override;
};

////////////////////////////////////////////////////////////////////

#endif
