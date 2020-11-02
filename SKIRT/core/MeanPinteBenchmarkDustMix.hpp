/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANPINTEBENCHMARKDUSTMIX_HPP
#define MEANPINTEBENCHMARKDUSTMIX_HPP

#include "SingleGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MeanPinteBenchmarkDustMix class represents a population of identical dust grains with a
    size of 1 micron used for the benchmark described by Pinte et al. 2009, A&A, 498, 967P. It
    includes scattering polarization data, assuming spherical grains, only at a single wavelength
    of 1 micron.

    The web site for the Pinte et al. 2009 benchmark is:
    http://ipag.osug.fr/~pintec/benchmark/index.shtml

    The basic optical properties for the dust mixture were downloaded from:
    http://ipag.osug.fr/~pintec/benchmark/data/opacity.dat

    The Mueller matrix coefficients for the dust mixture at 1 micron were downloaded from:
    http://ipag.osug.fr/~pintec/benchmark/data/mueller_matrix.txt
    */
class MeanPinteBenchmarkDustMix : public SingleGrainDustMix
{
    /** The enumeration type indicating the scattering mode. */
    ENUM_DEF(ScatteringType, HenyeyGreenstein, MaterialPhaseFunction, SphericalPolarization)
        ENUM_VAL(ScatteringType, HenyeyGreenstein, "use the Henyey-Greenstein phase function (unpolarized)")
        ENUM_VAL(ScatteringType, MaterialPhaseFunction,
                 "use the phase function derived from actual material properties (unpolarized)")
        ENUM_VAL(ScatteringType, SphericalPolarization, "support polarization through scattering by spherical grains")
    ENUM_END()

    ITEM_CONCRETE(MeanPinteBenchmarkDustMix, SingleGrainDustMix,
                  "a Pinte 2D benchmark dust mix (mean properties, optionally with polarization at 1 micron)")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MeanPinteBenchmarkDustMix, "Level2")

        PROPERTY_ENUM(scatteringType, ScatteringType, "the type of scattering to be implemented")
        ATTRIBUTE_DEFAULT_VALUE(scatteringType, "HenyeyGreenstein")

    ITEM_END()

    //======================== Capabilities =======================

public:
    /** This function returns the scattering mode supported by this material mix as configured by
        the user through the scatteringType property. */
    ScatteringMode scatteringMode() const override;

    /** This function returns a flag indicating whether the material mix supports polarization
        during scattering events or not. For this dust mix, the function returns true if the user
        configured the SphericalPolarization scattering type, and false otherwise. */
    bool hasPolarizedScattering() const override;

    //======================== Other Functions =======================

protected:
    /** This function returns the name of the stored table resource tabulating the basic optical
        properties for this dust mix. */
    string resourceNameForOpticalProps() const override;

    /** This function returns the name of the stored table resource tabulating the Mueller matrix
        elements for this dust mix. */
    string resourceNameForMuellerMatrix() const override;
};

////////////////////////////////////////////////////////////////////

#endif
