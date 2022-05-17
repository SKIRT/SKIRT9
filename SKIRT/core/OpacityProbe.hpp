/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OPACITYPROBE_HPP
#define OPACITYPROBE_HPP

#include "MaterialWavelengthRangeInterface.hpp"
#include "SpatialGridWhenFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** OpacityProbe probes the opacity of the medium at a given wavelength, as discretized on the
    spatial grid of the simulation. In the context of this probe, the opacity \f$k\f$ of a medium
    (at a given position and wavelength) is defined as the product of the local density and
    extinction cross section, \f$k=n\varsigma=\rho\kappa\f$. The probe can be used with any Form
    subclass. When associated with a form that projects the quantity along a path \f$s\f$, the
    probe outputs optical depth, \f$\tau=\int_s k(s) \mathrm{d}s\f$.

    Because probing is performed without the context of a photon packet, default values are used
    for any relevant incoming photon packet properties. For example, the effects of kinematics are
    ignored and the radiation is assumed to be unpolarized.

    The user can select the aggregation level, i.e. whether to produce an output file per medium
    component, per medium type (dust, electrons, gas), or for the complete medium system. There is
    also an option to decide whether the probe should be performed after setup or after the full
    simulation run. The latter option is meaningful if the density and/or cross section of the
    media may change during the simulation.

    This probe implements the MaterialWavelengthRangeInterface to indicate that
    wavelength-dependent material properties may be required for the configured wavelength. */
class OpacityProbe : public SpatialGridWhenFormProbe, public MaterialWavelengthRangeInterface
{
    /** The enumeration type indicating how to aggregate the output: per medium component, per
        medium type (dust, electrons, gas), or for the complete medium system. */
    ENUM_DEF(Aggregation, Component, Type, System)
        ENUM_VAL(Aggregation, Component, "per medium component")
        ENUM_VAL(Aggregation, Type, "per medium type (dust, electrons, gas)")
        ENUM_VAL(Aggregation, System, "for the complete medium system")
    ENUM_END()

    ITEM_CONCRETE(OpacityProbe, SpatialGridWhenFormProbe, "internal spatial grid: opacity of the medium")
        ATTRIBUTE_TYPE_DISPLAYED_IF(OpacityProbe, "Medium&SpatialGrid")

        PROPERTY_DOUBLE(wavelength, "the wavelength at which to determine the opacity")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "0.55 micron")
        ATTRIBUTE_DISPLAYED_IF(wavelength, "Level2")

        PROPERTY_ENUM(aggregation, Aggregation, "how to aggregate the opacity")
        ATTRIBUTE_DEFAULT_VALUE(aggregation, "Type")
        ATTRIBUTE_DISPLAYED_IF(aggregation, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function returns a wavelength range corresponding to the user-configured wavelength,
        indicating that wavelength-dependent material properties may be required for this
        wavelength. */
    Range wavelengthRange() const override;

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
