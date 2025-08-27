/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OPACITYPROBE_HPP
#define OPACITYPROBE_HPP

#include "DisjointWavelengthGrid.hpp"
#include "MaterialWavelengthRangeInterface.hpp"
#include "SpatialGridWhenFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** OpacityProbe probes the opacity of the medium on a given wavelength grid, as discretized on the
    spatial grid of the simulation. In the context of this probe, the opacity \f$k\f$ of a medium
    (at a given position and wavelength) is defined as the product of the local density and
    extinction cross section, \f$k=n\varsigma=\rho\kappa\f$. The probe can be used with any Form
    subclass. When associated with a form that projects the quantity along a path \f$s\f$, the
    probe outputs optical depth, \f$\tau=\int_s k(s) \mathrm{d}s\f$.

    The resulting opacities are evaluated at the characteristic wavelengths of the provided 
    wavelength grid. If none is provided, the default wavelength grid is used.

    Because probing is performed without the context of a photon packet, default values are used
    for any relevant incoming photon packet properties. For example, the effects of kinematics are
    ignored and the radiation is assumed to be unpolarized.

    The user can select the aggregation level, i.e. whether to produce an output file per medium
    component, per medium type (dust, electrons, gas), or for the complete medium system. If one or
    more medium components in the simulation are equipped with a FragmentDustMixDecorator, the
    probe can provide information for each of the dust grain populations represented by the
    decorator. Depending on the value of the \em fragmentSizeBins flag on the decorator, there are
    fragments for each of the grain material types or even for each of the grain size bins defined
    by the underlying dust mixture. The probed information is written in a separate file for each
    fragment, identified by a zero-based fragment index in addition to the zero-based component
    index.

    There is also an option to decide whether the probe should be performed after setup or after
    the full simulation run. The latter option is meaningful if the density and/or cross section of
    the media may change during the simulation.

    This probe implements the MaterialWavelengthRangeInterface to indicate that
    wavelength-dependent material properties may be required for the configured wavelength. */
class OpacityProbe : public SpatialGridWhenFormProbe, public MaterialWavelengthRangeInterface
{
    /** The enumeration type indicating how to aggregate the output: per medium component, per
        medium type (dust, electrons, gas), or for the complete medium system. */
    ENUM_DEF(Aggregation, Fragment, Component, Type, System)
        ENUM_VAL(Aggregation, Fragment, "per fragment (dust grain material type and/or size bin)")
        ENUM_VAL(Aggregation, Component, "per medium component")
        ENUM_VAL(Aggregation, Type, "per medium type (dust, electrons, gas)")
        ENUM_VAL(Aggregation, System, "for the complete medium system")
    ENUM_END()

    ITEM_CONCRETE(OpacityProbe, SpatialGridWhenFormProbe, "internal spatial grid: opacity of the medium")
        ATTRIBUTE_TYPE_DISPLAYED_IF(OpacityProbe, "Medium&SpatialGrid")

        PROPERTY_ITEM(wavelengthGrid, DisjointWavelengthGrid, "the wavelength grid for this probe")
        ATTRIBUTE_RELEVANT_IF(wavelengthGrid, "Panchromatic")
        ATTRIBUTE_REQUIRED_IF(wavelengthGrid, "!DefaultInstrumentWavelengthGrid")
        ATTRIBUTE_DISPLAYED_IF(wavelengthGrid, "Level2")

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

    /** This function returns a pointer to the user-configured wavelength grid for this probe, if
        any, indicating that wavelength-dependent material properties may be required for these
        wavelengths. */
    WavelengthGrid* materialWavelengthGrid() const override;

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
