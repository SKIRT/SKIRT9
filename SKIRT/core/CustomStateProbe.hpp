/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CUSTOMSTATEPROBE_HPP
#define CUSTOMSTATEPROBE_HPP

#include "SpatialGridWhenFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** CustomStateProbe probes the custom medium state quantities for medium components configured
    with a material mix that requests such variables. Examples include the dust fragment weights
    used by the FragmentDustMixDecorator, the neutral and atomic hydrogen fractions required by the
    SpinFlipHydrogenGasMix, and the energy level population densities used by the
    MolecularLineGasMix. These custom quantities are always discretized on the simulation's
    spatial grid. The probe outputs a compound quantity with a value for each of the custom medium
    state quantities. It can be used with any Form subclass, taking into account the caveat
    discussed below. When associated with a form that projects the quantity along a path, the value
    is density-weighted.

    A given material mix (and thus its associated medium component) may use custom state quantities
    of different types and hence different output units. This is no problem for text column output
    files because the quantity types and units of the different columns are listed in the file
    header. The situation is somewhat different for FITS output files. Each quantity is stored in
    its own frame ("image") of the data cube using the appropriate units, but there is no mechanism
    to advertise the units for each of the frames. The file header lists just the units of the
    quantity stored in the first frame.

    The \em indices option allows specifying a subset of the custom state quantities to be probed.
    This allows restricting the probe to quantities with the same units or otherwise grouping or
    reordering quantities. If the \em indices string is empty, all custom quantities are included.
    Otherwise, the probe expects a comma-separated list of zero-based indices or index ranges
    indicated with a dash. Examples of valid strings include: "5,1,2" or "0-9,17". Incorrectly
    formatted string segments and out-of-range indices are ignored.

    The probe produces output for each medium component that has custom quantities, or if the \em
    indices string is nonempty, for each medium component that has custom quantities for at least
    one of the specified indices.

    There is also an option to decide whether the probe should be performed after setup or after
    the full simulation run. The latter option is meaningful if the custom medium state may change
    during the simulation. */
class CustomStateProbe : public SpatialGridWhenFormProbe
{
    ITEM_CONCRETE(CustomStateProbe, SpatialGridWhenFormProbe, "internal spatial grid: custom medium state quantities")
        ATTRIBUTE_TYPE_DISPLAYED_IF(CustomStateProbe, "Level3&CustomMediumState")

        PROPERTY_STRING(indices, "zero-based indices or index ranges of quantities to probe, or empty for all")
        ATTRIBUTE_DEFAULT_VALUE(indices, "")
        ATTRIBUTE_REQUIRED_IF(indices, "false")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
