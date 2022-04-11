/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIATIONFIELDPROBE_HPP
#define RADIATIONFIELDPROBE_HPP

#include "SpatialGridFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** RadiationFieldProbe probes the mean radiation field intensity as calculated by the simulation
    and discretized on the spatial grid of the simulation. The probe outputs a compound quantity,
    \f$J_\ell\f$, with \f$\ell=1 \ldots N_\mathrm{rf}\f$, i.e. a value for each of the wavelength
    bins in the radiation field wavelength grid configured for the simulation. The spectral flavor
    and ordering of the components and the output units conform to the \em wavelengthOutputStyle,
    \em fluxOutputStyle, and unit system configured for the simulation.

    The probe can be used with any Form subclass. When associated with a form that projects the
    quantity along a path, the value is density-weighted.

    The probe offers an option to output a separate text column file with details on the radiation
    field wavelength grid. For each wavelength bin, the file lists the characteristic wavelength,
    the wavelength bin width, and the left and right borders of the bin. */
class RadiationFieldProbe : public SpatialGridFormProbe
{
    ITEM_CONCRETE(RadiationFieldProbe, SpatialGridFormProbe, "the mean radiation field intensity")
        ATTRIBUTE_TYPE_DISPLAYED_IF(RadiationFieldProbe, "Level2&SpatialGrid&RadiationField")

        PROPERTY_BOOL(writeWavelengthGrid, "output a text file with the radiation field wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(writeWavelengthGrid, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns an enumeration indicating that probing must be performed at the end
        of the simulation. */
    When when() const override;

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
