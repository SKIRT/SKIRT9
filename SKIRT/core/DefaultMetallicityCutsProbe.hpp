/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTMETALLICITYCUTSPROBE_HPP
#define DEFAULTMETALLICITYCUTSPROBE_HPP

#include "StateProbe.hpp"

////////////////////////////////////////////////////////////////////

/** For each medium in the simulation that has a metallicity state variable (as requested by the
    associated material mix), DefaultMetallicityCutsProbe outputs FITS files with cuts through the
    metallicity along the coordinate planes. The field of view of each cut covers the extent of the
    spatial grid in the simulation in the relevant directions. Each cut has 1024 x 1024 pixels.

    The files are named <tt>prefix_probe_Z_XX_N.fits</tt> where XX indicates the orientation of the
    cut and N is replaced with the zero-based index of the medium in the configuration (i.e. in the
    ski file). Each file contains a single image frame.

    Note that dust components do not store the imported metallicity; for those components,
    metallicity is simply used as a multiplier to calculate the mass density. */
class DefaultMetallicityCutsProbe : public StateProbe
{
    ITEM_CONCRETE(DefaultMetallicityCutsProbe, StateProbe, "cuts of the metallicity along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultMetallicityCutsProbe, "Level2&Gas&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs the probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
