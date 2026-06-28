/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INTEGRATEDSECONDARYLINELUMINOSITYPROBE_HPP
#define INTEGRATEDSECONDARYLINELUMINOSITYPROBE_HPP

#include "SpecialtyWhenProbe.hpp"

////////////////////////////////////////////////////////////////////

/** IntegratedSecondaryLineLuminosityProbe outputs a text column file with the spatially
    integrated luminosity in each emission line for every gas medium component that supports line
    emission. For each line, the file lists the central wavelength and the total luminosity summed
    over all spatial cells.

    The probe is intended as a global diagnostic counterpart to SecondaryLineLuminosityProbe,
    which writes per-cell luminosities through a spatial form. Use this probe when only the
    integrated per-line totals are needed (e.g. for line-flux comparisons against 1D models). */
class IntegratedSecondaryLineLuminosityProbe : public SpecialtyWhenProbe
{
    ITEM_CONCRETE(IntegratedSecondaryLineLuminosityProbe, SpecialtyWhenProbe,
                  "secondary: spatially integrated luminosities per emission line")
        ATTRIBUTE_TYPE_DISPLAYED_IF(IntegratedSecondaryLineLuminosityProbe, "GasEmission")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
