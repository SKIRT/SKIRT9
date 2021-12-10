/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COMPOSITEWAVELENGTHGRID_HPP
#define COMPOSITEWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** CompositeWavelengthGrid is a subclass of the DisjointWavelengthGrid class representing
    wavelength grids composited from a list of disjoint wavelength grids. */
class CompositeWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(CompositeWavelengthGrid, DisjointWavelengthGrid,
                  "a wavelength grid composited from a list of wavelength grids")
        ATTRIBUTE_TYPE_DISPLAYED_IF(CompositeWavelengthGrid, "Level2")

        PROPERTY_ITEM_LIST(wavelengthGrids, DisjointWavelengthGrid, "the wavelength grids to be composited")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthGrids, "LogWavelengthGrid")

        PROPERTY_BOOL(log, "use logarithmic scale")
        ATTRIBUTE_DEFAULT_VALUE(log, "true")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs the composite wavelength grid from the configured list of
        wavelength grids. */
    void setupSelfAfter() override;
};

////////////////////////////////////////////////////////////////////

#endif
