/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COMPOSITEWAVELENGTHGRID_HPP
#define COMPOSITEWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** The CompositeWavelengthGrid class aggregates a number of "child" wavelength grids configured by
    the user into a single, composite wavelength grid. All involved wavelength grids (i.e. all
    child grids and the composited grid) are disjoint wavelength grids, i.e. they inherit from the
    DisjointWavelengthGrid class. These grids have non-overlapping but possibly adjacent wavelength
    bins with constant maximum transmission within the bins and zero transmission outside of the
    bins. Refer to the DisjointWavelengthGrid class description for a more formal definition. Many
    of the wavelength grids frequently used in SKIRT have these properties. Refer to the list of
    DisjointWavelengthGrid subclasses for more information.

    The \em wavelengthGrids property specifies a nonempty list of disjoint wavelength grids of any
    type and with arbitrary configuration options. To construct a single combined grid, each of the
    child grids is processed one by one in order of occurrence in the user configuration. After
    initializing the result buffer to an empty grid, the procedure composites each child grid in
    turn into the current contents of the result buffer. Specifically, the procedure replaces any
    bins in the result buffer that are overlapped by child bins by those child bins, splitting
    partially overlapped bins where needed. The characteristic wavelengths of the bins are
    preserved where possible. If needed, a new characteristic wavelength is calculated for a
    shortened (partially overlapped) bin from the new bin borders using linear or logarithmic
    interpolation as indicated by the \em log configuration option.

    It is worth noting the following:

    - If there is just one child, the composite wavelength grid is merely a copy of that child.

    - As long as all child wavelength ranges are mutually disjoint, the order in which the children
    are specified is irrelevant because the wavelength bins will be automatically sorted by the
    compositing process. As soon as there is overlap, however, the ordering of the children
    determines which child’s wavelength bins will be preserved in the composite wavelength grid
    for the overlapping ranges. It is therefore usually preferable to sort the child wavelength
    grids from lower to higher resolution.

    - It is possible to nest a composite wavelength grid inside another composite wavelength grid.
    This allows joining multiple groups of child wavelength grids each with a different setting of
    the \em log option. However, this seems useful only in pathetic cases.

    */
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
