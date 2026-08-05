/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOXCELLDENSITYMIXIN_HPP
#define BOXCELLDENSITYMIXIN_HPP

#include "Box.hpp"
#include "DensityInCellInterface.hpp"
class MassInBoxInterface;
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** The BoxCellDensityMixIn class implements the DensityInCellInterface for spatial grids that have
    spatial cells in the form of axes-aligned bounding boxes, and that reside in a simulation where
    all medium components offer the MassInBoxInterface (e.g., imported smoothed particles).

    Eligible spatial grid classes should inherit ("mix-in") this class so that the medium density
    in a cell can be calculated directly without the need for random sampling. */
class BoxCellDensityMixIn : public DensityInCellInterface
{
    //============= Construction - Setup - Destruction =============

protected:
    /** The default constructor initializes the data members to indicate that the
        DensityInCellInterface is disabled. */
    BoxCellDensityMixIn();

    /** This function must be called from the inheriting class during setup. It verifies whether
        the DensityInCellInterface can actually be enabled and if so, caches some information to
        accelerate the operation later on. Specifically, the interface can be enabled only if all
        medium components in the simulation offer the MassInBoxInterface.

        The argument is used to locate the medium components in the simulation; it should point to
        the inheriting class. */
    void setup(SimulationItem* item);

    /** This function must be called from the inheriting class to implement its offersInterface()
        function. It returns true if the DensityInCellInterface is enabled, and false if not. The
        function assumes that setup() has properly been called. */
    bool offersInterface() const;

    //======================== Other Functions =======================

public:
    /** This function implements the DensityInCellInterface. It returns the number density for
        medium component \em h in the spatial grid cell with index \em m. */
    double numberDensity(int h, int m) const override;

protected:
    /** This function must be implemented in the inheriting class. It returns the axis-aligned
        bounding box defining the cell with index \f$m\f$, */
    virtual Box cellBox(int m) const = 0;

    //======================== Data Members ========================

private:
    bool _enabled{false};               // true if all medium component offer MassInBoxInterface
    vector<MassInBoxInterface*> _mibv;  // MassInBoxInterface ptr for each medium component
};

////////////////////////////////////////////////////////////////////

#endif
