/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALSTATE_HPP
#define MATERIALSTATE_HPP

#include "MediumSystem.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the MaterialState class provides access to the medium state for a particular
    spatial cell and medium component, as determined during construction. This includes the common
    state variables for the cell as well as the specific state variables for the medium component.

    This is a temporary implementation of this class. */
class MaterialState
{
    //============= Construction - Setup - Destruction =============

public:
    /** The constructor arguments specify the spatial cell and medium component to be represented.
        */
    MaterialState(const MediumSystem* ms, int m, int h)
    {
        _ms = ms;
        _m = m;
        _h = h;
    }

    //======================== Other Functions =======================

public:
    /** This function returns the volume \f$V\f$ of the spatial cell. */
    double volume() const { return _ms->volume(_m); }

    /** This function returns the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in
        the spatial cell. */
    Vec bulkVelocity() const { return _ms->bulkVelocity(_m); }

    /** This function returns the magnetic field \f${\boldsymbol{B}}\f$ in the spatial cell. */
    Vec magneticField() const { return _ms->magneticField(_m); }

    /** This function returns the gas temperature \f$T\f$ of the medium component in the spatial
        cell. */
    double temperature() const { return _ms->temperature(_m, _h); }

    /** This function returns the number density of the medium component in the spatial cell. */
    double numberDensity() const { return _ms->numberDensity(_m, _h); }

    //======================== Data Members ========================

private:
    const MediumSystem* _ms{nullptr};
    int _m{-1};
    int _h{-1};
};

////////////////////////////////////////////////////////////////

#endif
