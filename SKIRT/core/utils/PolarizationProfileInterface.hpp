/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POLARIZATIONPROFILEINTERFACE_HPP
#define POLARIZATIONPROFILEINTERFACE_HPP

#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

/** PolarizationProfileInterface is a pure interface to obtain the polarization state of the
    radiation emitted by a source into a given direction \f$(\theta,\phi)\f$. */
class PolarizationProfileInterface
{
protected:
    /** The empty constructor for the interface. */
    PolarizationProfileInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~PolarizationProfileInterface() {}

    /** This function returns the Stokes vector defining the polarization state of the radiation
        emitted into the given direction \f$(\theta,\phi)\f$. For unpolarized emission, this
        function would return a default-constructed StokesVector instance. */
    virtual StokesVector polarizationForDirection(Direction bfk) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
