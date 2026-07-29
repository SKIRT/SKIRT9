/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MASSINBOXINTERFACE_HPP
#define MASSINBOXINTERFACE_HPP

#include "Basics.hpp"
class Box;

////////////////////////////////////////////////////////////////////

/** MassInBoxInterface is a pure interface implemented by medium components that are capable of
    efficiently and accurately calculating the medium mass in a given axis-aligned bounding box
    (without random sampling). */
class MassInBoxInterface
{
protected:
    /** The empty constructor for the interface. */
    MassInBoxInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~MassInBoxInterface() {}

    /** This function returns the mass for the medium component in the specified axis-aligned
        bounding box. */
    virtual double massInBox(const Box& box) const = 0;

    /** This function returns the volume-integrated number density for the medium component in the
        specified axis-aligned bounding box. */
    virtual double numberInBox(const Box& box) const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
