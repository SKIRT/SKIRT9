/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MESH_HPP
#define MESH_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"

//////////////////////////////////////////////////////////////////////

/** The Mesh class is an abstract base class that characterizes different types of one-dimensional
    meshes over the unit interval \f$[0,1]\f$. A mesh is essentially a partition of this interval
    into a number of \f$N\f$ finite bins. Internally, a mesh consists of an ordered array of
    \f$N+1\f$ mesh points \f$\{t_i\}\f$, with \f$t_0=0\f$ and \f$t_N=1\f$. The different subclasses
    of the Mesh class indicate different mesh point distributions, such as linear distributions,
    etc. */
class Mesh : public SimulationItem
{
    ITEM_ABSTRACT(Mesh, SimulationItem, "a mesh")

        PROPERTY_INT(numBins, "the number of bins in the mesh")
        ATTRIBUTE_MIN_VALUE(numBins, "1")
        ATTRIBUTE_MAX_VALUE(numBins, "100000")
        ATTRIBUTE_DEFAULT_VALUE(numBins, "100")

    ITEM_END()

    //======== Setters for Discoverable Attributes =======

protected:
    /** Sets the number of bins in the mesh. Can be used by subclasses that wish to override the
        number of bins configured by the corresponding property in this base class. */
    void setNumBins(int value);

    //======================== Other Functions =======================

public:
    /** This pure virtual function returns an array containing the \f$N+1\f$ mesh points
        \f$\{t_i\}\f$ in ascending order, with \f$t_0=0\f$ and \f$t_N=1\f$. */
    virtual Array mesh() const = 0;
};

//////////////////////////////////////////////////////////////////////

#endif
