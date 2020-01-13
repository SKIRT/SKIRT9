/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POWMESH_HPP
#define POWMESH_HPP

#include "MoveableMesh.hpp"

//////////////////////////////////////////////////////////////////////

/** The PowMesh class represents meshes on the unit interval \f$[0,1]\f$ with a power-law
    distribution of the mesh points. This distribution is such that the bin sizes form a geometric
    series, i.e. each bin is a constant factor larger than the previous one. The actual
    distribution is characterized by a single parameter: the bin width ratio between the last and
    the first bin, \f[ {\cal{R}} = \frac{t_N-t_{N-1}}{t_1-t_0}.\f] This ratio can be larger than
    one (in which case the first bin is the smallest) or smaller than one (in which case the last
    bin is the smallest). If the mesh has only one bin, this single bin spans the complete interval
    \f$[0,1]\f$. */
class PowMesh : public MoveableMesh
{
    ITEM_CONCRETE(PowMesh, MoveableMesh, "a power-law mesh")

        PROPERTY_DOUBLE(ratio, "the bin width ratio between the last and the first bin")
        ATTRIBUTE_MIN_VALUE(ratio, "]0")
        ATTRIBUTE_DEFAULT_VALUE(ratio, "1")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns an array containing the mesh points. */
    Array mesh() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
