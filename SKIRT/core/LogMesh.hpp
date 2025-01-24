/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LOGMESH_HPP
#define LOGMESH_HPP

#include "Mesh.hpp"

//////////////////////////////////////////////////////////////////////

/** The LogMesh class represents an origin-anchored, logarithmic mesh. The first bin covers the
    interval \f$[0,t_\text{c}]\f$ and the widths of the remaining bins are distributed
    logarithmically over the interval \f$[t_\text{c},1]\f$, where \f$t_\text{c}\f$ is a
    user-configurable parameter. If the mesh has only one bin, the value of \f$t_\text{c}\f$ is
    ignored and the single bin spans the complete interval \f$[0,1]\f$. */
class LogMesh : public Mesh
{
    ITEM_CONCRETE(LogMesh, Mesh, "a logarithmic mesh")

        PROPERTY_DOUBLE(centralBinFraction, "the central bin width fraction")
        ATTRIBUTE_MIN_VALUE(centralBinFraction, "]0")
        ATTRIBUTE_MAX_VALUE(centralBinFraction, "1[")
        ATTRIBUTE_DEFAULT_VALUE(centralBinFraction, "1e-3")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns an array containing the mesh points. */
    Array mesh() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
