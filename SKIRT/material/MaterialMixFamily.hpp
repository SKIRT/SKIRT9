/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALMIXFAMILY_HPP
#define MATERIALMIXFAMILY_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"
#include "SnapshotParameter.hpp"
class MaterialMix;

//////////////////////////////////////////////////////////////////////

/** MaterialMixFamily is an abstract class for representing a family of material mixes of the same
    material type (i.e. either dust, electrons or gas; never a mixture of these). Each subclass
    implements a material mix family from which a particular material mix can be chosen based on
    one or more parameters. This base class offers a generic interface for obtaining information on
    the number and type of parameters appropriate for the subclass, and for actually selecting a
    material mix given specific parameter values.

    Because each material mix is represented by a distinct C++ object, a material mix family must
    by necessity consist of a discrete number of family members. As a result, each subclass must
    implement a mechanism to map specified parameter values from their continuous ranges to one of
    a discrete set of material mixes.

    Also, because a material mix object often consumes a significant amount of memory for holding
    material properties and precalculated values, the number of material mixes that can be offered
    by a material mix family is limited by the available memory resources. */
class MaterialMixFamily : public SimulationItem
{
    ITEM_ABSTRACT(MaterialMixFamily, SimulationItem, "a family of material mixes")
    ITEM_END()

public:
    /** This function returns the number and type of parameters used by this particular material
        mix family as a list of SnapshotParameter objects. Each of these objects specifies unit
        information and a human-readable descripton for the parameter. */
    virtual vector<SnapshotParameter> parameterInfo() const = 0;

    /** This function returns (a pointer to) the material mix corresponding to the given parameter
        values, or some default material mix if one or more of the parameter values are out of
        range. The material mix family retains ownership of the returned material mix, and
        guarantees that it will not be destroyed until the family itself is destroyed.

        The number and type of parameters must match the information returned by the
        parameterInfo() function; if not the behavior is undefined. */
    virtual const MaterialMix* mix(const Array& parameters) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
