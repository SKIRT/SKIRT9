/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIGRAINPOPULATIONINTERFACE_HPP
#define MULTIGRAINPOPULATIONINTERFACE_HPP

#include "Range.hpp"
class GrainSizeDistribution;

////////////////////////////////////////////////////////////////////

/** MultiGrainPopulationInterface is a pure interface exposing information about the individual
    grain populations making up a multi-grain dust mixture. It is implemented by the
    MultiGrainDustMix and CompositeDustMix classes so that probes can retrieve such details. */
class MultiGrainPopulationInterface
{
protected:
    /** The empty constructor for the interface. */
    MultiGrainPopulationInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~MultiGrainPopulationInterface() {}

    /** This function returns the number of dust grain populations (with indices \f$c\f$). Each
        grain population represents the combination of a grain composition, providing the optical
        and calorimetric properties of the grain material, and a grain size distribution with some
        normalization to specify the the amount of dust contained in the population. No grain size
        discretization has been applied to these populations. */
    virtual int numPopulations() const = 0;

    /** This function returns a brief human-readable identifier for the type of grain material
        represented by the population with index \f$c\f$. The identifier does not contain white
        space. */
    virtual string populationGrainType(int c) const = 0;

    /** This function returns the bulk mass density \f$\rho_\text{bulk}\f$ of the grain material
        represented by the population with index \f$c\f$. */
    virtual double populationBulkDensity(int c) const = 0;

    /** This function returns the minimum and maximum grain sizes \f$a_{\text{min},c},
        a_{\text{max},c}\f$ for the population with index \f$c\f$. */
    virtual Range populationSizeRange(int c) const = 0;

    /** This function returns the grain size distribution object for the population with index
        \f$c\f$. */
    virtual const GrainSizeDistribution* populationSizeDistribution(int c) const = 0;

    /** This function returns the dust mass \f$\mu_c\f$ per hydrogen atom for the population with
        index \f$c\f$. */
    virtual double populationMass(int c) const = 0;

    /** This function returns the total dust mass \f$\mu\f$ per hydrogen atom for all populations
        combined. */
    virtual double totalMass() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
