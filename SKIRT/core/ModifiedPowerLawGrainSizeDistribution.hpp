/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MODIFIEDPOWERLAWGRAINSIZEDISTRIBUTION_HPP
#define MODIFIEDPOWERLAWGRAINSIZEDISTRIBUTION_HPP

#include "RangeGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** ModifiedPowerLawGrainSizeDistribution is a GrainSizeDistribution subclass that represents a
    dust grain size distribution of the form \f[ \frac{\text{d}n_\text{D}}{\text{d}a} \propto
    a^{\alpha} \,f_\text{ed}(a) \,f_\text{cv}(a) \qquad\text{for}\quad a_\text{min} \leq a \leq
    a_\text{max}, \f] with an exponential decay term \f[ f_\text{ed}(a) = \begin{cases} 1 & \quad
    a\leq a_\text{t} \\ \exp\left(-[(a-a_\text{t})/a_\text{c}]^\gamma \right) & \quad a>a_\text{t}
    \end{cases} \f] and a curvature term \f[ f_\text{cv}(a) = \left[ 1+|\zeta|\,(a/a_u)^\eta
    \right]^{\text{sign}(\zeta)}. \f]

    The size range of the distribution can be configured in the RangeGrainSizeDistribution base
    class. The remaining seven parameters \f$\alpha\f$, \f$a_\text{t}\f$, \f$a_\text{c}\f$,
    \f$\gamma\f$, \f$a_\text{u}\f$, \f$\zeta\f$ and \f$\eta\f$ can be configured as attributes in
    this class. The function is scaled arbitrarily.

    The functional form for the grain size distribution implemented by this class is inspired by
    the DustEM code, which is described in Compiègne et al. 2011 (AA, 525, A103) and can be
    downloaded from http://www.ias.u-psud.fr/DUSTEM/. */
class ModifiedPowerLawGrainSizeDistribution : public RangeGrainSizeDistribution
{
    ITEM_CONCRETE(ModifiedPowerLawGrainSizeDistribution, RangeGrainSizeDistribution,
                  "a modified power-law dust grain size distribution")

        PROPERTY_DOUBLE(powerLawIndex, "the index α of the power law")
        ATTRIBUTE_MIN_VALUE(powerLawIndex, "[-10")
        ATTRIBUTE_MAX_VALUE(powerLawIndex, "0[")
        ATTRIBUTE_DEFAULT_VALUE(powerLawIndex, "-3.5")

        PROPERTY_DOUBLE(turnOffPoint, "the turn-off point a_t in the exponential decay term")
        ATTRIBUTE_QUANTITY(turnOffPoint, "grainsize")
        ATTRIBUTE_MIN_VALUE(turnOffPoint, "[0")
        ATTRIBUTE_MAX_VALUE(turnOffPoint, "1 mm]")
        ATTRIBUTE_DEFAULT_VALUE(turnOffPoint, "0.1 micron")

        PROPERTY_DOUBLE(scaleExponentialDecay, "the scale a_c in the exponential decay term")
        ATTRIBUTE_QUANTITY(scaleExponentialDecay, "grainsize")
        ATTRIBUTE_MIN_VALUE(scaleExponentialDecay, "[1 Angstrom")
        ATTRIBUTE_MAX_VALUE(scaleExponentialDecay, "1 mm]")
        ATTRIBUTE_DEFAULT_VALUE(scaleExponentialDecay, "0.1 micron")

        PROPERTY_DOUBLE(exponentExponentialDecay, "the exponent γ in the exponential decay term")
        ATTRIBUTE_MIN_VALUE(exponentExponentialDecay, "[0")
        ATTRIBUTE_MAX_VALUE(exponentExponentialDecay, "10]")
        ATTRIBUTE_DEFAULT_VALUE(exponentExponentialDecay, "3")

        PROPERTY_DOUBLE(scaleCurvature, "the scale a_u in the curvature term")
        ATTRIBUTE_QUANTITY(scaleCurvature, "grainsize")
        ATTRIBUTE_MIN_VALUE(scaleCurvature, "[1 Angstrom")
        ATTRIBUTE_MAX_VALUE(scaleCurvature, "1 mm]")
        ATTRIBUTE_DEFAULT_VALUE(scaleCurvature, "0.1 micron")

        PROPERTY_DOUBLE(strengthCurvature, "the strength ζ in the curvature term")
        ATTRIBUTE_MIN_VALUE(strengthCurvature, "[-10")
        ATTRIBUTE_MAX_VALUE(strengthCurvature, "10]")
        ATTRIBUTE_DEFAULT_VALUE(strengthCurvature, "0.3")

        PROPERTY_DOUBLE(exponentCurvature, "the exponent η in the curvature term")
        ATTRIBUTE_MIN_VALUE(exponentCurvature, "[0")
        ATTRIBUTE_MAX_VALUE(exponentCurvature, "10]")
        ATTRIBUTE_DEFAULT_VALUE(exponentCurvature, "1")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        grain size distribution object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), its
        properties have been initialized to the specified values, and its setup() function has been
        called. */
    explicit ModifiedPowerLawGrainSizeDistribution(SimulationItem* parent, double minSize, double maxSize,
                                                   double powerLawIndex, double turnOffPoint,
                                                   double scaleExponentialDecay, double exponentExponentialDecay,
                                                   double scaleCurvature, double strengthCurvature,
                                                   double exponentCurvature);

    //======================== Other Functions =======================

public:
    /** This function returns the value of \f$\frac{\text{d}n_\text{D}}{\text{d}a}\f$ as described
        in the header for this class (with an arbitrary proportionality factor of one). */
    double dnda(double a) const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation
    const double& _alpha{_powerLawIndex};
    const double& _at{_turnOffPoint};
    const double& _ac{_scaleExponentialDecay};
    const double& _gamma{_exponentExponentialDecay};
    const double& _au{_scaleCurvature};
    const double& _zeta{_strengthCurvature};
    const double& _eta{_exponentCurvature};
};

////////////////////////////////////////////////////////////////////

#endif
