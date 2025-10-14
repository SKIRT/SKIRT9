/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AbsorptionOnlyMaterialMixDecorator.hpp"

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType AbsorptionOnlyMaterialMixDecorator::materialType() const
{
    return materialMix()->materialType();
}

////////////////////////////////////////////////////////////////////

bool AbsorptionOnlyMaterialMixDecorator::hasPolarizedScattering() const
{
    return materialMix()->hasPolarizedScattering();
}

////////////////////////////////////////////////////////////////////

bool AbsorptionOnlyMaterialMixDecorator::hasPolarizedAbsorption() const
{
    return materialMix()->hasPolarizedAbsorption();
}

////////////////////////////////////////////////////////////////////

bool AbsorptionOnlyMaterialMixDecorator::hasResonantScattering() const
{
    return materialMix()->hasResonantScattering();
}

////////////////////////////////////////////////////////////////////

bool AbsorptionOnlyMaterialMixDecorator::hasNegativeExtinction() const
{
    return materialMix()->hasNegativeExtinction();
}

////////////////////////////////////////////////////////////////////

bool AbsorptionOnlyMaterialMixDecorator::hasExtraSpecificState() const
{
    return materialMix()->hasExtraSpecificState();
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType AbsorptionOnlyMaterialMixDecorator::hasDynamicMediumState() const
{
    return materialMix()->hasDynamicMediumState();
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> AbsorptionOnlyMaterialMixDecorator::parameterInfo() const
{
    return materialMix()->parameterInfo();
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> AbsorptionOnlyMaterialMixDecorator::specificStateVariableInfo() const
{
    return materialMix()->specificStateVariableInfo();
}

////////////////////////////////////////////////////////////////////

void AbsorptionOnlyMaterialMixDecorator::initializeSpecificState(MaterialState* state, double metallicity,
                                                                 double temperature, const Array& params) const
{
    return materialMix()->initializeSpecificState(state, metallicity, temperature, params);
}

////////////////////////////////////////////////////////////////////

UpdateStatus AbsorptionOnlyMaterialMixDecorator::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    return materialMix()->updateSpecificState(state, Jv);
}

////////////////////////////////////////////////////////////////////

bool AbsorptionOnlyMaterialMixDecorator::isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged,
                                                                  MaterialState* currentAggregate,
                                                                  MaterialState* previousAggregate) const
{
    return materialMix()->isSpecificStateConverged(numCells, numUpdated, numNotConverged, currentAggregate,
                                                   previousAggregate);
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::mass() const
{
    return materialMix()->mass();
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::sectionAbs(double lambda) const
{
    return materialMix()->sectionAbs(lambda);
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::sectionExt(double lambda) const
{
    return materialMix()->sectionAbs(lambda);
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::opacityAbs(double lambda, const MaterialState* state,
                                                      const PhotonPacket* pp) const
{
    return materialMix()->opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::opacitySca(double /*lambda*/, const MaterialState* /*state*/,
                                                      const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::opacityExt(double lambda, const MaterialState* state,
                                                      const PhotonPacket* pp) const
{
    return materialMix()->opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::indicativeTemperature(const MaterialState* state, const Array& Jv) const
{
    return materialMix()->indicativeTemperature(state, Jv);
}

////////////////////////////////////////////////////////////////////
