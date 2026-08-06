/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AbsorptionOnlyMaterialMixDecorator.hpp"
#include "FatalError.hpp"

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

bool AbsorptionOnlyMaterialMixDecorator::offersInterface(const std::type_info& interfaceTypeInfo) const
{
    if (interfaceTypeInfo == typeid(MultiGrainPopulationInterface))
    {
        return materialMix()->interface<MultiGrainPopulationInterface>(0, false) != nullptr;
    }
    return MaterialMix::offersInterface(interfaceTypeInfo);
}

////////////////////////////////////////////////////////////////////

int AbsorptionOnlyMaterialMixDecorator::numPopulations() const
{
    const auto* mgpi = materialMix()->interface<MultiGrainPopulationInterface>(0, false);
    if (mgpi) return mgpi->numPopulations();
    throw FATALERROR("This function should only be called for a multi-grain dust mix");
}

////////////////////////////////////////////////////////////////////

string AbsorptionOnlyMaterialMixDecorator::populationGrainType(int c) const
{
    const auto* mgpi = materialMix()->interface<MultiGrainPopulationInterface>(0, false);
    if (mgpi) return mgpi->populationGrainType(c);
    throw FATALERROR("This function should only be called for a multi-grain dust mix");
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::populationBulkDensity(int c) const
{
    const auto* mgpi = materialMix()->interface<MultiGrainPopulationInterface>(0, false);
    if (mgpi) return mgpi->populationBulkDensity(c);
    throw FATALERROR("This function should only be called for a multi-grain dust mix");
}

////////////////////////////////////////////////////////////////////

Range AbsorptionOnlyMaterialMixDecorator::populationSizeRange(int c) const
{
    const auto* mgpi = materialMix()->interface<MultiGrainPopulationInterface>(0, false);
    if (mgpi) return mgpi->populationSizeRange(c);
    throw FATALERROR("This function should only be called for a multi-grain dust mix");
}

////////////////////////////////////////////////////////////////////

const GrainSizeDistribution* AbsorptionOnlyMaterialMixDecorator::populationSizeDistribution(int c) const
{
    const auto* mgpi = materialMix()->interface<MultiGrainPopulationInterface>(0, false);
    if (mgpi) return mgpi->populationSizeDistribution(c);
    throw FATALERROR("This function should only be called for a multi-grain dust mix");
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::populationMass(int c) const
{
    const auto* mgpi = materialMix()->interface<MultiGrainPopulationInterface>(0, false);
    if (mgpi) return mgpi->populationMass(c);
    throw FATALERROR("This function should only be called for a multi-grain dust mix");
}

////////////////////////////////////////////////////////////////////

double AbsorptionOnlyMaterialMixDecorator::totalMass() const
{
    const auto* mgpi = materialMix()->interface<MultiGrainPopulationInterface>(0, false);
    if (mgpi) return mgpi->totalMass();
    throw FATALERROR("This function should only be called for a multi-grain dust mix");
}

////////////////////////////////////////////////////////////////////
