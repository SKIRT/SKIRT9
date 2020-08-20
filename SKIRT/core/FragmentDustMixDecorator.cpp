/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FragmentDustMixDecorator.hpp"
#include "MaterialState.hpp"

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::setupSelfAfter()
{
    MaterialMix::setupSelfAfter();
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType FragmentDustMixDecorator::materialType() const
{
    return MaterialType::Dust;
}

////////////////////////////////////////////////////////////////////

bool FragmentDustMixDecorator::hasPolarizedScattering() const
{
    return _dustMix->hasPolarizedScattering();
}

////////////////////////////////////////////////////////////////////

bool FragmentDustMixDecorator::hasStochasticDustEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FragmentDustMixDecorator::parameterInfo() const
{
    vector<SnapshotParameter> descriptors;

    int numPops = _dustMix->numPopulations();
    for (int c = 0; c != numPops; ++c)
    {
        string description = "number density fraction for grain population " + std::to_string(c) + " -- "
                             + _dustMix->populationGrainType(c);
        descriptors.emplace_back(description);
    }
    return descriptors;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> FragmentDustMixDecorator::specificStateVariableInfo() const
{
    vector<StateVariable> descriptors{StateVariable::numberDensity()};

    int numPops = _dustMix->numPopulations();
    for (int c = 0; c != numPops; ++c)
    {
        string description = "number density fraction for grain population " + std::to_string(c) + " -- "
                             + _dustMix->populationGrainType(c);
        descriptors.emplace_back(StateVariable::custom(c, description, string()));
    }
    return descriptors;
}

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::initializeSpecificState(MaterialState* state, double /*temperature*/,
                                                       const Array& params) const
{
    if (params.size())
    {
        double norm = params.sum();
        int numPops = _dustMix->numPopulations();
        for (int c = 0; c != numPops; ++c)
        {
            state->setCustom(c, params[c] / norm);
        }
    }
    else
    {
        // TO DO.
    }
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::mass() const
{
    return _dustMix->mass();
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::sectionAbs(double lambda) const
{
    return _dustMix->sectionAbs(lambda);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::sectionSca(double lambda) const
{
    return _dustMix->sectionSca(lambda);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::sectionExt(double lambda) const
{
    return _dustMix->sectionExt(lambda);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::asymmpar(double lambda) const
{
    return _dustMix->asymmpar(lambda);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return _dustMix->opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return _dustMix->opacitySca(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return _dustMix->opacityExt(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, double w,
                                                 Direction bfkobs, Direction bfky, const MaterialState* state,
                                                 PhotonPacket* pp) const
{
    return _dustMix->peeloffScattering(I, Q, U, V, lambda, w, bfkobs, bfky, state, pp);
}

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    return _dustMix->performScattering(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::indicativeTemperature(const MaterialState* state, const Array& Jv) const
{
    return _dustMix->indicativeTemperature(state, Jv);
}

////////////////////////////////////////////////////////////////////

Array FragmentDustMixDecorator::emissivity(const Array& Jv) const
{
    return _dustMix->emissivity(Jv);
}

////////////////////////////////////////////////////////////////////

int FragmentDustMixDecorator::numPopulations() const
{
    return _dustMix->numPopulations();
}

////////////////////////////////////////////////////////////////////

string FragmentDustMixDecorator::populationGrainType(int c) const
{
    return _dustMix->populationGrainType(c);
}

////////////////////////////////////////////////////////////////////

Range FragmentDustMixDecorator::populationSizeRange(int c) const
{
    return _dustMix->populationSizeRange(c);
}

////////////////////////////////////////////////////////////////////

const GrainSizeDistribution* FragmentDustMixDecorator::populationSizeDistribution(int c) const
{
    return _dustMix->populationSizeDistribution(c);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::populationMass(int c) const
{
    return _dustMix->populationMass(c);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::totalMass() const
{
    return _dustMix->totalMass();
}

////////////////////////////////////////////////////////////////////
