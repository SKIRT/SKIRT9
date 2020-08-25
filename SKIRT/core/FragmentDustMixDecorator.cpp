/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FragmentDustMixDecorator.hpp"
#include "DustMixFragment.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::setupSelfAfter()
{
    MaterialMix::setupSelfAfter();
    auto scatteringMode = _dustMix->scatteringMode();

    // loop over the grain populations in the dust mix being decorated
    int numPops = _dustMix->numPopulations();
    for (int c = 0; c != numPops; ++c)
    {
        const GrainPopulation* population = _dustMix->population(c);

        // create a dust mix fragment for each grain population
        if (!fragmentSizeBins())
        {
            _fragments.push_back(new DustMixFragment(this, scatteringMode, population));
        }

        // or create a dust mix fragment for each size bin in each grain population
        else
        {
            // construct the size bins (i.e. the bin border points) for this population
            int numSizes = population->numSizes();
            Array borders;
            NR::buildLogGrid(borders, population->sizeDistribution()->amin(), population->sizeDistribution()->amax(),
                             numSizes);

            // loop over the size bins for this population
            for (int s = 0; s != numSizes; ++s)
            {
                _fragments.push_back(new DustMixFragment(this, scatteringMode, population, borders[s], borders[s + 1]));
            }
        }
    }

    _numFrags = _fragments.size();
    if (!_numFrags) throw FATALERROR("There are no fragments in the dust mix to be fragmented");
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

    for (int f = 0; f != _numFrags; ++f)
    {
        string description = "weight for dust mix fragment " + std::to_string(f) + " -- " + populationGrainType(f);
        descriptors.emplace_back(description);
    }
    return descriptors;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> FragmentDustMixDecorator::specificStateVariableInfo() const
{
    vector<StateVariable> descriptors{StateVariable::numberDensity()};

    for (int f = 0; f != _numFrags; ++f)
    {
        string description = "weight for dust mix fragment " + std::to_string(f) + " -- " + populationGrainType(f);
        descriptors.emplace_back(StateVariable::custom(f, description, string()));
    }
    return descriptors;
}

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::initializeSpecificState(MaterialState* state, double /*temperature*/,
                                                       const Array& params) const
{
    for (int f = 0; f != _numFrags; ++f)
    {
        state->setCustom(f, params.size() ? params[f] : 1.);
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
    double opacity = 0.;
    for (int f = 0; f != _numFrags; ++f) opacity += state->custom(f) * _fragments[f]->opacityAbs(lambda, state, pp);
    return opacity;
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    double opacity = 0.;
    for (int f = 0; f != _numFrags; ++f) opacity += state->custom(f) * _fragments[f]->opacitySca(lambda, state, pp);
    return opacity;
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    double opacity = 0.;
    for (int f = 0; f != _numFrags; ++f) opacity += state->custom(f) * _fragments[f]->opacityExt(lambda, state, pp);
    return opacity;
}

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, double w,
                                                 Direction bfkobs, Direction bfky, const MaterialState* state,
                                                 PhotonPacket* pp) const
{
    // calculate the weights corresponding to the scattering opacities for each fragment and their sum
    Array wv(_numFrags);
    double sum = 0.;
    for (int f = 0; f != _numFrags; ++f)
    {
        wv[f] = state->custom(f) * _fragments[f]->opacitySca(lambda, state, pp);
        sum += wv[f];
    }

    // normalize the weights
    if (sum > 0.)
    {
        wv /= sum;

        // perform the peel-off for each fragment with adjusted weights
        for (int f = 0; f != _numFrags; ++f)
            _dustMix->peeloffScattering(I, Q, U, V, lambda, wv[f] * w, bfkobs, bfky, state, pp);
    }
}

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // build the cumulative distribution corresponding to the scattering opacities for each fragment
    Array Xv;
    NR::cdf(Xv, _numFrags, [this, lambda, state, pp](int f) {
        return state->custom(f) * _fragments[f]->opacitySca(lambda, state, pp);
    });

    // randomly select a fragment
    int f = NR::locateClip(Xv, random()->uniform());

    // actually perform the scattering event for this fragment
    _fragments[f]->performScattering(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::indicativeTemperature(const MaterialState* state, const Array& Jv) const
{
    double sumwT = 0.;
    double sumw = 0.;

    for (int f = 0; f != _numFrags; ++f)
    {
        double w = state->custom(f);
        double T = _fragments[f]->indicativeTemperature(state, Jv);
        sumwT += w * T;
        sumw += w;
    }
    return sumw > 0. ? sumwT / sumw : 0.;
}

////////////////////////////////////////////////////////////////////

Array FragmentDustMixDecorator::emissivity(const Array& Jv) const
{
    // TO DO: NEED MATERIAL STATE HERE
    Array ev = _fragments[0]->emissivity(Jv);
    for (int f = 1; f != _numFrags; ++f) ev += _fragments[f]->emissivity(Jv);
    return ev;
}

////////////////////////////////////////////////////////////////////

int FragmentDustMixDecorator::numPopulations() const
{
    return _fragments.size();
}

////////////////////////////////////////////////////////////////////

string FragmentDustMixDecorator::populationGrainType(int f) const
{
    return _fragments[f]->populationGrainType(0);
}

////////////////////////////////////////////////////////////////////

Range FragmentDustMixDecorator::populationSizeRange(int f) const
{
    return _fragments[f]->populationSizeRange(0);
}

////////////////////////////////////////////////////////////////////

const GrainSizeDistribution* FragmentDustMixDecorator::populationSizeDistribution(int f) const
{
    return _fragments[f]->populationSizeDistribution(0);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::populationMass(int f) const
{
    return _fragments[f]->populationMass(0);
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::totalMass() const
{
    return _dustMix->totalMass();
}

////////////////////////////////////////////////////////////////////
