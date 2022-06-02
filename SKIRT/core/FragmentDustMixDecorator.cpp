/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FragmentDustMixDecorator.hpp"
#include "Configuration.hpp"
#include "DustMixFragment.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"
#include "Random.hpp"
#include "ShortArray.hpp"

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
                _fragments.push_back(new DustMixFragment(this, scatteringMode, population, borders[s], borders[s + 1],
                                                         _dustMix->populationNormalization(c)));
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

bool FragmentDustMixDecorator::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool FragmentDustMixDecorator::hasContinuumEmission() const
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
        descriptors.push_back(SnapshotParameter::custom(description));
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

    if (hasDynamicDensities())
    {
        for (int f = 0; f != _numFrags; ++f)
        {
            string description =
                "dynamic density fraction for dust mix fragment " + std::to_string(f) + " -- " + populationGrainType(f);
            descriptors.emplace_back(StateVariable::custom(_numFrags + f, description, string()));
        }
    }
    return descriptors;
}

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::initializeSpecificState(MaterialState* state, double /*metallicity*/,
                                                       double /*temperature*/, const Array& params) const
{
    for (int f = 0; f != _numFrags; ++f)
    {
        state->setCustom(f, params.size() ? params[f] : 1.);
    }

    if (hasDynamicDensities())
    {
        for (int f = 0; f != _numFrags; ++f)
        {
            state->setCustom(_numFrags + f, initialDensityFraction());
        }
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

// define macro for fragment weight combining initial weight and dynamic density fraction, if available
#define WEIGHT(f) (_hasDynamicDensities ? state->custom(f) * state->custom(_numFrags + f) : state->custom(f))

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    double opacity = 0.;
    for (int f = 0; f != _numFrags; ++f) opacity += WEIGHT(f) * _fragments[f]->opacityAbs(lambda, state, pp);
    return opacity;
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    double opacity = 0.;
    for (int f = 0; f != _numFrags; ++f) opacity += WEIGHT(f) * _fragments[f]->opacitySca(lambda, state, pp);
    return opacity;
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    double opacity = 0.;
    for (int f = 0; f != _numFrags; ++f) opacity += WEIGHT(f) * _fragments[f]->opacityExt(lambda, state, pp);
    return opacity;
}

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda,
                                                 Direction bfkobs, Direction bfky, const MaterialState* state,
                                                 const PhotonPacket* pp) const
{
    // calculate the weights corresponding to the scattering opacities for each fragment and their sum
    ShortArray wv(_numFrags);
    double sum = 0.;
    for (int f = 0; f != _numFrags; ++f)
    {
        wv[f] = WEIGHT(f) * _fragments[f]->opacitySca(lambda, state, pp);
        sum += wv[f];
    }

    // perform the peel-off for each fragment with corresponding normalized weights
    if (sum > 0.)
    {
        for (int f = 0; f != _numFrags; ++f)
        {
            double If = 0., Qf = 0., Uf = 0., Vf = 0.;
            _dustMix->peeloffScattering(If, Qf, Uf, Vf, lambda, bfkobs, bfky, state, pp);
            double wf = wv[f] / sum;
            I += If * wf;
            Q += Qf * wf;
            U += Uf * wf;
            V += Vf * wf;
        }
    }
}

////////////////////////////////////////////////////////////////////

void FragmentDustMixDecorator::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // build the cumulative distribution corresponding to the scattering opacities for each fragment
    Array Xv;
    NR::cdf(Xv, _numFrags,
            [this, lambda, state, pp](int f) { return WEIGHT(f) * _fragments[f]->opacitySca(lambda, state, pp); });

    // randomly select a fragment
    int f = NR::locateClip(Xv, random()->uniform());

    // actually perform the scattering event for this fragment
    _fragments[f]->performScattering(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

DisjointWavelengthGrid* FragmentDustMixDecorator::emissionWavelengthGrid() const
{
    return find<Configuration>()->dustEmissionWLG();
}

////////////////////////////////////////////////////////////////////

Array FragmentDustMixDecorator::emissivity(const Array& Jv) const
{
    return _dustMix->emissivity(Jv);
}

////////////////////////////////////////////////////////////////////

Array FragmentDustMixDecorator::emissionSpectrum(const MaterialState* state, const Array& Jv) const
{
    Array ev = WEIGHT(0) * _fragments[0]->emissivity(Jv);
    for (int f = 1; f != _numFrags; ++f) ev += WEIGHT(f) * _fragments[f]->emissivity(Jv);
    return state->numberDensity() * ev;
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::indicativeTemperature(const MaterialState* state, const Array& Jv) const
{
    double sumwT = 0.;
    double sumw = 0.;

    for (int f = 0; f != _numFrags; ++f)
    {
        double w = WEIGHT(f);
        double T = _fragments[f]->indicativeTemperature(nullptr, Jv);  // material state is not used by dust mixes
        sumwT += w * T;
        sumw += w;
    }
    return sumw > 0. ? sumwT / sumw : 0.;
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

double FragmentDustMixDecorator::populationBulkDensity(int f) const
{
    return _fragments[f]->populationBulkDensity(0);
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

double FragmentDustMixDecorator::populationTemperature(int f, const Array& Jv) const
{
    return _fragments[f]->indicativeTemperature(nullptr, Jv);  // material state is not used by dust mixes
}

////////////////////////////////////////////////////////////////////

bool FragmentDustMixDecorator::populationIsGraphite(int f) const
{
    return _fragments[f]->isGraphite();
}

////////////////////////////////////////////////////////////////////

double FragmentDustMixDecorator::populationGrainRadius(int f) const
{
    return _fragments[f]->grainRadius();
}

////////////////////////////////////////////////////////////////////
