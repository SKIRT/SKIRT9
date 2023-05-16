/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumSystem.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DensityInCellInterface.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "LockFree.hpp"
#include "Log.hpp"
#include "LyaUtils.hpp"
#include "MaterialMix.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PathSegmentGenerator.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "ShortArray.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of cell densities calculated between two invocations of infoIfElapsed()
    const size_t logProgressChunkSize = 10000;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This helper class performs spatial sampling of the medium properties within a given cell and
    // for a given medium component. The number of samples to be taken can be specified separately
    // for the density and for all other properties (velocity, magnetic field, metallicity, ...).
    // The main reason for bundling this functionality in its own class/object is to allow caching
    // the random sampling positions with the corresponding sampled density values in a cell.
    // The constructor takes arguments that depend only on the simulation configuration and thus
    // do not vary between cells. This allows consecutively reusing the same object for multiple
    // cells, avoiding memory allocation/deallocation for each cell.
    class PropertySampler
    {
    private:
        vector<Medium*>& _media;
        SpatialGrid* _grid;
        int _numDensitySamples;
        int _numPropertySamples;
        int _numMedia;
        int _cellIndex;
        Position _center;
        vector<Position> _positions;  // indexed on n
        Table<2> _densities;          // indexed on h and n

    public:
        // The constructor allocates the required memory depending on the sampling configuration.
        // The number of samples for each type (density or other properties) can be:
        //   0: the sampling function(s) will never be called
        //   1: "sample" at the cell center only
        //  >1: sample at the specified nr of random points in the cell
        PropertySampler(vector<Medium*>& media, SpatialGrid* grid, int numDensitySamples, int numPropertySamples)
            : _media(media), _grid(grid), _numDensitySamples(numDensitySamples),
              _numPropertySamples(numPropertySamples), _numMedia(_media.size()), _cellIndex(0)
        {
            if (_numDensitySamples > 1 || _numPropertySamples > 1)
            {
                int numSamples = max(_numDensitySamples, _numPropertySamples);
                _positions.resize(numSamples);
                _densities.resize(_numMedia, numSamples);
            }
        }

        // This function should be called with a spatial cell index before calling any of the other
        // functions for that cell. It initializes the required sampling positions and obtains the
        // corresponding sampled density values.
        void prepareForCell(int cellIndex)
        {
            _cellIndex = cellIndex;

            // if we need center "sampling", get the cell center
            if (_numDensitySamples == 1 || _numPropertySamples == 1)
            {
                _center = _grid->centralPositionInCell(_cellIndex);
            }

            // if we need random sampling, draw the positions and get the corresponding densities
            int numSamples = _positions.size();
            for (int n = 0; n != numSamples; ++n)
            {
                _positions[n] = _grid->randomPositionInCell(_cellIndex);
                for (int h = 0; h != _numMedia; ++h) _densities(h, n) = _media[h]->numberDensity(_positions[n]);
            }
        }

        // Returns the central or average density in the cell for the given medium component.
        double density(int h)
        {
            if (_numDensitySamples == 1) return _media[h]->numberDensity(_center);
            double sum = 0.;
            for (int n = 0; n != _numDensitySamples; ++n) sum += _densities(h, n);
            return sum / _numDensitySamples;
        }

        // Returns the central or density-averaged velocity in the cell for the given medium component.
        Vec bulkVelocity(int h)
        {
            if (_numPropertySamples == 1) return _media[h]->bulkVelocity(_center);
            double nsum = 0;
            Vec vsum;
            for (int n = 0; n != _numPropertySamples; ++n)
            {
                nsum += _densities(h, n);
                vsum += _densities(h, n) * _media[h]->bulkVelocity(_positions[n]);
            }
            return nsum > 0. ? vsum / nsum : Vec();
        }

        // Returns the central or average magnetic field in the cell for the given medium component.
        // Note the magnetic field average is NOT weighed with density.
        Vec magneticField(int h)
        {
            if (_numPropertySamples == 1) return _media[h]->magneticField(_center);
            Vec Bsum;
            for (int n = 0; n != _numPropertySamples; ++n) Bsum += _media[h]->magneticField(_positions[n]);
            return Bsum / _numPropertySamples;
        }

        // Returns the central or density-averaged metallicity in the cell for the given medium component.
        double metallicity(int h)
        {
            if (_numPropertySamples == 1) return _media[h]->metallicity(_center);
            double nsum = 0;
            double Zsum = 0.;
            for (int n = 0; n != _numPropertySamples; ++n)
            {
                nsum += _densities(h, n);
                Zsum += _densities(h, n) * _media[h]->metallicity(_positions[n]);
            }
            return nsum > 0. ? Zsum / nsum : 0.;
        }

        // Returns the central or density-averaged temperature in the cell for the given medium component.
        double temperature(int h)
        {
            if (_numPropertySamples == 1) return _media[h]->temperature(_center);
            double nsum = 0;
            double Tsum = 0.;
            for (int n = 0; n != _numPropertySamples; ++n)
            {
                nsum += _densities(h, n);
                Tsum += _densities(h, n) * _media[h]->temperature(_positions[n]);
            }
            return nsum > 0. ? Tsum / nsum : 0.;
        }

        // Returns the central or density-averaged custom parameter values in the cell for the given medium component.
        void parameters(int h, Array& params)
        {
            // always sample at the center to make sure that the output array has the proper size
            _media[h]->parameters(_center, params);
            // if center-sampling is requested, we're done
            if (_numPropertySamples == 1) return;
            // otherwise, clear the output values and accumulate the samples
            params = 0.;
            double nsum = 0;
            Array sample;
            for (int n = 0; n != _numPropertySamples; ++n)
            {
                nsum += _densities(h, n);
                _media[h]->parameters(_positions[n], sample);
                params += _densities(h, n) * sample;
            }
            params /= nsum;
            return;
        }
    };
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();
    _config = find<Configuration>();

    _numCells = _grid->numCells();
    if (_numCells < 1) throw FATALERROR("The spatial grid must have at least one cell");
    _numMedia = _media.size();
    size_t allocatedBytes = 0;

    // ----- allocate memory for the medium state -----

    // basic configuration
    _state.initConfiguration(_numCells, _numMedia);

    // common state variables
    vector<StateVariable> variables;
    variables.emplace_back(StateVariable::volume());
    if (_config->hasMovingMedia()) variables.emplace_back(StateVariable::bulkVelocity());
    if (_config->hasMagneticField()) variables.emplace_back(StateVariable::magneticField());
    _state.initCommonStateVariables(variables);

    // specific state variables
    for (auto medium : _media) _state.initSpecificStateVariables(medium->mix()->specificStateVariableInfo());

    // finalize
    allocatedBytes += _state.initAllocate() * sizeof(double);

    // ----- allocate memory for the radiation field -----

    if (_config->hasRadiationField())
    {
        _wavelengthGrid = _config->radiationFieldWLG();
        _rf1.resize(_numCells, _wavelengthGrid->numBins());
        allocatedBytes += _rf1.size() * sizeof(double);

        if (_config->hasSecondaryRadiationField())
        {
            _rf2.resize(_numCells, _wavelengthGrid->numBins());
            _rf2c.resize(_numCells, _wavelengthGrid->numBins());
            allocatedBytes += 2 * _rf2.size() * sizeof(double);
        }
    }

    // ----- cache info on the dust emission wavelength grid -----

    if (_config->hasDustEmission())
    {
        _numDustEmissionWavelengths = _config->dustEmissionWLG()->extlambdav().size();
    }

    // ----- obtain the material mix pointers -----

    if (_config->hasVariableMedia())
    {
        // we need a separate material mix pointer for each spatial cell (per medium component)
        _mixv.resize(_numCells * _numMedia);
        _mixPerCell = true;
        for (int m = 0; m != _numCells; ++m)
        {
            Position bfr = _grid->centralPositionInCell(m);
            for (int h = 0; h != _numMedia; ++h) _mixv[m * _numMedia + h] = _media[h]->mix(bfr);
        }
    }
    else
    {
        // the material mix pointer is identical for all spatial cells (per medium component)
        _mixv.resize(_numMedia);
        for (int h = 0; h != _numMedia; ++h) _mixv[h] = _media[h]->mix();
    }
    allocatedBytes += _mixv.size() * sizeof(MaterialMix*);

    // cache a list of medium component indices for each material type
    for (int h = 0; h != _numMedia; ++h)
    {
        switch (mix(0, h)->materialType())
        {
            case MaterialMix::MaterialType::Dust: _dust_hv.push_back(h); break;
            case MaterialMix::MaterialType::Gas: _gas_hv.push_back(h); break;
            case MaterialMix::MaterialType::Electrons: _elec_hv.push_back(h); break;
        }
        switch (mix(0, h)->hasDynamicMediumState())
        {
            case MaterialMix::DynamicStateType::None: break;
            case MaterialMix::DynamicStateType::Primary: _pdms_hv.push_back(h); break;
            case MaterialMix::DynamicStateType::Secondary: _sdms_hv.push_back(h); break;
            case MaterialMix::DynamicStateType::PrimaryIfMergedIterations:
                if (_config->hasMergedIterations())
                    _pdms_hv.push_back(h);
                else
                    _sdms_hv.push_back(h);
                break;
        }
    }

    // ----- inform user about allocated memory -----

    log->info(typeAndName() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");

    // ----- calculate medium properties parallelized on spatial cells -----

    log->info("Calculating medium properties for " + std::to_string(_numCells) + " cells...");
    auto dic = _grid->interface<DensityInCellInterface>(0, false);  // optional fast-track interface for densities
    log->infoSetElapsed(_numCells);
    parfac->parallelDistributed()->call(_numCells, [this, log, dic](size_t firstIndex, size_t numIndices) {
        // construct a property sampler to be shared by all cells handled in the loop below
        bool hasProperty = _config->hasMovingMedia() || _config->hasMagneticField();
        for (int h = 0; h != _numMedia; ++h)
            if (_media[h]->hasMetallicity() || _media[h]->hasTemperature() || _media[h]->hasParameters())
                hasProperty = true;
        PropertySampler sampler(_media, _grid, dic ? 0 : _config->numDensitySamples(),
                                hasProperty ? _config->numPropertySamples() : 0);

        // loop over cells
        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
            {
                // prepare the sampler for this cell (i.e. store the relevant positions and density samples)
                sampler.prepareForCell(m);

                // volume
                _state.setVolume(m, _grid->volume(m));

                // density: use optional fast-track interface or sample within the cell
                if (dic)
                {
                    for (int h = 0; h != _numMedia; ++h) _state.setNumberDensity(m, h, dic->numberDensity(h, m));
                }
                else
                {
                    for (int h = 0; h != _numMedia; ++h) _state.setNumberDensity(m, h, sampler.density(h));
                }

                // bulk velocity: aggregate values sampled for each component according to specfied policy
                if (_config->hasMovingMedia())
                {
                    // accumulate total density in the cell
                    double n = 0.;
                    for (int h = 0; h != _numMedia; ++h) n += _state.numberDensity(m, h);

                    // leave bulk velocity at zero if cell has no material
                    if (n > 0.)
                    {
                        Vec v;
                        switch (_samplingOptions->aggregateVelocity())
                        {
                            case SamplingOptions::AggregatePolicy::Average:
                                for (int h = 0; h != _numMedia; ++h)
                                {
                                    v += _state.numberDensity(m, h) * sampler.bulkVelocity(h);
                                }
                                v /= n;
                                break;
                            case SamplingOptions::AggregatePolicy::Maximum:
                                for (int h = 0; h != _numMedia; ++h)
                                {
                                    if (_state.numberDensity(m, h) > 0)
                                    {
                                        Vec w = sampler.bulkVelocity(h);
                                        if (w.norm2() > v.norm2()) v = w;
                                    }
                                }
                                break;
                            case SamplingOptions::AggregatePolicy::First:
                                for (int h = 0; h != _numMedia; ++h)
                                {
                                    if (_media[h]->hasVelocity())
                                    {
                                        if (_state.numberDensity(m, h) > 0) v = sampler.bulkVelocity(h);
                                        break;
                                    }
                                    break;
                                }
                        }
                        _state.setBulkVelocity(m, v);
                    }
                }

                // magnetic field: retrieve value sampled from medium component that specifies it
                if (_config->hasMagneticField())
                {
                    _state.setMagneticField(m, sampler.magneticField(_config->magneticFieldMediumIndex()));
                }

                // specific state variables other than density:
                // retrieve value sampled from corresponding medium component
                for (int h = 0; h != _numMedia; ++h)
                {
                    double Z = _media[h]->hasMetallicity() ? sampler.metallicity(h) : -1.;
                    double T = _media[h]->hasTemperature() ? sampler.temperature(h) : -1.;
                    Array params;
                    if (_media[h]->hasParameters()) sampler.parameters(h, params);
                    MaterialState mst(_state, m, h);
                    mix(m, h)->initializeSpecificState(&mst, Z, T, params);
                }
            }
            log->infoIfElapsed("Calculated medium properties: ", currentChunkSize);
            firstIndex += currentChunkSize;
            numIndices -= currentChunkSize;
        }
    });

    // communicate the calculated states across multiple processes, if needed
    _state.initCommunicate();

    log->info("Done calculating cell densities");
}

////////////////////////////////////////////////////////////////////

int MediumSystem::dimension() const
{
    int result = 1;
    for (auto medium : _media) result = max(result, medium->dimension());
    return result;
}

////////////////////////////////////////////////////////////////////

int MediumSystem::gridDimension() const
{
    return _grid->dimension();
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numMedia() const
{
    return _numMedia;
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numCells() const
{
    return _numCells;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* MediumSystem::mix(int m, int h) const
{
    return _mixPerCell ? _mixv[m * _numMedia + h] : _mixv[h];
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::hasMaterialType(MaterialMix::MaterialType type) const
{
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) return true;
    return false;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::isMaterialType(MaterialMix::MaterialType type, int h) const
{
    return mix(0, h)->materialType() == type;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustMassDensity(Position bfr) const
{
    double result = 0.;
    for (int h : _dust_hv) result += _media[h]->massDensity(bfr);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::electronNumberDensity(Position bfr) const
{
    double result = 0.;
    for (int h : _elec_hv) result += _media[h]->numberDensity(bfr);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::gasNumberDensity(Position bfr) const
{
    double result = 0.;
    for (int h : _gas_hv) result += _media[h]->numberDensity(bfr);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::volume(int m) const
{
    return _state.volume(m);
}

////////////////////////////////////////////////////////////////////

Vec MediumSystem::bulkVelocity(int m) const
{
    return _state.bulkVelocity(m);
}

////////////////////////////////////////////////////////////////////

Vec MediumSystem::magneticField(int m) const
{
    return _state.magneticField(m);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::massDensity(int m) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h) result += massDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustMassDensity(int m) const
{
    double result = 0.;
    for (int h : _dust_hv) result += massDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::electronNumberDensity(int m) const
{
    double result = 0.;
    for (int h : _elec_hv) result += numberDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::gasNumberDensity(int m) const
{
    double result = 0.;
    for (int h : _gas_hv) result += numberDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::numberDensity(int m, int h) const
{
    return _state.numberDensity(m, h);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::massDensity(int m, int h) const
{
    return _state.numberDensity(m, h) * mix(m, h)->mass();
}

////////////////////////////////////////////////////////////////////

double MediumSystem::metallicity(int m, int h) const
{
    return _state.metallicity(m, h);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::temperature(int m, int h) const
{
    return _state.temperature(m, h);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::custom(int m, int h, int i) const
{
    return _state.custom(m, h, i);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacityAbs(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacitySca(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacityExt(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) result += opacityAbs(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) result += opacityExt(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h) result += opacityExt(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, const PhotonPacket* pp) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        result += mix(m, h)->opacityAbs(lambda, &mst, pp);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, const PhotonPacket* pp) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        result += mix(m, h)->opacitySca(lambda, &mst, pp);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, const PhotonPacket* pp) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        result += mix(m, h)->opacityExt(lambda, &mst, pp);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::perceivedWavelengthForScattering(const PhotonPacket* pp) const
{
    if (_config->hasMovingMedia())
        return pp->perceivedWavelength(_state.bulkVelocity(pp->interactionCellIndex()),
                                       _config->hubbleExpansionRate() * pp->interactionDistance());
    else
        return pp->wavelength();
}

////////////////////////////////////////////////////////////////////

double MediumSystem::albedoForScattering(const PhotonPacket* pp) const
{
    double lambda = perceivedWavelengthForScattering(pp);
    int m = pp->interactionCellIndex();
    if (m < 0) throw FATALERROR("Cannot locate photon packet interaction point");

    double ksca = 0.;
    double kext = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        ksca += mix(m, h)->opacitySca(lambda, &mst, pp);
        kext += mix(m, h)->opacityExt(lambda, &mst, pp);
    }
    return kext > 0. ? ksca / kext : 0.;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::weightsForScattering(ShortArray& wv, double lambda, const PhotonPacket* pp) const
{
    // resize the target array
    wv.resize(_numMedia);

    // for a single component, the weight is trivial
    if (_numMedia == 1)
    {
        wv[0] = 1.;
        return true;
    }

    // locate the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the weights and their sum
    double sum = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        wv[h] = mix(m, h)->opacitySca(lambda, &mst, pp);
        sum += wv[h];
    }

    // normalize the weights
    if (sum > 0.)
    {
        for (int h = 0; h != _numMedia; ++h) wv[h] /= sum;
        return true;
    }

    // none of the media scatter
    return false;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::peelOffScattering(const ShortArray& wv, double lambda, Direction bfkobs, Direction bfky,
                                     PhotonPacket* pp, PhotonPacket* ppp) const
{
    // get the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the weighted sum of the effects on the Stokes vector for all media
    // we assume the wavelength does not change
    double I = 0., Q = 0., U = 0., V = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        // skip media that don't scatter this photon packet
        if (wv[h] > 0.)
        {
            double Ih = 0., Qh = 0., Uh = 0., Vh = 0.;
            MaterialState mst(_state, m, h);
            pp->setScatteringComponent(h);
            mix(m, h)->peeloffScattering(Ih, Qh, Uh, Vh, lambda, bfkobs, bfky, &mst, pp);
            I += Ih * wv[h];
            Q += Qh * wv[h];
            U += Uh * wv[h];
            V += Vh * wv[h];
        }
    }

    // pass the result to the peel-off photon packet
    ppp->launchScatteringPeelOff(pp, bfkobs, _state.bulkVelocity(m), lambda, I);
    if (_config->hasPolarization()) ppp->setPolarized(I, Q, U, V, pp->normal());
}

////////////////////////////////////////////////////////////////////

void MediumSystem::peelOffScattering(int h, double w, double lambda, Direction bfkobs, Direction bfky, PhotonPacket* pp,
                                     PhotonPacket* ppp) const
{
    // get the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the effects on the Stokes vector and on the wavelength for this medium component
    double I = 0., Q = 0., U = 0., V = 0.;
    MaterialState mst(_state, m, h);
    pp->setScatteringComponent(h);
    mix(m, h)->peeloffScattering(I, Q, U, V, lambda, bfkobs, bfky, &mst, pp);

    // pass the result to the peel-off photon packet
    ppp->launchScatteringPeelOff(pp, bfkobs, _state.bulkVelocity(m), lambda, I * w);
    if (_config->hasPolarization()) ppp->setPolarized(I, Q, U, V, pp->normal());
}

////////////////////////////////////////////////////////////////////

void MediumSystem::simulateScattering(Random* random, PhotonPacket* pp) const
{
    // locate the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the perceived wavelength in the cell
    double lambda = perceivedWavelengthForScattering(pp);

    // select a medium component within that cell
    int h = 0;
    if (_numMedia > 1)
    {
        // build the cumulative distribution corresponding to the scattering opacities
        Array Xv;
        NR::cdf(Xv, _numMedia, [this, lambda, pp, m](int h) {
            MaterialState mst(_state, m, h);
            return mix(m, h)->opacitySca(lambda, &mst, pp);
        });

        // randomly select an index
        h = NR::locateClip(Xv, random->uniform());
    }

    // actually perform the scattering event for this cell and medium component
    MaterialState mst(_state, m, h);
    pp->setScatteringComponent(h);
    mix(m, h)->performScattering(lambda, &mst, pp);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function returns a thread-local instance of the path segment generator for the specified grid
    // that is initialized to the starting position and direction of the specified path.
    // Providing a thread-local instance avoids creating a new generator for each use.
    PathSegmentGenerator* getPathSegmentGenerator(SpatialGrid* grid, const SpatialGridPath* path)
    {
        thread_local SpatialGrid* t_grid{nullptr};
        thread_local std::unique_ptr<PathSegmentGenerator> t_generator;

        if (grid != t_grid)
        {
            t_grid = grid;
            t_generator = grid->createPathSegmentGenerator();
        }
        t_generator->start(path);
        return t_generator.get();
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setExtinctionOpticalDepths(PhotonPacket* pp) const
{
    // determine and store the path segments in the photon packet
    auto generator = getPathSegmentGenerator(_grid, pp);
    pp->clear();
    while (generator->next())
    {
        pp->addSegment(generator->m(), generator->ds());
    }

    // calculate the cumulative optical depth and store it in the photon packet for each path segment
    double tau = 0.;

    // single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double section = mix(0, 0)->sectionExt(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0) tau += section * _state.numberDensity(segment.m(), 0) * segment.ds();
            segment.setOpticalDepth(tau);
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = mix(0, h)->sectionExt(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
                for (int h = 0; h != _numMedia; ++h)
                    tau += sectionv[h] * _state.numberDensity(segment.m(), h) * segment.ds();
            segment.setOpticalDepth(tau);
        }
    }

    // spatially variable cross sections
    else
    {
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(segment.m()),
                                                        _config->hubbleExpansionRate() * segment.s());
                tau += opacityExt(lambda, segment.m(), pp) * segment.ds();
            }
            segment.setOpticalDepth(tau);
        }
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setScatteringAndAbsorptionOpticalDepths(PhotonPacket* pp) const
{
    // determine and store the path segments in the photon packet
    auto generator = getPathSegmentGenerator(_grid, pp);
    pp->clear();
    while (generator->next())
    {
        pp->addSegment(generator->m(), generator->ds());
    }

    // calculate the cumulative optical depths and store them in the photon packet for each path segment
    double tauSca = 0.;
    double tauAbs = 0.;

    // single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double sectionSca = mix(0, 0)->sectionSca(pp->wavelength());
        double sectionAbs = mix(0, 0)->sectionAbs(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
            {
                double ns = _state.numberDensity(segment.m(), 0) * segment.ds();
                tauSca += sectionSca * ns;
                tauAbs += sectionAbs * ns;
            }
            segment.setOpticalDepth(tauSca, tauAbs);
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionScav(_numMedia);
        ShortArray sectionAbsv(_numMedia);
        for (int h = 0; h != _numMedia; ++h)
        {
            sectionScav[h] = mix(0, h)->sectionSca(pp->wavelength());
            sectionAbsv[h] = mix(0, h)->sectionAbs(pp->wavelength());
        }
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
                for (int h = 0; h != _numMedia; ++h)
                {
                    double ns = _state.numberDensity(segment.m(), h) * segment.ds();
                    tauSca += sectionScav[h] * ns;
                    tauAbs += sectionAbsv[h] * ns;
                }
            segment.setOpticalDepth(tauSca, tauAbs);
        }
    }

    // spatially variable cross sections
    else
    {
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(segment.m()),
                                                        _config->hubbleExpansionRate() * segment.s());
                tauSca += opacitySca(lambda, segment.m(), pp) * segment.ds();
                tauAbs += opacityAbs(lambda, segment.m(), pp) * segment.ds();
            }
            segment.setOpticalDepth(tauSca, tauAbs);
        }
    }
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::setInteractionPointUsingExtinction(PhotonPacket* pp, double tauinteract) const
{
    auto generator = getPathSegmentGenerator(_grid, pp);
    double tau = 0.;
    double s = 0.;

    // loop over the segments of the path until the interaction optical depth is reached or the path ends

    // --> single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double section = mix(0, 0)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0) tau += section * _state.numberDensity(m, 0) * generator->ds();
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tau)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tau0, tau, s0, s));
                return true;
            }
        }
    }

    // --> multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = mix(0, h)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
                for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * _state.numberDensity(m, h) * ds;
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tau)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tau0, tau, s0, s));
                return true;
            }
        }
    }

    // --> spatially variable cross sections
    else
    {
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(m), _config->hubbleExpansionRate() * s);
                tau += opacityExt(lambda, m, pp) * ds;
            }
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tau)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tau0, tau, s0, s));
                return true;
            }
        }
    }

    // the interaction point is outside of the path
    return false;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::setInteractionPointUsingScatteringAndAbsorption(PhotonPacket* pp, double tauinteract) const
{
    auto generator = getPathSegmentGenerator(_grid, pp);
    double tauSca = 0.;
    double tauAbs = 0.;
    double s = 0.;

    // loop over the segments of the path until the interaction optical depth is reached or the path ends

    // --> single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double sectionSca = mix(0, 0)->sectionSca(pp->wavelength());
        double sectionAbs = mix(0, 0)->sectionAbs(pp->wavelength());
        while (generator->next())
        {
            // remember the cumulative scattering optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tauSca0 = tauSca;
            double s0 = s;

            // calculate the cumulative optical depths and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0) tauSca += sectionSca * _state.numberDensity(m, 0) * generator->ds();
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tauSca)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tauSca0, tauSca, s0, s),
                                        tauinteract * sectionAbs / sectionSca);
                return true;
            }
        }
    }

    // --> multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionScav(_numMedia);
        ShortArray sectionAbsv(_numMedia);
        for (int h = 0; h != _numMedia; ++h)
        {
            sectionScav[h] = mix(0, h)->sectionSca(pp->wavelength());
            sectionAbsv[h] = mix(0, h)->sectionAbs(pp->wavelength());
        }
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tauSca0 = tauSca;
            double tauAbs0 = tauAbs;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                for (int h = 0; h != _numMedia; ++h)
                {
                    double ns = _state.numberDensity(m, h) * ds;
                    tauSca += sectionScav[h] * ns;
                    tauAbs += sectionAbsv[h] * ns;
                }
            }
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tauSca)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tauSca0, tauSca, s0, s),
                                        NR::interpolateLinLin(tauinteract, tauSca0, tauSca, tauAbs0, tauAbs));
                return true;
            }
        }
    }

    // --> spatially variable cross sections
    else
    {
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tauSca0 = tauSca;
            double tauAbs0 = tauAbs;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(m), _config->hubbleExpansionRate() * s);
                tauSca += opacitySca(lambda, m, pp) * ds;
                tauAbs += opacityAbs(lambda, m, pp) * ds;
            }
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tauSca)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tauSca0, tauSca, s0, s),
                                        NR::interpolateLinLin(tauinteract, tauSca0, tauSca, tauAbs0, tauAbs));
                return true;
            }
        }
    }

    // the interaction point is outside of the path
    return false;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::getExtinctionOpticalDepth(const PhotonPacket* pp, double distance) const
{
    // abort if the packet's contribution is zero to begin with
    double L = pp->luminosity();
    if (L <= 0) return std::numeric_limits<double>::infinity();

    // if extinction is always positive, determine the optical depth at which the packet's contribution becomes zero
    double taumax = _config->hasNegativeExtinction() ? std::numeric_limits<double>::infinity() : std::log(L) + 745;

    // determine the geometric details of the path and calculate the optical depth at the same time
    auto generator = getPathSegmentGenerator(_grid, pp);
    double tau = 0.;
    double s = 0.;

    // single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double section = mix(0, 0)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            if (generator->m() >= 0)
            {
                tau += section * _state.numberDensity(generator->m(), 0) * generator->ds();
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += generator->ds();
            if (s > distance) break;
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = mix(0, h)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * _state.numberDensity(m, h) * ds;
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += ds;
            if (s > distance) break;
        }
    }

    // spatially variable cross sections
    else
    {
        while (generator->next())
        {
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(m), _config->hubbleExpansionRate() * s);
                tau += opacityExt(lambda, m, pp) * ds;
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += ds;
            if (s > distance) break;
        }
    }

    return tau;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::getExtinctionOpticalDepth(const SpatialGridPath* path, double lambda,
                                               MaterialMix::MaterialType type) const
{
    // determine the geometric details of the path and calculate the optical depth at the same time
    auto generator = getPathSegmentGenerator(_grid, path);
    double tau = 0.;
    while (generator->next())
    {
        if (generator->m() >= 0) tau += opacityExt(lambda, generator->m(), type) * generator->ds();
    }
    return tau;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::clearRadiationField(bool primary)
{
    if (primary)
    {
        _rf1.setToZero();
        if (_rf2.size()) _rf2.setToZero();
    }
    else
    {
        _rf2c.setToZero();
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::storeRadiationField(bool primary, int m, int ell, double Lds)
{
    if (primary)
        LockFree::add(_rf1(m, ell), Lds);
    else
        LockFree::add(_rf2c(m, ell), Lds);
}

////////////////////////////////////////////////////////////////////

void MediumSystem::communicateRadiationField(bool primary)
{
    if (primary)
        ProcessManager::sumToAll(_rf1.data());
    else
    {
        ProcessManager::sumToAll(_rf2c.data());
        _rf2 = _rf2c;
    }
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> MediumSystem::totalDustAbsorbedLuminosity() const
{
    double Labs1 = 0.;
    double Labs2 = 0.;
    int numWavelengths = _wavelengthGrid->numBins();
    for (int ell = 0; ell != numWavelengths; ++ell)
    {
        double lambda = _wavelengthGrid->wavelength(ell);
        for (int m = 0; m != _numCells; ++m)
        {
            double opacity = opacityAbs(lambda, m, MaterialMix::MaterialType::Dust);
            Labs1 += opacity * _rf1(m, ell);
            Labs2 += opacity * _rf2(m, ell);
        }
    }
    return std::make_pair(Labs1, Labs2);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::radiationField(int m, int ell) const
{
    double rf = 0.;
    if (_rf1.size()) rf += _rf1(m, ell);
    if (_rf2.size()) rf += _rf2(m, ell);
    return rf;
}

////////////////////////////////////////////////////////////////////

Array MediumSystem::meanIntensity(int m) const
{
    int numWavelengths = _wavelengthGrid->numBins();
    Array Jv(numWavelengths);
    double factor = 1. / (4. * M_PI * _state.volume(m));
    for (int ell = 0; ell < numWavelengths; ell++)
    {
        Jv[ell] = radiationField(m, ell) * factor / _wavelengthGrid->effectiveWidth(ell);
    }
    return Jv;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeTemperature(int m, int h) const
{
    // get the radiation field, if available
    Array Jv;
    if (_config->hasRadiationField()) Jv = meanIntensity(m);

    // obtain the temperature for the requested component and cell
    if (massDensity(m, h) > 0.)
    {
        MaterialState mst(_state, m, h);
        return mix(m, h)->indicativeTemperature(&mst, Jv);
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeTemperature(int m, MaterialMix::MaterialType type) const
{
    // get the radiation field, if available
    Array Jv;
    if (_config->hasRadiationField()) Jv = meanIntensity(m);

    // obtain the temperature and weight for each component of the requested type
    double sumRhoT = 0.;
    double sumRho = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        if (isMaterialType(type, h))
        {
            double rho = massDensity(m, h);
            if (rho > 0.)
            {
                MaterialState mst(_state, m, h);
                double T = mix(m, h)->indicativeTemperature(&mst, Jv);
                sumRhoT += rho * T;
                sumRho += rho;
            }
        }
    }

    // return the average, if there is one
    return sumRho > 0. ? sumRhoT / sumRho : 0.;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeDustTemperature(int m) const
{
    return indicativeTemperature(m, MaterialMix::MaterialType::Dust);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeElectronTemperature(int m) const
{
    return indicativeTemperature(m, MaterialMix::MaterialType::Electrons);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeGasTemperature(int m) const
{
    return indicativeTemperature(m, MaterialMix::MaterialType::Gas);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustLuminosity(int m) const
{
    double Labs = 0.;
    int numWavelengths = _wavelengthGrid->numBins();
    for (int ell = 0; ell < numWavelengths; ell++)
    {
        Labs +=
            opacityAbs(_wavelengthGrid->wavelength(ell), m, MaterialMix::MaterialType::Dust) * radiationField(m, ell);
    }
    return Labs;
}

////////////////////////////////////////////////////////////////////

Array MediumSystem::dustEmissionSpectrum(int m) const
{
    const Array& Jv = meanIntensity(m);
    Array ev(_numDustEmissionWavelengths);
    for (int h : _dust_hv)
    {
        MaterialState mst(_state, m, h);
        ev += mix(m, h)->emissionSpectrum(&mst, Jv);
    }
    return ev;
}

////////////////////////////////////////////////////////////////////

Array MediumSystem::continuumEmissionSpectrum(int m, int h) const
{
    const Array& Jv = meanIntensity(m);
    MaterialState mst(_state, m, h);
    return mix(m, h)->emissionSpectrum(&mst, Jv);
}

////////////////////////////////////////////////////////////////////

Array MediumSystem::lineEmissionSpectrum(int m, int h) const
{
    const Array& Jv = meanIntensity(m);
    MaterialState mst(_state, m, h);
    return mix(m, h)->lineEmissionSpectrum(&mst, Jv);
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::updateDynamicStateRecipes()
{
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();
    auto& recipes = dynamicStateOptions()->recipes();

    // tell all recipes to begin the update cycle
    for (auto recipe : recipes) recipe->beginUpdate(_numCells);

    // update status for each cell
    std::vector<UpdateStatus> flags(_numCells);

    // loop over the spatial cells in parallel
    log->info("Updating medium state for " + std::to_string(_numCells) + " cells...");
    log->infoSetElapsed(_numCells);
    parfac->parallelDistributed()->call(_numCells, [this, log, &recipes, &flags](size_t firstIndex, size_t numIndices) {
        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
            {
                const Array& Jv = meanIntensity(m);

                // tell all recipes to update the state for all media in this cell
                for (auto recipe : recipes)
                {
                    for (int h = 0; h != _numMedia; ++h)
                    {
                        MaterialState mst(_state, m, h);
                        flags[m].update(recipe->update(&mst, Jv));
                    }
                }
            }
            log->infoIfElapsed("Updated medium state: ", currentChunkSize);
            firstIndex += currentChunkSize;
            numIndices -= currentChunkSize;
        }
    });

    // synchronize the updated state between processes
    int numUpdated, numNotConverged;
    std::tie(numUpdated, numNotConverged) = _state.synchronize(flags);

    // log statistics
    log->info("  Updated cells: " + std::to_string(numUpdated) + " out of " + std::to_string(_numCells) + " ("
              + StringUtils::toString(100. * numUpdated / _numCells, 'f', 2) + " %)");
    log->info("  Not converged: " + std::to_string(numNotConverged) + " out of " + std::to_string(_numCells) + " ("
              + StringUtils::toString(100. * numNotConverged / _numCells, 'f', 2) + " %)");

    // tell all recipes to end the update cycle and collect convergence info
    bool converged = true;
    for (auto recipe : recipes) converged &= recipe->endUpdate(_numCells, numUpdated, numNotConverged);
    return converged;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::updateDynamicStateMedia(bool primary)
{
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();

    // update status for each cell
    std::vector<UpdateStatus> flags(_numCells);

    // loop over the spatial cells in parallel
    log->info("Updating medium state for " + std::to_string(_numCells) + " cells...");
    log->infoSetElapsed(_numCells);
    parfac->parallelDistributed()->call(_numCells, [this, primary, log, &flags](size_t firstIndex, size_t numIndices) {
        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
            {
                const Array& Jv = meanIntensity(m);
                for (int h : (primary ? _pdms_hv : _sdms_hv))
                {
                    MaterialState mst(_state, m, h);
                    flags[m].update(mix(m, h)->updateSpecificState(&mst, Jv));
                }
            }
            log->infoIfElapsed("Updated medium state: ", currentChunkSize);
            firstIndex += currentChunkSize;
            numIndices -= currentChunkSize;
        }
    });

    // synchronize the updated state between processes
    int numUpdated, numNotConverged;
    std::tie(numUpdated, numNotConverged) = _state.synchronize(flags);

    // log statistics
    log->info("  Updated cells: " + std::to_string(numUpdated) + " out of " + std::to_string(_numCells) + " ("
              + StringUtils::toString(100. * numUpdated / _numCells, 'f', 2) + " %)");
    if (numNotConverged)
        log->info("  Not converged: " + std::to_string(numNotConverged) + " out of " + std::to_string(_numCells) + " ("
                  + StringUtils::toString(100. * numNotConverged / _numCells, 'f', 2) + " %)");

    // collect convergence info
    bool converged = true;
    for (int h : (primary ? _pdms_hv : _sdms_hv))
        converged &= mix(0, h)->isSpecificStateConverged(_numCells, numUpdated, numNotConverged);
    return converged;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::updatePrimaryDynamicMediumState()
{
    bool converged = true;
    if (_config->hasDynamicStateRecipes()) converged &= updateDynamicStateRecipes();
    if (_config->hasPrimaryDynamicStateMedia()) converged &= updateDynamicStateMedia(true);
    return converged;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::updateSecondaryDynamicMediumState()
{
    bool converged = true;
    if (_config->hasSecondaryDynamicStateMedia()) converged &= updateDynamicStateMedia(false);
    return converged;
}

////////////////////////////////////////////////////////////////////
