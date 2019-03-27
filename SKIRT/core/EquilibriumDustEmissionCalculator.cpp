/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EquilibriumDustEmissionCalculator.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "NR.hpp"
#include "PlanckFunction.hpp"

////////////////////////////////////////////////////////////////////

void EquilibriumDustEmissionCalculator::precalculate(SimulationItem* item,
                                                     const Array& lambdav, const Array& sigmaabsv)
{
    // get the index of the bin being added
    int b = _rfsigmaabsvv.size();

    // perform initialization that needs to happen only once
    if (!b)
    {
        auto config = item->find<Configuration>();

        // obtain the simulation's radiation field wavelength grid
        auto radiationFieldWLG = config->radiationFieldWLG();
        radiationFieldWLG->setup();
        _rflambdav.resize(radiationFieldWLG->numBins());
        _rfdlambdav.resize(radiationFieldWLG->numBins());
        for (int k=0; k!=radiationFieldWLG->numBins(); ++k)
        {
            _rflambdav[k] = radiationFieldWLG->wavelength(k);
            _rfdlambdav[k] = radiationFieldWLG->effectiveWidth(k);
        }

        // obtain the simulation's dust emission wavelength grid, if there is one
        if (config->hasDustEmission())
        {
            auto dustEmissionWLG = config->dustEmissionWLG();
            dustEmissionWLG->setup();
            _emlambdav.resize(dustEmissionWLG->numBins());
            for (int ell=0; ell!=dustEmissionWLG->numBins(); ++ell)
            {
                _emlambdav[ell] = dustEmissionWLG->wavelength(ell);
            }
        }

        // build the temperature grid on which we store the Planck-integrated absorption cross sections
        NR::buildPowerLawGrid(_Tv, 0., 5000., 1000, 500.);
    }

    // interpolate the absorption cross sections on the radiation field wavelength grid
    _rfsigmaabsvv.emplace_back(NR::resample<NR::interpolateLogLog>(_rflambdav, lambdav, sigmaabsv));

    // interpolate the absorption cross sections on the dust emission field wavelength grid, if there is one
    if (_emlambdav.size())
    {
        _emsigmaabsvv.emplace_back(NR::resample<NR::interpolateLogLog>(_emlambdav, lambdav, sigmaabsv));
    }

    // calculate the Planck-integrated absorption cross sections on the temperature grid
    int numLambda = lambdav.size();
    int numT = _Tv.size();
    _planckabsvv.emplace_back(numT);
    for (int p=1; p!=numT; ++p)   // leave value for p==0 at zero
    {
        PlanckFunction B(_Tv[p]);
        double planckabs = 0.;
        for (int j=1; j!=numLambda; ++j) // skip the first wavelength so we can determine a bin width
        {
            double lambda = lambdav[j];
            double dlambda = lambdav[j] - lambdav[j-1];
            planckabs += sigmaabsv[j] * B(lambda) * dlambda;
        }
        _planckabsvv[b][p] = planckabs;
    }
}

////////////////////////////////////////////////////////////////////

size_t EquilibriumDustEmissionCalculator::allocatedBytes() const
{
    size_t allocatedSize = 0;
    allocatedSize += _rflambdav.size();
    allocatedSize += _rfdlambdav.size();
    allocatedSize += _emlambdav.size();
    allocatedSize += _Tv.size();
    if (!_rfsigmaabsvv.empty()) allocatedSize += _rfsigmaabsvv.size() * _rfsigmaabsvv[0].size();
    if (!_emsigmaabsvv.empty()) allocatedSize += _emsigmaabsvv.size() * _emsigmaabsvv[0].size();
    if (!_planckabsvv.empty()) allocatedSize += _planckabsvv.size() * _planckabsvv[0].size();
    return allocatedSize * sizeof(double);
}

////////////////////////////////////////////////////////////////////

int EquilibriumDustEmissionCalculator::numBins() const
{
    return _rfsigmaabsvv.size();
}

////////////////////////////////////////////////////////////////////

double EquilibriumDustEmissionCalculator::equilibriumTemperature(int b, const Array& Jv) const
{
    // integrate the input side of the energy balance equation
    double inputabs = (_rfsigmaabsvv[b] * Jv * _rfdlambdav).sum();

    // find the temperature corresponding to this amount of emission on the output side of the equation
    if (inputabs > 0.) return NR::clampedValue<NR::interpolateLinLin>(inputabs, _planckabsvv[b], _Tv);
    else return 0.;
}

////////////////////////////////////////////////////////////////////

Array EquilibriumDustEmissionCalculator::emissivity(const Array& Jv) const
{
    int numWavelengths = _emlambdav.size();
    int numBins = _rfsigmaabsvv.size();

    Array ev(numWavelengths);
    for (int b=0; b!=numBins; ++b)
    {
        double T = equilibriumTemperature(b, Jv);
        PlanckFunction B(T);
        for (int ell=0; ell<numWavelengths; ell++)
        {
            ev[ell] += _emsigmaabsvv[b][ell] * B(_emlambdav[ell]);
        }
    }
    return ev;
}

////////////////////////////////////////////////////////////////////
