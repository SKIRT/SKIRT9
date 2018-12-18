/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EquilibriumDustTemperatureCalculator.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "NR.hpp"
#include "PlanckFunction.hpp"

////////////////////////////////////////////////////////////////////

void EquilibriumDustTemperatureCalculator::precalculate(SimulationItem* item,
                                                        const Array& lambdav, const Array& sigmaabsv)
{
    // get the index of the bin being added
    int b = _rfsigmaabsvv.size();

    // perform initialization that needs to happen only once
    if (!b)
    {
        // obtain the simulation's radiation field wavelength grid
        _radiationFieldWLG = item->find<Configuration>()->radiationFieldWLG();
        // and prepare it for immediate use
        _radiationFieldWLG->setup();

        // build the temperature grid on which we store the Planck-integrated absorption cross sections
        NR::buildPowerLawGrid(_Tv, 0., 5000., 1000, 500.);
    }

    // allocate memory for the new bin
    int numWavelengths = _radiationFieldWLG->numBins();
    _rfsigmaabsvv.emplace_back(numWavelengths);
    int numTemperatures = _Tv.size();
    _planckabsvv.emplace_back(numTemperatures);

    // interpolate the absorption cross sections on the radiation field wavelength grid
    for (int ell=0; ell!=numWavelengths; ++ell)
    {
        double lambda = _radiationFieldWLG->wavelength(ell);
        _rfsigmaabsvv[b][ell] = NR::value<NR::interpolateLogLog>(lambda, lambdav, sigmaabsv);
    }

    // calculate the Planck-integrated absorption cross sections on the temperature grid
    int numLambda = lambdav.size();
    for (int p=1; p!=numTemperatures; ++p)   // leave value for p==0 at zero
    {
        PlanckFunction B(_Tv[p]);
        double planckabs = 0.;
        for (int ell=1; ell!=numLambda; ++ell) // skip the first wavelength so we can determine a bin width
        {
            double lambda = lambdav[ell];
            double dlambda = lambdav[ell] - lambdav[ell-1];
            planckabs += sigmaabsv[ell] * B(lambda) * dlambda;
        }
        _planckabsvv[b][p] = planckabs;
    }
}

////////////////////////////////////////////////////////////////////

size_t EquilibriumDustTemperatureCalculator::allocatedBytes() const
{
    size_t allocatedSize = 0;
    allocatedSize += _Tv.size();
    if (!_rfsigmaabsvv.empty()) allocatedSize += _rfsigmaabsvv.size() * _rfsigmaabsvv[0].size();
    if (!_planckabsvv.empty()) allocatedSize += _planckabsvv.size() * _planckabsvv[0].size();
    return allocatedSize * sizeof(double);
}

////////////////////////////////////////////////////////////////////

int EquilibriumDustTemperatureCalculator::numBins() const
{
    return _rfsigmaabsvv.size();
}

////////////////////////////////////////////////////////////////////

double EquilibriumDustTemperatureCalculator::equilibriumTemperature(int b, const Array& Jv) const
{
    // integrate the input side of the energy balance equation
    int numWavelengths = _radiationFieldWLG->numBins();
    double inputabs = 0.;
    for (int ell=0; ell!=numWavelengths; ++ell)
    {
        inputabs += _rfsigmaabsvv[b][ell] * Jv[ell] * _radiationFieldWLG->effectiveWidth(ell);
    }

    // find the temperature corresponding to this amount of emission on the output side of the equation
    return NR::clampedValue<NR::interpolateLinLin>(inputabs, _planckabsvv[b], _Tv);
}

////////////////////////////////////////////////////////////////////
