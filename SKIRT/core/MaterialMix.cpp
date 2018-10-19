/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MaterialMix.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "NR.hpp"
#include "PlanckFunction.hpp"
#include "Random.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

void MaterialMix::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _random = find<Random>();
}

////////////////////////////////////////////////////////////////////

void MaterialMix::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // precalculate information to accelerate solving the energy balance equation for the temperature;
    // this calculation is relevant only if the simulation tracks the radiation field
    if (find<Configuration>()->hasRadiationField())
    {
        // energy input side
        {
            // cache the simulation's radiation field wavelength grid
            _radiationFieldWLG = find<Configuration>()->radiationFieldWLG();

            // prepare it for immediate use
            _radiationFieldWLG->setup();

            // cache the absorption cross sections on the above wavelength grid
            int numWavelengths = _radiationFieldWLG->numBins();
            _sigmaabsv.resize(numWavelengths);
            for (int ell=0; ell!=numWavelengths; ++ell)
                _sigmaabsv[ell] = sectionAbs(_radiationFieldWLG->wavelength(ell));
        }

        // energy output side
        {
            // the temperature grid on which we store the Planck-integrated absorption cross sections
            const int numTemperatures = 1001;
            NR::buildPowerLawGrid(_Tv, 0., 5000., numTemperatures-1, 500.);

            // the wavelength grid over which to perform the integration
            const int numWavelengths = 3001;
            Array lambdav;
            NR::buildLogGrid(lambdav, 0.1e-6, 10000e-6, numWavelengths-1);

            // the absorption cross sections of this material on the above wavelength grid
            Array sigmaabsv(numWavelengths);
            for (int ell=0; ell!=numWavelengths; ++ell) sigmaabsv[ell] = sectionAbs(lambdav[ell]);

            // the Planck-integrated absorption cross sections on the above temperature grid
            _planckabsv.resize(numTemperatures);
            for (int p=1; p!=numTemperatures; ++p)   // leave value for p==0 at zero
            {
                PlanckFunction B(_Tv[p]);
                double planckabs = 0.;
                for (int ell=1; ell!=numWavelengths; ++ell) // skip the first wavelength so we can determine a bin width
                {
                    double lambda = lambdav[ell];
                    double dlambda = lambdav[ell] - lambdav[ell-1];
                    planckabs += sigmaabsv[ell] * B(lambda) * dlambda;
                }
                _planckabsv[p] = planckabs;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

double MaterialMix::asymmpar(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::phaseFunctionValueForCosine(double /*lambda*/, double /*costheta*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::generateCosineFromPhaseFunction(double /*lambda*/) const
{
    return 2.*random()->uniform() - 1.;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::phaseFunctionValue(double /*lambda*/, double /*theta*/, double /*phi*/,
                                       const StokesVector* /*sv*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> MaterialMix::generateAnglesFromPhaseFunction(double /*lambda*/,
                                                                       const StokesVector* /*sv*/) const
{
    return std::make_pair(0.,0.);
}

////////////////////////////////////////////////////////////////////

void MaterialMix::applyMueller(double /*lambda*/, double /*theta*/, StokesVector* /*sv*/) const
{
}

////////////////////////////////////////////////////////////////////

double MaterialMix::equilibriumTemperature(const Array& Jv) const
{
    // integrate the input side of the energy balance equation
    int numWavelengths = _radiationFieldWLG->numBins();
    double inputabs = 0.;
    for (int ell=0; ell!=numWavelengths; ++ell)
    {
        inputabs += _sigmaabsv[ell] * Jv[ell] * _radiationFieldWLG->effectiveWidth(ell);
    }

    // find the temperature corresponding to this amount of emission on the output side of the equation
    int p = NR::locateClip(_planckabsv, inputabs);
    return NR::interpolateLinLin(inputabs, _planckabsv[p], _planckabsv[p+1], _Tv[p], _Tv[p+1]);
}

////////////////////////////////////////////////////////////////////
