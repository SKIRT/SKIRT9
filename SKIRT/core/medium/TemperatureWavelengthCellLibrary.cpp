/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TemperatureWavelengthCellLibrary.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

int TemperatureWavelengthCellLibrary::numEntries() const
{
    return _numTemperatures * _numWavelengths;
}

////////////////////////////////////////////////////////////////////

namespace
{
    double indicativeDustWavelength(int m, const MediumSystem* ms, const DisjointWavelengthGrid* wavelengthGrid)
    {
        const Array& Jv = ms->meanIntensity(m);
        double sumtop = 0.;
        double sumbot = 0.;
        int numWavelengths = wavelengthGrid->numBins();
        for (int ell = 0; ell != numWavelengths; ++ell)
        {
            double lambda = wavelengthGrid->wavelength(ell);
            double dlambda = wavelengthGrid->effectiveWidth(ell);
            double opacity = ms->opacityAbs(lambda, m, MaterialMix::MaterialType::Dust);
            double product = opacity * Jv[ell] * dlambda;
            sumtop += product * lambda;
            sumbot += product;
        }
        if (sumbot > 0.)
            return sumtop / sumbot;
        else
            return 0.;
    }
}

////////////////////////////////////////////////////////////////////

vector<int> TemperatureWavelengthCellLibrary::mapping(const Array& bv) const
{
    // get the radiation field wavelength grid and the medium system
    auto wavelengthGrid = find<Configuration>()->radiationFieldWLG();
    auto ms = find<MediumSystem>();
    int numCells = ms->numCells();

    // calculate the indicative temperature and wavelength for all spatial cells
    // and track their minimum and maximum values
    double Tmin = DBL_MAX;
    double Tmax = 0.0;
    double lambdamin = DBL_MAX;
    double lambdamax = 0.0;
    Array Tv(numCells);
    Array lambdav(numCells);
    for (int m = 0; m != numCells; ++m)
    {
        // ignore cells that won't be used by the caller
        if (bv[m])
        {
            double T = ms->indicativeDustTemperature(m);
            double lambda = indicativeDustWavelength(m, ms, wavelengthGrid);

            // ignore cells with meaningless property values
            if (T > 0. && lambda > 0.)
            {
                Tv[m] = T;
                lambdav[m] = lambda;

                Tmin = min(Tmin, T);
                Tmax = max(Tmax, T);
                lambdamin = min(lambdamin, lambda);
                lambdamax = max(lambdamax, lambda);
            }
        }
    }

    // log the property ranges
    auto log = find<Log>();
    auto units = find<Units>();
    log->info("  Indicative temperatures vary"
              " from T = "
              + StringUtils::toString(units->otemperature(Tmin), 'g', 5) + " " + units->utemperature()
              + " to T = " + StringUtils::toString(units->otemperature(Tmax), 'g', 5) + " " + units->utemperature());
    log->info("  Indicative wavelengths vary"
              " from λ = "
              + StringUtils::toString(units->owavelength(lambdamin), 'g', 5) + " " + units->uwavelength()
              + " to λ = " + StringUtils::toString(units->owavelength(lambdamax), 'g', 5) + " " + units->uwavelength());

    // determine for every dust cell m the corresponding library entry n
    double dT = (Tmax - Tmin) / _numTemperatures;
    double loglambdamin = log10(lambdamin);
    double loglambdamax = log10(lambdamax);
    double dloglambda = (loglambdamax - loglambdamin) / _numWavelengths;
    vector<int> nv(numCells);
    for (int m = 0; m != numCells; ++m)
    {
        if (Tv[m] > 0. && lambdav[m] > 0.)
        {
            double T = Tv[m];
            int i = max(0, min(_numTemperatures - 1, static_cast<int>((T - Tmin) / dT)));

            double lambda = lambdav[m];
            double loglambda = log10(lambda);
            int j = max(0, min(_numWavelengths - 1, static_cast<int>((loglambda - loglambdamin) / dloglambda)));

            nv[m] = i + _numTemperatures * j;
        }
        else
        {
            nv[m] = -1;
        }
    }

    return nv;
}

////////////////////////////////////////////////////////////////////
