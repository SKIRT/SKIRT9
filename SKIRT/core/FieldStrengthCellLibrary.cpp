/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FieldStrengthCellLibrary.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

int FieldStrengthCellLibrary::numEntries() const
{
    return _numFieldStrengths;
}

////////////////////////////////////////////////////////////////////

vector<int> FieldStrengthCellLibrary::mapping(const Array& bv) const
{
    // get the radiation field wavelength grid and the medium system
    auto wavelengthGrid = find<Configuration>()->radiationFieldWLG();
    auto ms = find<MediumSystem>();
    int numCells = ms->numCells();

    // the local radiation field in the Milky Way (Mathis et al. 1983) integrated over all wavelengths
    double JtotMW = 1.7623e-06;

    // calculate the field strengths for all spatial cells and track the minimum and maximum values
    double Umin = DBL_MAX;
    double Umax = 0.0;
    Array Uv(numCells);
    for (int m = 0; m != numCells; ++m)
    {
        // ignore cells that won't be used by the caller
        if (bv[m])
        {
            double U = (ms->meanIntensity(m) * wavelengthGrid->dlambdav()).sum() / JtotMW;
            // ignore cells with extremely small radiation fields (compared to the average in the Milky Way)
            // to avoid wasting library grid points on fields that won't change simulation results anyway
            if (U > 1e-6)
            {
                Uv[m] = U;
                Umin = min(Umin, U);
                Umax = max(Umax, U);
            }
        }
    }

    // log the field strength range
    find<Log>()->info("  Radiation field strengths vary from U = " + StringUtils::toString(Umin, 'e', 4)
                      + " to U = " + StringUtils::toString(Umax, 'e', 4));

    // determine for every dust cell m the corresponding library entry n
    double logUmin = log10(Umin);
    double logUmax = log10(Umax);
    double dlogU = (logUmax - logUmin) / _numFieldStrengths;
    vector<int> nv(numCells);
    for (int m = 0; m != numCells; ++m)
    {
        double U = Uv[m];
        if (U > 0.)
        {
            double logU = log10(U);
            nv[m] = max(0, min(_numFieldStrengths - 1, static_cast<int>((logU - logUmin) / dlogU)));
        }
        else
        {
            nv[m] = -1;
        }
    }

    return nv;
}

////////////////////////////////////////////////////////////////////
