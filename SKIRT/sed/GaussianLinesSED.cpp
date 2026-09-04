/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GaussianLinesSED.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // intrinsic half-range of a single line in dispersion units (see LyaGaussianSED for a
    // discussion of this limit: the Gaussian value at 9 dispersion units is negligible)
    const double intrinsicRange = 9;

    // number of wavelength points per dispersion unit in the private grid for a single line
    const int numWavelengthsPerDispersionUnit = 100;

    // relative wavelength difference below which two grid points are considered coincident
    const double mergeTolerance = 1e-6;

    // Gaussian centered on 0 with dispersion of 1, evaluated at x
    double unitGaussian(double x)
    {
        return (0.5 * M_SQRT1_2 * M_2_SQRTPI) * exp(-0.5 * x * x);
    }
}

////////////////////////////////////////////////////////////////////

void GaussianLinesSED::getWavelengthsAndLuminosities(Array& lambdav, Array& pv) const
{
    // ask the subclass for the line properties
    Array inlambdav, indispersionv, inLv;
    getLineProperties(inlambdav, indispersionv, inLv);

    // verify the number of lines and the value of their properties
    size_t numLines = inlambdav.size();
    if (!numLines || indispersionv.size() != numLines || inLv.size() != numLines)
        throw FATALERROR("Gaussian lines SED must have a nonzero and equal number of wavelengths, "
                         "dispersions and luminosities");
    for (size_t i = 0; i != numLines; ++i)
    {
        if (indispersionv[i] <= 0) throw FATALERROR("Gaussian line velocity dispersion must be positive");
        if (inLv[i] < 0) throw FATALERROR("Gaussian line luminosity cannot be negative");
    }

    // convert to a list of lines with the dispersion expressed in wavelength units, and sort by
    // center wavelength so that the lines can be specified in any order
    struct Line
    {
        double lambda;  // center wavelength
        double sigma;   // dispersion in wavelength units
        double L;       // relative luminosity
    };
    vector<Line> linev(numLines);
    for (size_t i = 0; i != numLines; ++i)
        linev[i] = Line{inlambdav[i], inlambdav[i] * indispersionv[i] / Constants::c(), inLv[i]};
    std::sort(linev.begin(), linev.end(), [](const Line& l1, const Line& l2) { return l1.lambda < l2.lambda; });

    // verify that at least one line falls inside the source wavelength range, so that normalization
    // (performed by the TabulatedSED base class) does not silently fail
    auto sourceRange = find<Configuration>()->sourceWavelengthRange();
    bool hasLineInRange = false;
    for (const auto& line : linev)
        if (sourceRange.containsFuzzy(line.lambda)) hasLineInRange = true;
    if (!hasLineInRange)
        throw FATALERROR("Gaussian lines SED must have at least one line in the source wavelength range");

    // for each line, build a private linear grid limited to its own intrinsic range,
    // and pool all grid points together
    vector<double> pooledv;
    for (const auto& line : linev)
    {
        double lo = line.lambda - intrinsicRange * line.sigma;
        double hi = line.lambda + intrinsicRange * line.sigma;
        int n = static_cast<int>(2 * intrinsicRange * numWavelengthsPerDispersionUnit);
        Array grid;
        NR::buildLinearGrid(grid, lo, hi, n);
        for (double x : grid) pooledv.push_back(x);
    }

    // sort the pooled points and merge points that end up extremely close together, so that no
    // two final grid points are (nearly) coincident
    NR::sort(pooledv);
    vector<double> mergedv;
    mergedv.reserve(pooledv.size());
    for (double x : pooledv)
        if (mergedv.empty() || x > mergedv.back() * (1. + mergeTolerance)) mergedv.push_back(x);

    // tabulate the sum of all Gaussian line profiles on the merged grid
    NR::assign(lambdav, mergedv);
    size_t numPoints = lambdav.size();
    pv.resize(numPoints);
    for (size_t k = 0; k != numPoints; ++k)
    {
        double lambda = lambdav[k];
        double sum = 0.;
        for (const auto& line : linev) sum += line.L * unitGaussian((lambda - line.lambda) / line.sigma) / line.sigma;
        pv[k] = sum;
    }
}

////////////////////////////////////////////////////////////////////
