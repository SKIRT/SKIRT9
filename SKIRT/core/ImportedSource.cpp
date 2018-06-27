/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSource.hpp"
#include "Constants.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "RedshiftInterface.hpp"
#include "Snapshot.hpp"
#include "WavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // An instance of this class offers the redshift interface for the bulk velocity specified in the constructor
    class BulkVelocity : public RedshiftInterface
    {
    private:
        Vec _bfv;
    public:
        BulkVelocity(Vec bfv) : _bfv(bfv) { }
        double redshiftForDirection(Direction bfk) const override { return -Vec::dot(_bfv, bfk) / Constants::c(); }
    };
}

//////////////////////////////////////////////////////////////////////

void ImportedSource::setupSelfBefore()
{
    Source::setupSelfBefore();

    // create the snapshot with preconfigured spatial columns
    _snapshot = createAndOpenSnapshot();

    // add optional columns if applicable
    if (_importVelocity) _snapshot->importVelocity();
    //_snapshot->importParameters(1);     // XXXX

    // read the data from file
    _snapshot->readAndClose();

    // construct a vector with the luminosity for each entity
    int M = _snapshot->numEntities();
    _Lv.resize(M);
    for (int m=0; m!=M; ++m)
    {
        _Lv[m] = 1. * Constants::Lsun();    // XXX
    }
    _L = _Lv.sum();
    _Lv /= _L;
}

//////////////////////////////////////////////////////////////////////

ImportedSource::~ImportedSource()
{
    delete _snapshot;
}

//////////////////////////////////////////////////////////////////////

int ImportedSource::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

double ImportedSource::luminosity() const
{
   return _L;
}

//////////////////////////////////////////////////////////////////////

double ImportedSource::specificLuminosity(double wavelength) const
{
    if (!interface<WavelengthRangeInterface>()->wavelengthRange().contains(wavelength)) return 0.;
  //  return _sed->specificLuminosity(wavelength) * luminosity();   XXX
    return 0;
}

//////////////////////////////////////////////////////////////////////

void ImportedSource::prepareForLaunch(double sourceBias, size_t firstIndex, size_t numIndices)
{
    int M = _snapshot->numEntities();

    // calculate the launch weight for each entity, normalized to unity
    _Wv = (1-sourceBias)*_Lv + sourceBias/M;

    // determine the first history index for each entity
    _Iv.resize(M+1);
    _Iv[0] = firstIndex;
    for (int m=0; m!=M; ++m)
    {
        // limit first index to last index to avoid run-over due to rounding errors
        _Iv[m+1] = min(firstIndex+numIndices, _Iv[m] + static_cast<size_t>(std::round(_Wv[m] * numIndices)));
    }
}

//////////////////////////////////////////////////////////////////////

void ImportedSource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // select the entity corresponding to this history index
    auto m = std::upper_bound(_Iv.cbegin(), _Iv.cend(), historyIndex) - _Iv.cbegin() - 1;

    // calculate the weight related to biased source selection
    double ws = _Lv[m] / _Wv[m];

    // get the SED for this entity
    // XXXX

    // generate a random wavelength from the SED and/or from the bias distribution
    double lambda, w;
    double xi = wavelengthBias();
    if (!xi)
    {
        // no biasing -- simply use the intrinsic spectral distribution
        lambda = _sed->generateWavelength();
        w = 1.;
    }
    else
    {
        // biasing -- use one or the other distribution
        if (random()->uniform() > xi) lambda = sed()->generateWavelength();
        else lambda = wavelengthBiasDistribution()->generateWavelength();

        // calculate the compensating weight factor
        double s = sed()->specificLuminosity(lambda);
        if (!s)
        {
            // if the wavelength can't occur in the intrinsic distribution,
            // the weight factor is zero regardless of the probability in the bias distribution
            // (handling this separately also avoids NaNs in the pathological case s=d=0)
            w = 0.;
        }
        else
        {
            // regular composite bias weight
            double b = wavelengthBiasDistribution()->probability(lambda);
            w = s / ((1-xi)*s + xi*b);
        }
    }

    // generate a random position for this entity
    Position bfr = _snapshot->generatePosition(m);

    // provide a redshift interface for the appropriate velocity, if enabled
    // XXXX

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, L*w*ws, bfr, random()->direction());
}

//////////////////////////////////////////////////////////////////////
