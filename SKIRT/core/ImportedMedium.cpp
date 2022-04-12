/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedMedium.hpp"
#include "Configuration.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedMedium::setupSelfAfter()
{
    Medium::setupSelfAfter();

    // create the snapshot with preconfigured mass or density column
    _snapshot = createAndOpenSnapshot();

    // get list of requested snapshot parameters, if any
    vector<SnapshotParameter> parameterInfo =
        _importVariableMixParams ? _materialMixFamily->parameterInfo() : _materialMix->parameterInfo();

    // define function to discover standard snapshot parameters
    using Identifier = SnapshotParameter::Identifier;
    auto hasParam = [&parameterInfo](Identifier identifier) {
        for (const auto& param : parameterInfo)
            if (param.identifier() == identifier) return true;
        return false;
    };

    // add optional columns if applicable; skip metallicity or temperature if also requested by the material mix
    if (_importMetallicity && !hasParam(Identifier::Metallicity)) _snapshot->importMetallicity();
    if (_importTemperature && !hasParam(Identifier::Temperature)) _snapshot->importTemperature();
    if (hasVelocity()) _snapshot->importVelocity();
    if (_importMagneticField) _snapshot->importMagneticField();

    // add snapshot parameter import columns if applicable
    if (!parameterInfo.empty()) _snapshot->importParameters(parameterInfo);

    // set the density policy
    if (mix()->isDust())
    {
        // for dust, use temperature cutoff and metallicity multiplier
        _snapshot->setMassDensityPolicy(_massFraction, _importTemperature ? _maxTemperature : 0., _importMetallicity);
    }
    else
    {
        // for gas and electrons, do not use temperature or metallicity
        _snapshot->setMassDensityPolicy(_massFraction, 0., false);
    }

    // read the data from file
    _snapshot->readAndClose();
}

//////////////////////////////////////////////////////////////////////

ImportedMedium::~ImportedMedium()
{
    delete _snapshot;
}

//////////////////////////////////////////////////////////////////////

int ImportedMedium::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

const MaterialMix* ImportedMedium::mix(Position bfr) const
{
    if (_importVariableMixParams)
    {
        Array params;
        _snapshot->parameters(bfr, params);
        return _materialMixFamily->mix(params);
    }
    else
        return _materialMix;
}

//////////////////////////////////////////////////////////////////////

const MaterialMix* ImportedMedium::mix() const
{
    if (_importVariableMixParams)
    {
        Array params(_materialMixFamily->parameterInfo().size());
        return _materialMixFamily->mix(params);
    }
    else
        return _materialMix;
}

//////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasVariableMix() const
{
    return _importVariableMixParams;
}

//////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasVelocity() const
{
    if (_importVelocity)
    {
        // refuse velocity for oligochromatic simulations
        // (this function is called from Configure so we cannot precompute this during setup)
        auto config = find<Configuration>();
        if (!config->oligochromatic()) return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////

Vec ImportedMedium::bulkVelocity(Position bfr) const
{
    return hasVelocity() ? _snapshot->velocity(bfr) : Vec();
}

//////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasMagneticField() const
{
    return _importMagneticField;
}

//////////////////////////////////////////////////////////////////////

Vec ImportedMedium::magneticField(Position bfr) const
{
    return _importMagneticField ? _snapshot->magneticField(bfr) : Vec();
}

////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasMetallicity() const
{
    return !_importVariableMixParams && _importMetallicity && !materialMix()->isDust();
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::metallicity(Position bfr) const
{
    if (hasMetallicity()) return _snapshot->metallicity(bfr);
    return 0.;
}

////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasTemperature() const
{
    return !_importVariableMixParams && _importTemperature && !materialMix()->isDust();
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::temperature(Position bfr) const
{
    if (hasTemperature()) return _snapshot->temperature(bfr);
    return 0.;
}

////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasParameters() const
{
    return !_importVariableMixParams && _snapshot->hasParameters();
}

////////////////////////////////////////////////////////////////////

void ImportedMedium::parameters(Position bfr, Array& params) const
{
    if (hasParameters())
        _snapshot->parameters(bfr, params);
    else
        params.resize(0);
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::numberDensity(Position bfr) const
{
    double result = _snapshot->density(bfr);
    if (!_snapshot->holdsNumber()) result /= mix(bfr)->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::number() const
{
    double result = _snapshot->mass();
    if (!_snapshot->holdsNumber()) result /= mix()->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::massDensity(Position bfr) const
{
    double result = _snapshot->density(bfr);
    if (_snapshot->holdsNumber()) result *= mix(bfr)->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::mass() const
{
    double result = _snapshot->mass();
    if (_snapshot->holdsNumber()) result *= mix()->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::opticalDepthX(double lambda) const
{
    double result = _snapshot->SigmaX() * mix()->sectionExt(lambda);
    if (!_snapshot->holdsNumber()) result /= mix()->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::opticalDepthY(double lambda) const
{
    double result = _snapshot->SigmaY() * mix()->sectionExt(lambda);
    if (!_snapshot->holdsNumber()) result /= mix()->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::opticalDepthZ(double lambda) const
{
    double result = _snapshot->SigmaZ() * mix()->sectionExt(lambda);
    if (!_snapshot->holdsNumber()) result /= mix()->mass();
    return result;
}

//////////////////////////////////////////////////////////////////////

Position ImportedMedium::generatePosition() const
{
    return _snapshot->generatePosition();
}

//////////////////////////////////////////////////////////////////////

int ImportedMedium::numSites() const
{
    return _snapshot->numEntities();
}

//////////////////////////////////////////////////////////////////////

Position ImportedMedium::sitePosition(int index) const
{
    return _snapshot->position(index);
}

//////////////////////////////////////////////////////////////////////
