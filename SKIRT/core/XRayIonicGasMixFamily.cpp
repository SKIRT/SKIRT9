/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayIonicGasMixFamily.hpp"
#include "Atoms.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"

XRayIonicGasMixFamily::~XRayIonicGasMixFamily()
{
    for (XRayIonicGasMix* mix : _mixes) delete mix;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMixFamily::setupSelfBefore()
{
    MaterialMixFamily::setupSelfBefore();

    setup();
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> XRayIonicGasMixFamily::parameterInfo() const
{
    vector<SnapshotParameter> descriptors;

    for (string ionName : _ionNames) descriptors.push_back(SnapshotParameter::custom(ionName));

    return descriptors;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* XRayIonicGasMixFamily::mix(double /*Z*/, double T, const Array& parameters)
{
    // convert Array to vector
    vector<double> abundances(begin(parameters), end(parameters));

    // look for duplicates
    for (const XRayIonicGasMix* mix : _mixes)
    {
        if (mix->abundances() == abundances) return mix;
    }

    XRayIonicGasMix* mix = new XRayIonicGasMix(this, _ions, abundances, T, _boundElectrons, _resonantScattering, true);
    _mixes.push_back(mix);

    return mix;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMixFamily::setup()
{
    if (_setupDone) return;
    _setupDone = true;

    // read ions
    // _numIons = StringUtils::split(ions(), ",").size();
    // parse all required ions from the ions property
    string ionString = StringUtils::squeeze(ions());
    if (ionString.empty()) throw FATALERROR("No ions specified");
    for (string ion : StringUtils::split(ionString, ","))
    {
        // verify ion is in correct format
        Atoms::parseIon(ion);
        // store isolated ion name
        _ionNames.push_back(ion);
    }

    // convert enum
    switch (electronScattering())
    {
        case ElectronScattering::None: _boundElectrons = XRayIonicGasMix::ElectronScattering::None; break;
        case ElectronScattering::Free: _boundElectrons = XRayIonicGasMix::ElectronScattering::Free; break;
        case ElectronScattering::FreeWithPolarization:
            _boundElectrons = XRayIonicGasMix::ElectronScattering::FreeWithPolarization;
            break;
        case ElectronScattering::Good: _boundElectrons = XRayIonicGasMix::ElectronScattering::Good; break;
        case ElectronScattering::Exact: _boundElectrons = XRayIonicGasMix::ElectronScattering::Exact; break;
    }
}

////////////////////////////////////////////////////////////////////

const MaterialMix* XRayIonicGasMixFamily::mix()
{
    if (_mixes.size() > 0) return _mixes[0];

    // create a default mix in case this is called before the setupSelfBefore()
    if (!_defaultMix)
    {
        setup();

        vector<double> abundances(_ionNames.size(), 0.);
        _defaultMix = new XRayIonicGasMix(this, _ions, abundances, 0., _boundElectrons, _resonantScattering, true);
    }
    return _defaultMix;
}
