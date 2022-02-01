/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayAtomicGasMix.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // support the elements with atomic number up to 30
    // H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
    constexpr size_t numAtoms = 30;

    // default abundancies taken from Table 2 of Anders & Grevesse (1989), the default abundance table in Xspec
    constexpr std::initializer_list<double> defaultAbundancies = {
        1.00E+00, 9.77E-02, 1.45E-11, 1.41E-11, 3.98E-10, 3.63E-04, 1.12E-04, 8.51E-04, 3.63E-08, 1.23E-04,
        2.14E-06, 3.80E-05, 2.95E-06, 3.55E-05, 2.82E-07, 1.62E-05, 3.16E-07, 3.63E-06, 1.32E-07, 2.29E-06,
        1.26E-09, 9.77E-08, 1.00E-08, 4.68E-07, 2.45E-07, 4.68E-05, 8.32E-08, 1.78E-06, 1.62E-08, 3.98E-08};
    static_assert(numAtoms == defaultAbundancies.size(), "Incorrect number of default abundancies");
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    // verify the number of abundancies; if the list is empty, use our default list
    if (_abundancies.empty())
        _abundancies = defaultAbundancies;
    else if (_abundancies.size() != numAtoms)
        throw FATALERROR("The abundancies list must have exactly " + std::to_string(numAtoms) + " values");
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType XRayAtomicGasMix::materialType() const
{
    return MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool XRayAtomicGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayAtomicGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity()};
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::sectionAbs(double lambda) const
{
    (void)lambda;
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::sectionSca(double lambda) const
{
    (void)lambda;
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::sectionExt(double lambda) const
{
    (void)lambda;
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionAbs(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionSca(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionExt(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/, double& /*lambda*/,
                                         double /*w*/, Direction /*bfkobs*/, Direction /*bfky*/,
                                         const MaterialState* /*state*/, const PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::performScattering(double /*lambda*/, const MaterialState* /*state*/, PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::indicativeTemperature(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    return temperature();
}

////////////////////////////////////////////////////////////////////
