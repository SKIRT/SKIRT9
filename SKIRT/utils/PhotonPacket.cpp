/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PhotonPacket.hpp"
#include "AngularDistributionInterface.hpp"
#include "Constants.hpp"
#include "PolarizationProfileInterface.hpp"
#include "VelocityInterface.hpp"

////////////////////////////////////////////////////////////////////

PhotonPacket::PhotonPacket()
{
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launch(size_t historyIndex, double lambda, double L, Position bfr, Direction bfk,
                          VelocityInterface* bvi,
                          AngularDistributionInterface* adi,
                          PolarizationProfileInterface* ppi)
{
    _lambda = lambda;
    _W = L * lambda;
    _lambda0 = lambda;
    _bvi = bvi;
    _adi = adi;
    _ppi = ppi;
    _compIndex = 0;
    _historyIndex = historyIndex;
    _nscatt = 0;
    setPosition(bfr);
    setDirection(bfk);
    if (bvi) _lambda = shiftedEmissionWavelength(lambda, bfk, bvi->velocity());
    if (ppi) setPolarized(ppi->polarizationForDirection(bfk));
    else setUnpolarized();
    _hasObservedOpticalDepth = false;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::setPrimaryOrigin(int sourceCompIndex)
{
    _compIndex = sourceCompIndex + 1;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::setSecondaryOrigin(int mediumCompIndex)
{
    _compIndex = - (mediumCompIndex + 1);
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launchEmissionPeelOff(const PhotonPacket* pp, Direction bfk)
{
    _lambda = pp->_lambda;
    _W = pp->_W;
    _lambda0 = pp->_lambda0;
    _compIndex = pp->_compIndex;
    _historyIndex = pp->_historyIndex;
    _nscatt = 0;
    setPosition(pp->position());
    setDirection(bfk);
    if (pp->_bvi) _lambda = shiftedEmissionWavelength(_lambda0, bfk, pp->_bvi->velocity());
    if (pp->_adi) applyBias(pp->_adi->probabilityForDirection(bfk));
    if (pp->_ppi) setPolarized(pp->_ppi->polarizationForDirection(bfk));
    else setUnpolarized();
    _hasObservedOpticalDepth = false;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launchScatteringPeelOff(const PhotonPacket* pp, Direction bfk, Vec bfv, double w)
{
    if (bfv.isNull()) _lambda = pp->_lambda;
    else _lambda = shiftedEmissionWavelength(shiftedReceptionWavelength(pp->_lambda, pp->direction(), bfv), bfk, bfv);
    _W = pp->_W * w;
    _lambda0 = pp->_lambda0;
    _compIndex = pp->_compIndex;
    _historyIndex = pp->_historyIndex;
    _nscatt = pp->_nscatt + 1;
    setPosition(pp->position());
    setDirection(bfk);
    setUnpolarized();
    _hasObservedOpticalDepth = false;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::propagate(double s)
{
    propagatePosition(s);
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::scatter(Direction bfk, Vec bfv)
{
    if (!bfv.isNull())
        _lambda = shiftedEmissionWavelength(shiftedReceptionWavelength(_lambda, direction(), bfv), bfk, bfv);
    _nscatt++;
    setDirection(bfk);
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::applyBias(double w)
{
    _W *= w;
}

////////////////////////////////////////////////////////////////////

double PhotonPacket::shiftedEmissionWavelength(double sourceWavelength, Direction photonDirection, Vec sourceVelocity)
{
    return sourceWavelength * (1 - Vec::dot(photonDirection,sourceVelocity)/Constants::c());
}

////////////////////////////////////////////////////////////////////

double PhotonPacket::shiftedReceptionWavelength(double photonWavelength, Direction photonDirection, Vec receiverVelocity)
{
    return photonWavelength / (1 - Vec::dot(photonDirection,receiverVelocity)/Constants::c());
}

////////////////////////////////////////////////////////////////////

double PhotonPacket::perceivedWavelength(Vec receiverVelocity) const
{
    return shiftedReceptionWavelength(_lambda, direction(), receiverVelocity);
}

////////////////////////////////////////////////////////////////////
