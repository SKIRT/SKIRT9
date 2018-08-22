/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PhotonPacket.hpp"
#include "Constants.hpp"
#include "AngularDistributionInterface.hpp"
#include "PolarizationProfileInterface.hpp"
#include "BulkVelocityInterface.hpp"

////////////////////////////////////////////////////////////////////

PhotonPacket::PhotonPacket()
{
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launch(size_t historyIndex, double lambda, double L, Position bfr, Direction bfk,
                          BulkVelocityInterface* bvi,
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
    _bfr = bfr;
    _bfk = bfk;
    if (bvi) _lambda = shiftedEmissionWavelength(lambda, bfk, bvi->bulkVelocity());
    if (ppi) setPolarized(ppi->polarizationForDirection(bfk));
    else setUnpolarized();
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
    _bfr = pp->_bfr;
    _bfk = bfk;
    if (pp->_bvi) _lambda = shiftedEmissionWavelength(_lambda0, bfk, _bvi->bulkVelocity());
    if (pp->_adi) applyBias(pp->_adi->probabilityForDirection(bfk));
    if (pp->_ppi) setPolarized(pp->_ppi->polarizationForDirection(bfk));
    else setUnpolarized();
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launchScatteringPeelOff(const PhotonPacket* pp, Direction bfk, Vec bfv, double w)
{
    if (bfv.isNull()) _lambda = pp->_lambda;
    else _lambda = shiftedEmissionWavelength(shiftedReceptionWavelength(pp->_lambda, pp->_bfk, bfv), bfk, bfv);
    _W = pp->_W * w;
    _lambda0 = pp->_lambda0;
    _compIndex = pp->_compIndex;
    _historyIndex = pp->_historyIndex;
    _nscatt = pp->_nscatt + 1;
    _bfr = pp->_bfr;
    _bfk = bfk;
    setUnpolarized();
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::propagate(double s)
{
    _bfr += s*_bfk;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::scatter(Direction bfk, Vec bfv)
{
    if (!bfv.isNull()) _lambda = shiftedEmissionWavelength(shiftedReceptionWavelength(_lambda, _bfk, bfv), bfk, bfv);
    _nscatt++;
    _bfk = bfk;
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
    return shiftedReceptionWavelength(_lambda, _bfk, receiverVelocity);
}

////////////////////////////////////////////////////////////////////
