/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PhotonPacket.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

PhotonPacket::PhotonPacket()
{
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launch(size_t historyIndex, double lambda, double L, Position bfr, Direction bfk)
{
    _lambda = lambda;
    _nu = Constants::c() / lambda;
    _W = L / _nu;
    _nscatt = 0;
    _compIndex = 0;
    _historyIndex = historyIndex;
    _bfr = bfr;
    _bfk = bfk;
    setUnpolarized();
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

void PhotonPacket::launchEmissionPeelOff(const PhotonPacket* pp, Direction bfk, double w)
{
    _lambda = pp->_lambda;
    _nu = pp->_nu;
    _W = pp->_W * w;
    _nscatt = 0;
    _compIndex = pp->_compIndex;
    _historyIndex = pp->_historyIndex;
    _bfr = pp->_bfr;
    _bfk = bfk;
    setUnpolarized();
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launchScatteringPeelOff(const PhotonPacket* pp, Direction bfk, double w)
{
    _lambda = pp->_lambda;
    _nu = pp->_nu;
    _W = pp->_W * w;
    _nscatt = pp->_nscatt + 1;
    _compIndex = pp->_compIndex;
    _historyIndex = pp->_historyIndex;
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

void PhotonPacket::scatter(Direction bfk)
{
    _nscatt++;
    _bfk = bfk;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::applyBias(double w)
{
    _W *= w;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::applyRedshift(double z)
{
    _lambda *= 1+z;
    _nu /= 1+z;
}

////////////////////////////////////////////////////////////////////
