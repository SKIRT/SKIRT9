/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumState.hpp"
#include "FatalError.hpp"
#include "ProcessManager.hpp"

//////////////////////////////////////////////////////////////////////

void MediumState::initConfiguration(int numCells, int numComps)
{
    _numCells = numCells;
    _numComps = numComps;

    _off_dens.resize(_numComps);
    _off_temp.resize(_numComps);
    _off_cust.resize(_numComps);
}

//////////////////////////////////////////////////////////////////////

void MediumState::initCommonStateVariables(const vector<StateVariable>& variables)
{
    for (const StateVariable& variable : variables)
    {
        switch (variable.identifier())
        {
            case StateVariable::Identifier::Volume: _off_volu = _nextOffset++; break;
            case StateVariable::Identifier::BulkVelocity:
                _off_velo = _nextOffset;
                _nextOffset += 3;
                break;
            case StateVariable::Identifier::MagneticField:
                _off_mfld = _nextOffset;
                _nextOffset += 3;
                break;
            case StateVariable::Identifier::NumberDensity:
            case StateVariable::Identifier::Temperature:
            case StateVariable::Identifier::Custom:
                throw FATALERROR("Requesting common state variable of unsupported type");
        }
    }
}

//////////////////////////////////////////////////////////////////////

void MediumState::initSpecificStateVariables(const vector<StateVariable>& variables)
{
    for (const StateVariable& variable : variables)
    {
        switch (variable.identifier())
        {
            case StateVariable::Identifier::NumberDensity: _off_dens[_nextComponent] = _nextOffset++; break;
            case StateVariable::Identifier::Temperature: _off_temp[_nextComponent] = _nextOffset++; break;
            case StateVariable::Identifier::Custom:
                if (variable.customIndex() == 0) _off_cust[_nextComponent] = _nextOffset;
                _nextOffset++;
                break;
            case StateVariable::Identifier::Volume:
            case StateVariable::Identifier::BulkVelocity:
            case StateVariable::Identifier::MagneticField:
                throw FATALERROR("Requesting specific state variable of unsupported type");
        }
    }
    _nextComponent++;
}

//////////////////////////////////////////////////////////////////////

size_t MediumState::initAllocate()
{
    if (_nextComponent != _numComps) throw FATALERROR("Failed to request state variables for all medium components");
    _numVars = _nextOffset;

    size_t numAlloc = static_cast<size_t>(_numVars) * static_cast<size_t>(_numCells);
    _data.resize(numAlloc);
    return numAlloc;
}

//////////////////////////////////////////////////////////////////////

void MediumState::initCommunicate()
{
    ProcessManager::sumToAll(_data);
}

//////////////////////////////////////////////////////////////////////

void MediumState::setVolume(double value, int m)
{
    _data[_numVars * m + _off_volu] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setBulkVelocity(Vec value, int m)
{
    int i = _numVars * m + _off_velo;
    _data[i] = value.x();
    _data[i + 1] = value.y();
    _data[i + 2] = value.z();
}

//////////////////////////////////////////////////////////////////////

void MediumState::setMagneticField(Vec value, int m)
{
    int i = _numVars * m + _off_mfld;
    _data[i] = value.x();
    _data[i + 1] = value.y();
    _data[i + 2] = value.z();
}

//////////////////////////////////////////////////////////////////////

void MediumState::setNumberDensity(double value, int m, int h)
{
    _data[_numVars * m + _off_dens[h]] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setTemperature(double value, int m, int h)
{
    _data[_numVars * m + _off_temp[h]] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setCustom(double value, int i, int m, int h)
{
    _data[_numVars * m + _off_cust[h] + i] = value;
}

//////////////////////////////////////////////////////////////////////
