/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumState.hpp"
#include "FatalError.hpp"
#include "ProcessManager.hpp"

//////////////////////////////////////////////////////////////////////

void MediumState::initConfiguration(int numCells, int numMedia)
{
    _numCells = numCells;
    _numMedia = numMedia;

    _off_dens.resize(_numMedia);
    _off_temp.resize(_numMedia);
    _off_cust.resize(_numMedia);
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
    if (_nextComponent != _numMedia) throw FATALERROR("Failed to request state variables for all medium components");
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

int MediumState::synchronize(const std::vector<uint8_t>& cellFlags)
{
    int numUpdated = 0;

    if (ProcessManager::isMultiProc())
    {
        auto producer = [this, &cellFlags, &numUpdated](vector<double>& data) {
            for (int m = 0; m != _numCells; ++m)
            {
                if (cellFlags[m])
                {
                    data.push_back(m);
                    double* first = &_data[_numVars * m];
                    data.insert(data.end(), first, first + _numVars);
                    numUpdated++;
                }
            }
        };
        auto consumer = [this, &numUpdated](const vector<double>& data) {
            for (auto in = data.begin(); in != data.end(); in += (1 + _numVars))
            {
                int m = *in;
                std::copy(in + 1, in + 1 + _numVars, &_data[_numVars * m]);
                numUpdated++;
            }
        };
        ProcessManager::broadcastAllToAll(producer, consumer);
    }
    else
    {
        for (int m = 0; m != _numCells; ++m)
            if (cellFlags[m]) numUpdated++;
    }
    return numUpdated;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setVolume(int m, double value)
{
    _data[_numVars * m + _off_volu] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setBulkVelocity(int m, Vec value)
{
    int i = _numVars * m + _off_velo;
    _data[i] = value.x();
    _data[i + 1] = value.y();
    _data[i + 2] = value.z();
}

//////////////////////////////////////////////////////////////////////

void MediumState::setMagneticField(int m, Vec value)
{
    int i = _numVars * m + _off_mfld;
    _data[i] = value.x();
    _data[i + 1] = value.y();
    _data[i + 2] = value.z();
}

//////////////////////////////////////////////////////////////////////

void MediumState::setNumberDensity(int m, int h, double value)
{
    _data[_numVars * m + _off_dens[h]] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setTemperature(int m, int h, double value)
{
    _data[_numVars * m + _off_temp[h]] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setCustom(int m, int h, int i, double value)
{
    _data[_numVars * m + _off_cust[h] + i] = value;
}

//////////////////////////////////////////////////////////////////////
