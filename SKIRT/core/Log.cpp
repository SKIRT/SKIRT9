/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Log.hpp"
#include "ProcessCommunicator.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <chrono>

////////////////////////////////////////////////////////////////////

Log::Log()
{
}

////////////////////////////////////////////////////////////////////

void Log::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // get a pointer to the ProcessCommunicator and remember our rank
    auto comm = find<ProcessCommunicator>();
    if (comm->isMultiProc()) setRank(comm->rank());
}

////////////////////////////////////////////////////////////////////

void Log::setRank(int rank)
{
    if (_link) _link->setRank(rank);

    string rankName = std::to_string(rank);
    if (rankName.length() < 3) rankName.insert(0, 3-rankName.length(), '0');
    _procNameShort = "P" + rankName;
    _procNameLong = "[" + _procNameShort + "] ";
}

////////////////////////////////////////////////////////////////////

void Log::setLowestLevel(Log::Level level)
{
    _lowestLevel = level;
}

////////////////////////////////////////////////////////////////////

Log::Level Log::lowestLevel() const
{
    return _lowestLevel;
}

////////////////////////////////////////////////////////////////////

void Log::setLinkedLog(Log* log)
{
    destroyChild(_link);
    _link = log;
    addChild(_link) ;
}

////////////////////////////////////////////////////////////////////

Log* Log::linkedLog() const
{
    return _link;
}

////////////////////////////////////////////////////////////////////

void Log::setVerbose(bool value)
{
    _verbose = value;
    if (_link) _link->setVerbose(value);
}

////////////////////////////////////////////////////////////////////

bool Log::verbose() const
{
    return _verbose;
}

////////////////////////////////////////////////////////////////////

void Log::setMemoryLogging(bool value)
{
    _logmemory = value;
    if (_link) _link->setMemoryLogging(value);
}

////////////////////////////////////////////////////////////////////

bool Log::memoryLogging() const
{
    return _logmemory;
}

////////////////////////////////////////////////////////////////////

void Log::info(string message)
{
    // Pass the message to the linked log
    if (_link) _link->info(message);

    // Obtain a string denoting the amount of used memory, if requested
    string memory = _logmemory ? "(" + StringUtils::toMemSizeString(System::currentMemoryUsage()) + ") " : "";

    // Output the message
    if (verbose())
    {
        if (Level::Info >= _lowestLevel) output(_procNameLong + memory + message, Level::Info);
    }
    else if (ProcessManager::isRoot())
    {
        if (Level::Info >= _lowestLevel) output(memory + message, Level::Info);
    }
}

////////////////////////////////////////////////////////////////////

void Log::infoSetElapsed(int seconds)
{
    using namespace std::chrono;

    _interval = static_cast<uint64_t>(seconds) * steady_clock::period::den / steady_clock::period::num;
    _started = steady_clock::now().time_since_epoch().count();
}

////////////////////////////////////////////////////////////////////

void Log::infoIfElapsed(string message)
{
    using namespace std::chrono;

    // get a local copy of the shared start time
    uint64_t localStarted = _started;
    while (true)
    {
        // if the interval has not yet expired, return without issuing the message
        uint64_t now = steady_clock::now().time_since_epoch().count();
        if (now < localStarted + _interval) return;

        // the interval has expired; attempt to reset the start time to the current time
        if (_started.compare_exchange_weak(localStarted, now))
        {
            // if no other thread changed the start time, issue the message and return
            info(message);
            return;
        }

        // another thread intervened, or there was a spurious failure, so try again
    }
}

////////////////////////////////////////////////////////////////////

void Log::warning(string message)
{
    // Pass the message to the linked log
    if (_link) _link->warning(message);

    // Obtain a string denoting the amount of used memory, if requested
    string memory = _logmemory ? "(" + StringUtils::toMemSizeString(System::currentMemoryUsage()) + ") " : "";

    // Output the message
    if (Level::Warning >= _lowestLevel) output(_procNameLong + memory + message, Level::Warning);
}

////////////////////////////////////////////////////////////////////

void Log::success(string message)
{
    // Pass the message to the linked log
    if (_link) _link->success(message);

    // Obtain a string denoting the amount of used memory, if requested
    string memory = _logmemory ? "(" + StringUtils::toMemSizeString(System::currentMemoryUsage()) + ") " : "";

    // Output the message
    if (verbose())
    {
        if (Level::Success >= _lowestLevel) output(_procNameLong + memory + message, Level::Success);
    }
    else if (ProcessManager::isRoot())
    {
        if (Level::Success >= _lowestLevel) output(memory + message, Level::Success);
    }
}

////////////////////////////////////////////////////////////////////

void Log::error(string message)
{
    // Pass the message to the linked log
    if (_link) _link->error(message);

    // Output the message
    if (Level::Error >= _lowestLevel) output(_procNameLong + "*** Error: " + message, Level::Error);
}

////////////////////////////////////////////////////////////////////

string Log::processName()
{
    return _procNameShort;
}

////////////////////////////////////////////////////////////////////
