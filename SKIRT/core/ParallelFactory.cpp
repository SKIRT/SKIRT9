/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParallelFactory.hpp"
#include "FatalError.hpp"
#include "ProcessManager.hpp"

////////////////////////////////////////////////////////////////////

ParallelFactory::ParallelFactory()
{
}

////////////////////////////////////////////////////////////////////

ParallelFactory::ParallelFactory(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void ParallelFactory::setMaxThreadCount(int value)
{
    _maxThreadCount = std::max(1, value);
}

////////////////////////////////////////////////////////////////////

int ParallelFactory::maxThreadCount() const
{
    return _maxThreadCount;
}

////////////////////////////////////////////////////////////////////

int ParallelFactory::defaultThreadCount()
{
    int count = std::thread::hardware_concurrency();
    if (count < 1) return 1;    // the number of threads could not be determined
    if (count > 24) return 24;  // additional threads in single process do not increase performance
    return count;
}

////////////////////////////////////////////////////////////////////

Parallel* ParallelFactory::parallel(TaskMode mode, int maxThreadCount)
{
    // Verify that we're being called from our parent thread
    if (std::this_thread::get_id() != _parentThread)
        throw FATALERROR("Parallel not spawned from thread that constructed the factory");

    // Determine the number of threads
    //  - limited by both the factory maximum and the maximum specified here as an argument
    //  - reduced to one for Duplicated mode in multiprocessing environment
    int numThreads = maxThreadCount>0 ? std::min(maxThreadCount, _maxThreadCount) : _maxThreadCount;
    if (mode == TaskMode::Duplicated && ProcessManager::isMultiProc()) numThreads = 1;

    // Determine the Parallel subclass type (see class documentation for details)
    ParallelType type = numThreads == 1 ? ParallelType::Serial : ParallelType::MultiThread;
    if (ProcessManager::isMultiProc())
    {
        if (mode == TaskMode::Distributed)
            type = numThreads == 1 ? ParallelType::MultiProcess : ParallelType::Hybrid;
        else if (mode == TaskMode::RootOnly && !ProcessManager::isRoot())
            type = ParallelType::Null;
        // for the other cases, the type is already set correctly by the very first assignement
    }

    // Get or create a child of that type and with that number of threads
    auto& child = _children[std::make_pair(type, numThreads)];
    if (!child) child.reset( new Parallel(numThreads) );
    return child.get();
}

////////////////////////////////////////////////////////////////////
