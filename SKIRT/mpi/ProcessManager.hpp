/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROCESSMANAGER_HPP
#define PROCESSMANAGER_HPP

#include "Array.hpp"
#include <functional>

////////////////////////////////////////////////////////////////////

/** The ProcessManager class is an interface to the MPI library, and is the only place in the SKIRT
    code where explicit calls to this library are allowed. A single instance of the ProcessManager
    class must be constructed just after program startup, and before the command-line arguments are
    used by other parts of the program. The latter is important because the MPI initialization
    functions called by the ProcessManager constructor may adjust the command-line arguments to
    remove MPI specific options. The program should let the main() function run to normal
    completion and return an exit code rather than calling the std::exit() function. In case of
    abnormal termination, the program should call the ProcessManager::abort() function as opposed
    to the std::abort() function.

    The functions of this class should be called only from the main thread of the program.
    Violation of this rule causes undefined behavior that may differ between MPI implementations.
    */
class ProcessManager final
{
    //======== Construction - Destruction  ===========

private:
    /** This static function initializes the MPI library used for remote communication, if present.
        The function is called from the constructor of this class. */
    static void initialize(int* argc, char*** argv);

    /** This static function finalizes the MPI library used for remote communication, if present.
        The function is called from the destructor of this class. */
    static void finalize();

public:
    /** The constructor initializes the MPI library used for remote communication, if present. It
        is passed a reference to the command line arguments to provide the MPI library with the
        opportunity to use and/or adjust them (e.g. to remove any MPI related arguments). */
    ProcessManager(int* argc, char*** argv) { initialize(argc, argv); }

    /** The destructor finalizes the MPI library used for remote communication, if present. */
    ~ProcessManager() { finalize(); }

    /** The copy constructor is deleted because instances of this class should never be copied or
        moved. */
    ProcessManager(const ProcessManager&) = delete;

    /** The assignment operator is deleted because instances of this class should never be copied
        or moved. */
    ProcessManager& operator=(const ProcessManager&) = delete;

    /** This function should be called when a process experiences a fatal error (after the error
        was reported to the user). If there are two or more processes in the current run-time
        environment, this function aborts the complete MPI process group, including the calling
        process. If the MPI library is not present, or the program was invoked without MPI, or
        there is only one process, this function does nothing. */
    static void abort(int exitcode);

    //======== Logging  ===========

    /** This function sets the specified callback function as the usage logger, which is called
        with diagnostic messages bracketing each MPI operation. Initially, the logger is cleared.
        */
    static void setLogger(std::function<void(string)> logger);

    /** This function clears the usage logger so that it is no longer called. */
    static void clearLogger();

    //======== Environment info  ===========

    /** This function returns the number of processes in the current run-time environment. If the
        MPI library is not present, or the program was invoked without MPI, the function returns
        one. */
    static int size() { return _size; }

    /** This function returns true if the current run-time environment has two or more processes.
        If the MPI library is not present, or the program was invoked without MPI, the function
        returns false. */
    static bool isMultiProc() { return _size > 1; }

    /** This function returns the rank of the calling process in the current run-time environment,
        i.e. an integer in a range from zero to the number of processes minus one. If the MPI
        library is not present, or the program was invoked without MPI, the function always returns
        zero. */
    static int rank() { return _rank; }

    /** This function returns true if the calling process is considered to be the root process,
        i.e. its rank is zero. If the MPI library is not present, or the program was invoked
        without MPI, the function always returns true. */
    static bool isRoot() { return _rank == 0; }

    //======== Master-slave communication  ===========

    /** This function is part of the mechanism for dynamically allocating chunks of parallel
        tasks across multiple processes. It requests the next available chunk from the root process
        and waits for a response. When successful, the function places a chunk index range in its
        arguments and returns true. If no more chunks are available, the function returns false
        (and the output arguments are both set to zero). If there is only one process, or if the
        function is invoked from the root process, a fatal error is thrown. */
    static bool requestChunk(size_t& firstIndex, size_t& numIndices);

    /** This function is part of the mechanism for dynamically allocating chunks of parallel
        tasks across multiple processes. It waits for a chunk request from any of the processes in
        the MPI group and returns the rank of the requesting process. If there is only one process,
        or if the function is invoked from any process other than the root process, a fatal error
        is thrown. */
    static int waitForChunkRequest();

    /** This function is part of the mechanism for dynamically allocating chunks of parallel
        tasks across multiple processes. It communicates the index range of the next available
        chunk request to the process in the MPI group with the specified rank. If there is only one
        process, or if the function is invoked from any process other than the root process, a
        fatal error is thrown. */
    static void serveChunkRequest(int rank, size_t firstIndex, size_t numIndices);

    //======== Collective Communication  ===========

    /** This function causes the calling process to block until all other processes have invoked it
        as well. If there is only one process, the function does nothing. */
    static void wait();

    /** This function adds the floating point values of an array element-wise across the different
        processes. The resulting sums are then stored in the same Array passed to this function on
        each individual process. All processes must call this function for the communication to
        proceed. If there is only one process, or if the array has zero size, the function does
        nothing. */
    static void sumToAll(Array& arr);

    /** This function adds the floating point values of an array element-wise across the different
        processes. The resulting sums are then stored in the same Array passed to this function on
        the root process. The arrays on the other processes are left untouched. All processes must
        call this function for the communication to proceed. If there is only one process, or if
        the array has zero size, the function does nothing.

        If the \em wait argument is true, this function causes the processes to block until all
        other processes have invoked it as well. This is important in case the sumToRoot() call may
        be followed by a sequence of master-slave chunk requests without an intervening call to
        wait(). Indeed, because nonroot processes do not receive any data in this function, if the
        MPI implementation buffers the involved send operations, the sumToRoot() function may
        return in nonroot processes before the root process calls it. This causes problems in case
        the sumToRoot() call is followed by a sequence of master-slave chunk requests, because the
        nonroot process might steal a (terminating empty) chunk from the previous sequence,
        stalling some other non-root process indefinitely. In practice, this means \em wait should
        be set to true when sumToRoot() is called from probes and can be left to false when it is
        called from instruments. */
    static void sumToRoot(Array& arr, bool wait = false);

    /** This function broadcasts a separate sequence of floating point values from each process to
        the other processes. The chunk of data to be sent by the calling process must be generated
        by the provided call-back function \em producer. Similarly, the chunks of data reveived by
        the calling process from the other processes must be processed by the provided call-back
        function \em consumer. The chunks of data sent by each process may have different sizes.
        Using call-back functions (as apposed to passing the data directly) avoids the need for
        holding all data in memory at the same time, which can be a significant benefit in some
        cases. The following table lists the number of times each function is invoked.

        | nr of processes | \em producer invocations | \em consumer invocations |
        |--|--|--|
        | N=1 | 0 | 0   |
        | N>1 | 1 | N-1 |

        The \em producer function must store the data to be sent into the vector specified as its
        argument. The vector is guaranteed to be empty when the function is called (but it may have
        a nonzero memory allocation). Similarly, the \em consumer function can retrieve the
        received data from the vector specified as its argument.

        Within each process, the calls to the \em producer and \em consumer functions are
        serialized (i.e. they never occur in parallel threads). However, there is no ordering
        guarantee within or across processes other than that consuming data can happen only after
        it has been produced. This freedom allows a future implementation to use non-blocking
        communication primitives. */
    static void broadcastAllToAll(std::function<void(vector<double>& data)> producer,
                                  std::function<void(const vector<double>& data)> consumer);

    //======== Data members  ===========

private:
    static int _size;  // the number of processes in the run-time environment
    static int _rank;  // the rank of this process in the run-time environment
};

#endif
