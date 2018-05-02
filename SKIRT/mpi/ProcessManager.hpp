/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROCESSMANAGER_HPP
#define PROCESSMANAGER_HPP

#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** The ProcessManager class is an interface to the MPI library, and is the only place in the SKIRT
    code where explicit calls to this library are allowed. A single instance of the ProcessManager
    class must be constructed just after program startup, and before the command-line arguments are
    used by other parts of the program. The latter is important because the MPI initialization
    functions called by the ProcessManager constructor may adjust the command-line arguments to
    remove MPI specific options. The program should not use the exit() or abort() functions, but
    rather let the main() function run to normal completion and return an exit code.

    The functions of this class should be called only from the main thread of the program.
    Violation of this rule causes undefined behavior that may differ between MPI implementations.
    */
class ProcessManager final
{
    //======== Construction - Destruction  ===========

private:
    /** This static function initializes the MPI library used for remote communication, if present.
        The function is called from the constructor of this class. */
    static void initialize(int *argc, char ***argv);

    /** This static function finalizes the MPI library used for remote communication, if present.
        The function is called from the destructor of this class. */
    static void finalize();

public:
    /** The constructor initializes the MPI library used for remote communication, if present. It
        is passed a reference to the command line arguments to provide the MPI library with the
        opportunity to use and/or adjust them (e.g. to remove any MPI related arguments). */
    ProcessManager(int *argc, char ***argv) { initialize(argc, argv); }

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
    static bool isRoot() { return _rank==0; }

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
        proceed. If there is only one process, the function does nothing. */
    static void sumToAll(Array& arr);

    /** This function adds the floating point values of an array element-wise across the different
        processes. The resulting sums are then stored in the same Array passed to this function on
        the root process. The arrays on the other processes are left untouched. All processes must
        call this function for the communication to proceed. If there is only one process, the
        function does nothing. */
    static void sumToRoot(Array& arr);

    //======== Data members  ===========

private:
    static int _size;    // the number of processes in the run-time environment
    static int _rank;    // the rank of this process in the run-time environment
};

#endif
