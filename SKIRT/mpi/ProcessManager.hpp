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

    //======== Communication  ===========

    /** This function causes the calling process to block until all other processes have invoked it
        as well. If there is only one process, the function does nothing. */
    static void wait();

    /** This function adds the floating point values of an array element-wise across the different
        processes. The resulting sums are then stored in the same Array passed to this function on
        each individual process. All processes must call this function for the communication to
        proceed. */
    static void sumToAll(Array& arr);

    /** This function adds the floating point values of an array element-wise across the different
        processes. The resulting sums are then stored in the same Array passed to this function on
        the root process. The arrays on the other processes are left untouched. All processes must
        call this function for the communication to proceed. */
    void sumToRoot(Array& arr);

    //======== Data members  ===========

private:
    static int _size;    // the number of processes in the run-time environment
    static int _rank;    // the rank of this process in the run-time environment
};

#endif
