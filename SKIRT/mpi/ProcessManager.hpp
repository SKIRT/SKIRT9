/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROCESSMANAGER_HPP
#define PROCESSMANAGER_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The ProcessManager class is an interface to the MPI library, and is the only place in the SKIRT
    code where explicit calls to this library are allowed. The ProcessManager class acts as a
    resource manager, with the resource being MPI. This resource is only available once at a time.
    Therefore, the ProcessManager ensures that the usage rights for MPI are only acquired when it
    is not yet requested at that time.

    A single instance of the ProcessManager class must be constructed just after program startup,
    and before the command-line arguments are used by other parts of the program. The latter is
    important because the MPI initialization functions called by the ProcessManager constructor may
    adjust the command-line arguments to remove MPI specific options. The program should not use
    the exit() or abort() functions, but rather let the main() function run to normal completion
    and return an exit code. Except for the constructor and destructor, all public functions in
    this class are static and thread-safe. */
class ProcessManager final
{
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

    /** This function gives the resource of being able to use MPI to the requesting object, on the
        condition that it is the first request. The object receives the total number of processes
        (Nprocs) and the rank of the process (rank). If other requests for MPI have preceeded the
        call, this will result in the number of processes and the rank to be 1 and 0 respectively.
        The number of active requests for MPI is checked with the 'requests' data member. */
    static void acquireMPI(int& Nprocs, int& rank);

    /** This function is called when an object does not longer need MPI functionality. A release of
        the MPI resource results in the number of requests to be decremented. Only when all
        requests have been released again, a new request for MPI with acquireMPI is answered with
        the correct Nprocs and rank. */
    static void releaseMPI();

    /** If this function is called, a process remains idle until all other processes have called it
        too. */
    static void barrier();

    /** This function is used to send a buffer consisting of double values to another process, of
        which the rank is specified with the 'receiver' parameter. A tag must also be passed to
        this function, which can be used to reveal some information about the context or purpose of
        the sent message to the receiving process. Only the sending process calls this function. */
    static void sendDoubleBuffer(const double* buffer, size_t count, int receiver, int tag);

    /** This function is used to receive a buffer consisting of double values from another,
        arbitrary process. The last argument is an integer where the rank of the sending process
        should be stored. Only the receiving process calls this function. */
    static void receiveDoubleBuffer(double* buffer, size_t count, int& sender);

    /** This function is used to receive a buffer consisting of double values from a particular
        other process. The rank of this sending process is specified with the third argument. Also
        passed to this function is an integer where the tag of the received message should be
        stored. Only the receiving process calls this function. */
    static void receiveDoubleBuffer(double* buffer, size_t count, int sender, int& tag);

    /** This function gathers a number of doubles from all processes at a certain receiving rank.
        The user can specify a pattern that will determine where exactly the received doubles will
        be placed in the receive buffer. This pattern consists blocks of equal length, placed in
        certain positions that depend on the process the data was received from. The receival
        positions per sending process are given through the displacements argument. This argument
        should contain as many lists as there are processes involved in this communication, and
        each list should consist of the locations where the blocks should start, in units of the
        block length. All processes must call this function for the communication to proceed. */
    static void gatherWithPattern(const double* sendBuffer, size_t sendCount,
                                  double* recvBuffer, int recvRank, size_t recvLength,
                                  const vector<vector<int>>& recvDisplacements);

    /** This function lets all processes send and receive an amount of double values. The arguments
        provide the necessary flexibility to handle non-contiguous data. The user can specify any
        pattern within the buffers that indicates which data to send to and receive from each
        process individually. The count arguments determine the number of times the patterns will
        be repeated. The shape of such a pattern consists of blocks of which the length is
        determined by the length arguments. The starting point of each block is placed at a
        multiple of the block length, and the amount empty space between the blocks will therefore
        also be a multiple of the block length. The exact locations where the blocks need to be
        placed can be specified using the nested lists of displacements. Each of the two
        'displacements' arguments should contain a list of block displacements, in units of the
        block length, for each process individually. This means that the number of lists should be
        equal to the number of processes. The last arguments are the extents each of the patterns
        should have. These are especially important when the corresponding count argument is
        greater than 1, as then they will determine where the next iteration of a pattern will
        start. All processes must call this function for the communication to proceed. */
    static void displacedBlocksAllToAll(const double* sendBuffer, size_t sendCount, size_t sendLength,
                                        const vector<vector<int>>& sendDisplacements, size_t sendExtent,
                                        double* recvBuffer, size_t recvCount, size_t recvLength,
                                        const vector<vector<int>>& recvDisplacements, size_t recvExtent);

    /** The purpose of this function is to sum a particular array of double values element-wise
        across the different processes. The resulting values are stored in the array passed as the
        second argument 'result_array', only on the process that is assigned as root. The rank of
        this particular process is specified in the third argument. All processes must call this
        function for the communication to proceed. */
    static void sum(double* my_array, size_t nvalues, int root);

    /** The purpose of this function is to sum a particular array of double values element-wise across
        the different processes. The resulting values are stored in the original array passed to this
        function, on each individual process. All processes must call this function for the
        communication to proceed. */
    static void sumAll(double* my_array, size_t nvalues);

    /** This function performs a reduction of a given boolean, by applying the logical OR operator
        across all processes. The result will overwrite the original boolean to which a pointer was
        passed. All processes must call this function for the communication to proceed. */
    static void orAll(bool* boolean);

    /** This function is used to broadcast an array of double values from one process to all other
        processes. A pointer to the first value is passed as the first argument, the number of
        values as the second and the rank of the sending process as the final argument. All
        processes must call this function for the communication to proceed. The array passed to
        this function by the receiving processes gets overwritten during the communication with the
        values stored in the array passed by the root. */
    static void broadcast(double* my_array, size_t nvalues, int root);

    /** This function is used to broadcast a single integer value from one process to all other
        processes. A pointer to the value is passed as the first argument and the rank of the sending
        process as the second. All processes must call this function for the communication to proceed.
        The memory where the integer value that is passed to this function by the receiving processes
        is stored, will be overwritten by the value stored in the memory of the process with rank \c
        root during the communication. */
    static void broadcast(int* value, int root);

    /** This function returns a boolean indicating whether the process is assigned as root or not.
        The rank of the process is always the 'true' rank, irrespective of whether the object that
        calls this function has acquired the MPI resource or not. */
    static bool isRoot();

    /** This function returns a boolean indicating whether multiple MPI processes are present. The
        number of processes is directly inquired from the MPI library. Therefore, the result of
        this function is unaffected by whether the object that calls it has acquired the MPI
        resource or not. */
    static bool isMultiProc();
};

#endif
