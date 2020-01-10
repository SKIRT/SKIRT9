/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "Basics.hpp"
#include <fstream>

////////////////////////////////////////////////////////////////////

/** The System class provides a set of functions that interact with the operating system or with
    the file system, often in a platform-dependent manner. One important function of the class is
    translating internal UTF-8 encoded strings to and from the encoding expected by the host
    system.

    A single instance of the System class must be constructed just after program startup, and
    certainly before any parallel threads are started. A good place is the first line of the main()
    function. The program should not use the exit() or abort() functions, but rather let the main()
    function run to normal completion and return an exit code. Except for the constructor and
    destructor, all public functions in this class are static and thread-safe. */
class System final
{

    // ================== Initialization ==================

private:
    /** This function performs initialization for all system-specific facilities used by a client
        of this class. Specifically, it retrieves the command-line arguments either from the main()
        arguments passed on to this function or through a platform-specific method, sets the C and
        C++ locales to the portable "classic C" locale, and initializes the console for reading and
        writing. This function is \em not thread-safe; it should be invoked exactly once at program
        startup, before starting any parallel threads and before using any of the facilities
        offered by this class. */
    static void initialize(int argc, char** argv);

    /** This function performs clean-up for all system-specific facilities offered by this class.
        Specifically, it releases any resources and restores any external changes allocated or
        affected by the initialize() function. This function is \em not thread-safe; it should be
        invoked exactly once at program termination, after closing down any parallel threads and
        after the last use of the facilities offered by this class. */
    static void finalize();

public:
    /** This constructor performs initialization for all system-specific facilities used by a
        client of this class. For more information, refer to the description of the initialize()
        function. */
    System(int argc, char** argv) { initialize(argc, argv); }

    /** The destructor performs clean-up for all system-specific facilities offered by this class.
        For more information, refer to the description of the finalize() function. */
    ~System() { finalize(); }

    /** The copy constructor is deleted because instances of this class should never be copied or
        moved. */
    System(const System&) = delete;

    /** The assignment operator is deleted because instances of this class should never be copied
        or moved. */
    System& operator=(const System&) = delete;

    // ================== Environment ==================

    /** This function returns a list of the command line arguments, excluding the program's
        invocation path. */
    static vector<string> arguments();

    /** This function returns the absolute file path of the program's executable, or the empty
        string if the path cannot be determined. On Windows, Mac and Linux, the function works
        regardless of how the executable has been invoked. On other systems, the function requires
        the executable to be invoked with the absolute file path. */
    static string executablePath();

    /** This function returns the local host name, or the empty string if the local host name
        cannot be determined. */
    static string hostname();

    /** This function returns the name of the user running this executable, or the empty string if
        the user name cannot be determined. */
    static string username();

    /** This function returns the current local time with millisecond resolution as a string in one
        of two formats. If the specified flag is false or missing, the format is intended for human
        consumption ("DD/MM/YYYY hh:mm:ss.mmm"). If the flag is true, the format conforms to the
        ISO 8601 standard used for XML datetime fields ("YYYY-MM-DDThh:mm:ss.mmm"). */
    static string timestamp(bool iso8601 = false);

    // ================== Console ==================

    /** This enumeration includes a constant for each logging level, in increasing order of
        importance. */
    enum class LogLevel { Info, Warning, Success, Error };

    /** This function writes the specified message to the console, adding a time stamp and using a
        color that depends on the specified log level (if the console supports color). */
    static void log(string message, LogLevel level = LogLevel::Info);

    /** This function prompts for and returns user input on the console using the specified
        message. */
    static string prompt(string message);

    // ================== File System ==================

    /** This function returns an input file stream opened on the specified file path. On Windows
        the function replaces forward slashes in the file path by backward slashes. */
    static std::ifstream ifstream(string path);

    /** This function returns an output file stream opened on the specified file path. If a file
        already exists at the specified path, by default it is overwritten. However, if the \em
        append flag is specified and is true, new output will be appended to the existing file. On
        Windows the function replaces forward slashes in the file path by backward slashes. */
    static std::ofstream ofstream(string path, bool append = false);

    /** This function returns true if the specified path refers to an existing regular file. On
        Windows the function replaces forward slashes in the path by backward slashes. */
    static bool isFile(string path);

    /** This function returns true if the specified path refers to an existing directory. The empty
        string is interpreted as the current directory. On Windows the function replaces forward
        slashes in the path by backward slashes. */
    static bool isDir(string path);

    /** This function creates a new folder with the specified path, if it does not already exist.
        All path segments other than the last one should correspond to already existing
        directories. The function returns true if the directory already existed or was successfully
        created; otherwise it returns false. On Windows the function replaces forward slashes in
        the path by backward slashes. */
    static bool makeDir(string directory);

    /** This function removes the file with the specified path. It does nothing if the file does
        not exist or can't be removed. On Windows the function replaces forward slashes in the path
        by backward slashes. */
    static void removeFile(string path);

    /** This function returns the names for all regular files residing in the given directory,
        specified as an absolute or relative path without trailing slash, or the empty string for
        the current directory. On Windows the function replaces forward slashes in the path by
        backward slashes. The returned list is sorted alphabetically. If the given directory does
        not exist or can't be accessed, an empty list is returned. */
    static vector<string> filesInDirectory(string directory);

    /** This function returns the names for all subdirectories residing in the given directory,
        specified as an absolute or relative path without trailing slash, or the empty string for
        the current directory. On Windows the function replaces forward slashes in the path by
        backward slashes. The returned list is sorted alphabetically. If the given directory does
        not exist or can't be accessed, an empty list is returned. */
    static vector<string> dirsInDirectory(string directory);

    /** If the conversion process succeeds, this function returns a canonical absolute file path
        corresponding to the specified absolute or relative path, after resolving symbolic links
        and ".." or "." path segments (the empty string is interpreted as the current directory).
        This requires most if not all of the path segments to actually exist in the file system. If
        the conversion process fails, the function returns the original path. On Windows the
        function replaces forward slashes in the specified file path by backward slashes before
        processing the path. */
    static string canonicalPath(string path);

    // ================== Mapped File I/O ==================

    /** This function acquires a read-only memory map on the specified file. In other words, the
        contents of the file is mapped directly into the memory space of the calling process. If
        successful, the function returns the address and the length (in bytes) of the memory map.
        If the memory map cannot be acquired, the returned address is the null pointer and the
        returned length is zero. On Windows the function replaces forward slashes in the specified
        file path by backward slashes.

        It is allowed to acquire a memory map on the same file more than once (in the same or in
        different execution threads). As long as the previous memory map on the file has not been
        released, the same memory map is returned for subsequent acquisitions. */
    static std::pair<void*, size_t> acquireMemoryMap(string path);

    /** This function releases a previously acquired memory map on the specified file. Releasing a
        memory map invalidates any and all pointers into the memory range previously occupied by
        the memory map. On Windows the function replaces forward slashes in the specified file path
        by backward slashes.

        If multiple memory maps have been acquired on the same file, a matching number of release
        operations is needed to actually release the memory map. */
    static void releaseMemoryMap(string path);

    // ================== Debugging ==================

    /** This function returns a list of lines representing a stack trace to the current execution
        frame. The format should be semi-readable for humans. */
    static vector<string> stacktrace();

    // ================== Memory ==================

    /** Returns the size of physical memory (RAM) available on the system in bytes, or zero if the
        value cannot be determined. */
    static size_t availableMemory();

    /** Returns the peak (maximum so far) physical memory use for the current process in bytes, or
        zero if the value cannot be determined. */
    static size_t peakMemoryUsage();

    /** Returns the current physical memory use for the current process in bytes, or zero if the
        value cannot be determined. */
    static size_t currentMemoryUsage();
};

////////////////////////////////////////////////////////////////////

#endif
