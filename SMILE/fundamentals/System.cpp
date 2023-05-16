/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "System.hpp"

////////////////////////////////////////////////////////////////////

#include <chrono>
#include <clocale>
#include <ctime>
#include <locale>
#include <mutex>
#include <unordered_map>

#ifdef _WIN64
#    include <Windows.h>
#    include <Lmcons.h>   // for getting user name
#    include <Pathcch.h>  // for getting canonical path
#    ifdef _DEBUG
#        include <Dbghelp.h>  // for stack trace
#    endif
#else
#    include <dirent.h>    // for reading directories
#    include <execinfo.h>  // for stack trace
#    include <fcntl.h>     // for opening files (low-level)
#    include <sys/mman.h>  // for memory mapped files
#    include <sys/stat.h>  // for reading file status
#    include <unistd.h>    // for gethostname
#    include <iostream>
#endif

#if defined(__APPLE__) && defined(__MACH__)
#    include <CoreFoundation/CoreFoundation.h>
#endif

////////////////////////////////////////////////////////////////////

namespace
{
    // List of command line arguments retrieved in initialize()
    string argument0;
    vector<string> _arguments;

    // Mutex to guard console input/output operations
    std::mutex _consoleMutex;

    // The strings beginning and ending a message, indexed by level (Info, Warning, Success, Error, Prompt=Error+1)
    // NOTE: this depends on the order in the Level enum --> rather dirty
    const char* _messageBegin[] = {"   ", " ! ", " - ", " * ", " ? "};
    const char* _messageEnd[] = {"\n", "\n", "\n", "\n", ": "};

#ifdef _WIN64
    // The console output code page saved in initialize() and restored in finalize()
    UINT _oldOutputCodePage;

    // Windows character attributes for coloring, indexed by level (Info, Warning, Success, Error, Prompt=Error+1)
    // NOTE: this depends on the order in the Level enum --> rather dirty
    const WORD MY_WHITE = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE;
    const WORD MY_CYAN = FOREGROUND_GREEN | FOREGROUND_BLUE;
    const WORD MY_MAGENTA = FOREGROUND_RED | FOREGROUND_BLUE;
    WORD _colorBegin[] = {MY_WHITE, MY_MAGENTA, FOREGROUND_GREEN, FOREGROUND_RED | FOREGROUND_INTENSITY, MY_CYAN};
    WORD _colorEnd[] = {MY_WHITE, MY_WHITE, MY_WHITE, MY_WHITE, MY_WHITE};

    // The current directory retrieved at program startup (Windows only) and converted to UTF-8 string
    // (because getting the current directory or using relative paths is not thread-safe under Windows)
    string currentDir;
#else
    // Set to true in initialize() if the console supports ANSI escape sequences for coloring
    bool _colored = false;

    // ANSI escape sequences for coloring, indexed by level (Info, Warning, Success, Error, Prompt=Error+1)
    // NOTE: this depends on the order in the Level enum --> rather dirty
    const char* _colorBegin[] = {"", "\033[35m", "\033[32m", "\033[31m", "\033[34m"};
    const char* _colorEnd[] = {"", "\033[0m", "\033[0m", "\033[0m", "\033[0m"};
#endif

    // Arbitrary maximum path length; functions may fail or even crash if paths exceed this length
    constexpr int MAXPATHLEN = 8000;

    // Mutex to guard the acquisition and release of file memory maps
    std::mutex _mapMutex;

    // Data type for storing information on a currently acquired file memory map
    struct MapRecord
    {
        void* start{nullptr};  // pointer to start of memory map
        size_t length{0};      // length of memory map
        int count{0};          // current number of acquisitions for this map
#ifdef _WIN64
        HANDLE filehandle{INVALID_HANDLE_VALUE};
        HANDLE maphandle{NULL};
#else
        int filehandle{-1};
#endif
    };

    // Dictionary keeping track of all currently acquired file memory maps: <canonical_path, map_record>
    std::unordered_map<string, MapRecord> _maps;
}

////////////////////////////////////////////////////////////////////

#ifdef _WIN64
namespace
{
    // converts a zero-terminated 16-bit C-style string to a regular UTF-8 string object
    string toUTF8(wchar_t* input)
    {
        int utfLength = WideCharToMultiByte(CP_UTF8, 0, input, -1, 0, 0, 0, 0);
        auto utfInput = std::make_unique<char[]>(utfLength);
        WideCharToMultiByte(CP_UTF8, 0, input, -1, utfInput.get(), utfLength, 0, 0);
        return utfInput.get();
    }

    // converts a regular UTF-8 string object representing an absolute or relative file path
    // to a 16-bit C-style string encapsulated in unique pointer,
    // replacing forward slashes by backward slashes,
    // and converting a relative path to an absolute path using the current directory stored at initialization
    std::unique_ptr<wchar_t[]> toUTF16(string path)
    {
        // replace forward slashes by backward slashes
        replace(path.begin(), path.end(), '/', '\\');

        // if we have a current directory, convert a relative path to an absolute path
        if (!currentDir.empty())
        {
            if (path.empty())
                path = currentDir;
            else if (path.front() != '\\' && path.find(':') == string::npos)
                path = currentDir + '\\' + path;
        }

        // convert UTF-8 to UTF-16
        auto pathLength = path.length() + 1;  // include terminating zero
        auto widePath = std::make_unique<wchar_t[]>(pathLength);
        MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, widePath.get(), static_cast<int>(pathLength));
        return widePath;
    }
}
#endif

////////////////////////////////////////////////////////////////////

void System::initialize(int argc, char** argv)
{
    // Force standard locale so that snprintf and stream formatting always produces the same result
    std::locale::global(std::locale::classic());  // this may or may not affect C locale so do this first
    setlocale(LC_ALL, "C");

    // Becomes true if we're on Windows and the Unicode arguments have been successfully retrieved
    bool retrieved = false;

#ifdef _WIN64
    // Retrieve the Unicode version of the command line arguments
    int nArgs;
    LPWSTR* argList = CommandLineToArgvW(GetCommandLineW(), &nArgs);
    if (argList)
    {
        // The maximum required length for the UTF-8 buffer, for all arguments
        int maxLength = 0;

        // Test successful conversion of each argument to UTF-8, and keep track of max length
        for (int i = 1; i < nArgs; i++)
        {
            int thisLength = WideCharToMultiByte(CP_UTF8, 0, argList[i], -1, 0, 0, 0, 0);
            if (thisLength <= 0)
            {
                maxLength = -1;  // error -> do not retrieve Unicode arguments
                break;
            }
            maxLength = max(thisLength, maxLength);
        }

        // Convert each argument to UTF-8 and store the result
        if (maxLength >= 0)
        {
            auto utfArg = std::make_unique<char[]>(maxLength + 1);  // avoid constructing a zero-length array
            for (int i = 1; i < nArgs; i++)
            {
                WideCharToMultiByte(CP_UTF8, 0, argList[i], -1, utfArg.get(), maxLength, 0, 0);
                _arguments.emplace_back(utfArg.get());
            }
            retrieved = true;
        }

        // Free memory allocated by CommandLineToArgvW
        LocalFree(argList);
    }
#endif

    // If Windows Unicode arguments were not retrieved, store the C++ main() command line arguments
    if (!retrieved)
    {
        argument0 = argv[0];
        for (int i = 1; i < argc; ++i)
        {
            _arguments.emplace_back(argv[i]);
        }
    }

#ifdef _WIN64
    // Save the current output console code page for later restore
    // and set it to UTF-8 so that we can write UTF-8 to the console without conversion
    _oldOutputCodePage = GetConsoleOutputCP();
    SetConsoleOutputCP(CP_UTF8);

    // Get the current directory
    WCHAR buffer[MAXPATHLEN];
    DWORD size = GetCurrentDirectoryW(sizeof(buffer), buffer);
    if (size > 0 && size < sizeof(buffer) - 2) currentDir = toUTF8(buffer);
#else
    // We assume that coloring is supported if the TERM environment variable is defined
    _colored = getenv("TERM") != 0;
#endif
}

////////////////////////////////////////////////////////////////////

void System::finalize()
{
    // Release any remaining file memory mappings
    while (!_maps.empty()) System::releaseMemoryMap(_maps.begin()->first);

    // Clear and deallocate the list of command line arguments
    vector<string>().swap(_arguments);

#ifdef _WIN64
    // Restore the initial code page
    SetConsoleOutputCP(_oldOutputCodePage);
#endif
}

////////////////////////////////////////////////////////////////////

vector<string> System::arguments()
{
    return _arguments;
}

string System::executablePath()
{
#if defined(_WIN64)
    // Use Windows functions; there is no fallback because we don't store the zeroth command line argument,
    // and, besides, we would have problematic UTF-16/UTF-8 encoding issues
    WCHAR buffer[MAXPATHLEN];
    DWORD size = GetModuleFileNameW(0, buffer, sizeof(buffer));
    if (size > 0 && size < sizeof(buffer) - 2)
    {
        return toUTF8(buffer);
    }
    return string();

#else  // not windows
#    if defined(__APPLE__) && defined(__MACH__)
    // Try Mac OS X Core Foundation functions
    char buffer1[MAXPATHLEN];
    buffer1[0] = 0;
    CFURLRef bundleURL = CFBundleCopyExecutableURL(CFBundleGetMainBundle());
    if (bundleURL)
    {
        CFStringRef bundlePath = CFURLCopyFileSystemPath(bundleURL, kCFURLPOSIXPathStyle);
        if (bundlePath)
        {
            if (!CFStringGetCString(bundlePath, buffer1, sizeof(buffer1), kCFStringEncodingUTF8)) buffer1[0] = 0;
            CFRelease(bundlePath);
        }
        CFRelease(bundleURL);
    }
    char buffer2[MAXPATHLEN];
    if (realpath(buffer1, buffer2)) return buffer2;
#    elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    // Try Linux /proc/<pid>/exe symlink that points to the absolute path of the executable
    char buffer[MAXPATHLEN];
    buffer[0] = 0;
    string symlink = "/proc/" + std::to_string(getpid()) + "/exe";
    if (realpath(symlink.c_str(), buffer)) return buffer;
#    endif  // not Mac or Linux

    // If the Mac or Linux platform attempts failed, or if we're on another platform,
    // return the zeroth command line argument if it happens to be an absolute path
    if (argument0.size() && argument0[0] == '/') return argument0;
    return string();
#endif
}

////////////////////////////////////////////////////////////////////

string System::hostname()
{
#ifdef _WIN64
    // Initialize the winsock library that provides gethostname()
    WSAData wsadata;
    WSAStartup(MAKEWORD(2, 0), &wsadata);
#endif

    // Get the hostname; this function is strangely cross-platform!
    char buffer[256];
    if (gethostname(buffer, sizeof(buffer)) == 0)
        buffer[sizeof(buffer) - 1] = 0;
    else
        buffer[0] = 0;

#ifdef _WIN64
    // Release any winsock library resources
    WSACleanup();
#endif
    return buffer;
}

////////////////////////////////////////////////////////////////////

string System::username()
{
#ifdef _WIN64
    TCHAR buffer[UNLEN + 1];
    DWORD size = UNLEN + 1;
    if (GetUserName(buffer, &size) != 0) return buffer;
#else
    char buffer[256];
    if (getlogin_r(buffer, sizeof(buffer)) == 0) return buffer;
#endif

    // If the platform-specific function fails, try some environment variables
    char* pointer;
    if ((pointer = getenv("LOGNAME"))) return pointer;
    if ((pointer = getenv("USERNAME"))) return pointer;
    if ((pointer = getenv("USER"))) return pointer;
    return string();
}

////////////////////////////////////////////////////////////////////

// The thread-safe version of localtime is called differently on Unix and Windows,
// and the arguments are in reverse order. So define a macro...
#ifdef _WIN64
#    define TO_LOCAL_TIME(from, to) localtime_s(to, from)
#else
#    define TO_LOCAL_TIME(from, to) localtime_r(from, to)
#endif

string System::timestamp(bool iso8601)
{
    using namespace std::chrono;

    // setup the output format and buffer
    const char* format = iso8601 ? "%Y-%m-%dT%H:%M:%S.xxx" : "%d/%m/%Y %H:%M:%S.xxx";
    const size_t resultLength = 23;
    const size_t milliOffset = 20;
    char resultBuf[resultLength + 2];  // add 2 rather than 1 to avoid warning by GCC v8.1

    // get the current wall time
    system_clock::time_point now_tp = system_clock::now();

    // convert to calendar up to second resolution
    time_t now_tt = system_clock::to_time_t(now_tp);
    struct tm now_tm;
    TO_LOCAL_TIME(&now_tt, &now_tm);
    strftime(resultBuf, sizeof(resultBuf), format, &now_tm);

    // get the millisecond fraction, assuming that the epoch starts at integer-second time point
    system_clock::duration since_epoch = now_tp.time_since_epoch();
    since_epoch -= duration_cast<seconds>(since_epoch);
    milliseconds millis = duration_cast<milliseconds>(since_epoch);
    snprintf(resultBuf + milliOffset, sizeof(resultBuf) - milliOffset, "%03d", static_cast<int>(millis.count()));

    return resultBuf;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function outputs the specified message at the given log level index, adding all decorations
    // The function is NOT threadsafe, because the system console I/O is not threadsafe.
    // Thus: call this function only within the critical section setup for the console I/O.
    void outputMessage(string message, size_t level)
    {
#ifdef _WIN64
        SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), _colorBegin[level]);
        fputs(System::timestamp().c_str(), stdout);
        fputs(_messageBegin[level], stdout);
        fputs(message.c_str(), stdout);
        fputs(_messageEnd[level], stdout);
        fflush(stdout);
        SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), _colorEnd[level]);
#else
        std::cout << (_colored ? _colorBegin[level] : "") << System::timestamp() << _messageBegin[level] << message
                  << _messageEnd[level] << (_colored ? _colorEnd[level] : "") << std::flush;
#endif
    }
}

////////////////////////////////////////////////////////////////////

void System::log(string message, LogLevel level)
{
    std::unique_lock<std::mutex> lock(_consoleMutex);
    outputMessage(message, static_cast<size_t>(level));  // dirty cast
}

////////////////////////////////////////////////////////////////////

string System::prompt(string message)
{
    std::unique_lock<std::mutex> lock(_consoleMutex);
    outputMessage(message, static_cast<size_t>(LogLevel::Error) + 1);  // dirty cast

#ifdef _WIN64
    // Get user input and convert it from UTF-16 to UTF-8
    const DWORD MAXLEN = 2000;
    WCHAR input[MAXLEN];
    DWORD numRead = 0;
    BOOL success = ReadConsoleW(GetStdHandle(STD_INPUT_HANDLE), input, MAXLEN, &numRead, 0);
    if (!success || numRead < 2) return "";
    input[numRead - 2] = 0;
    return toUTF8(input);
#else
    string input;
    std::getline(std::cin, input);
    return input;
#endif
}

////////////////////////////////////////////////////////////////////

std::ifstream System::ifstream(string path)
{
#ifdef _WIN64
    return std::ifstream(toUTF16(path).get());
#else
    return std::ifstream(path);
#endif
}

////////////////////////////////////////////////////////////////////

std::ofstream System::ofstream(string path, bool append)
{
#ifdef _WIN64
    return std::ofstream(toUTF16(path).get(), append ? std::ios_base::app : std::ios_base::out);
#else
    return std::ofstream(path, append ? std::ios_base::app : std::ios_base::out);
#endif
}

////////////////////////////////////////////////////////////////////

bool System::isFile(string path)
{
#ifdef _WIN64
    DWORD attrs = GetFileAttributesW(toUTF16(path).get());
    if (attrs != INVALID_FILE_ATTRIBUTES && !(attrs & FILE_ATTRIBUTE_DIRECTORY)) return true;
#else
    struct stat st;
    if (!stat(path.c_str(), &st) && S_ISREG(st.st_mode)) return true;
#endif
    return false;
}

////////////////////////////////////////////////////////////////////

bool System::isDir(string path)
{
    // empty string means current directory
    if (path.empty()) path += '.';

#ifdef _WIN64
    DWORD attrs = GetFileAttributesW(toUTF16(path).get());
    if (attrs != INVALID_FILE_ATTRIBUTES && (attrs & FILE_ATTRIBUTE_DIRECTORY)) return true;
#else
    struct stat st;
    if (!stat(path.c_str(), &st) && S_ISDIR(st.st_mode)) return true;
#endif
    return false;
}

////////////////////////////////////////////////////////////////////

bool System::makeDir(string directory)
{
    if (isDir(directory)) return true;

#ifdef _WIN64
    return CreateDirectoryW(toUTF16(directory).get(), NULL) != 0;
#else
    return mkdir(directory.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == 0;
#endif
}

////////////////////////////////////////////////////////////////////

void System::removeFile(string path)
{
    // guard against removing something that is not a regular file, such as a directory
    if (!isFile(path)) return;

#ifdef _WIN64
    DeleteFileW(toUTF16(path).get());
#else
    remove(path.c_str());
#endif
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function returns the names for all regular files or directories residing in the given directory
    // according to the semantics described for the filesInDirectory() and dirsInDirectory() functions.
    void itemsInDirectory(vector<string>& result, string directory, bool files)
    {
        // empty string means current directory
        if (directory.empty()) directory += '.';

        // lock although functions are probably thread-safe in modern library implementations
        std::unique_lock<std::mutex> lock(_consoleMutex);

#ifdef _WIN64
        // open directory stream
        HANDLE dir;
        WIN32_FIND_DATAW filedata;
        if ((dir = FindFirstFileW(toUTF16(directory + "\\*").get(), &filedata)) != INVALID_HANDLE_VALUE)
        {
            // read entries
            do
            {
                // skip anything that is not a regular file or directory
                if (!(filedata.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) == files)
                {
                    string candidate = toUTF8(filedata.cFileName);
                    // skip special directories
                    if (!candidate.empty() && candidate[0] != '.') result.push_back(toUTF8(filedata.cFileName));
                }
            } while (FindNextFileW(dir, &filedata));

            FindClose(dir);
        }
#else
        // open directory stream
        DIR* dir = opendir(directory.c_str());
        if (dir)
        {
            // read the next entry
            struct dirent* ent;
            while ((ent = readdir(dir)))
            {
                // skip hidden files
                if (ent->d_name[0] != '.')
                {
                    string filename = ent->d_name;
                    string filepath = directory + "/" + filename;

                    // skip anything that is not a regular file or directory
                    struct stat st;
                    if (!stat(filepath.c_str(), &st) && (files ? S_ISREG(st.st_mode) : S_ISDIR(st.st_mode)))
                    {
                        result.push_back(filename);
                    }
                }
            }
            closedir(dir);
        }
#endif

        // sort the list alphabetically
        std::sort(result.begin(), result.end());
    }
}

////////////////////////////////////////////////////////////////////

vector<string> System::filesInDirectory(string directory)
{
    vector<string> result;
    itemsInDirectory(result, directory, true);
    return result;
}

////////////////////////////////////////////////////////////////////

vector<string> System::dirsInDirectory(string directory)
{
    vector<string> result;
    itemsInDirectory(result, directory, false);
    return result;
}

////////////////////////////////////////////////////////////////////

string System::canonicalPath(string path)
{
    // empty string means current directory
    if (path.empty()) path += '.';

#if defined(_WIN64)
    WCHAR buffer[MAXPATHLEN];
    if (S_OK == PathCchCanonicalize(buffer, sizeof(buffer), toUTF16(path).get())) return toUTF8(buffer);
#else
    char buffer[MAXPATHLEN];
    if (realpath(path.c_str(), buffer)) return buffer;
#endif
    return path;
}

////////////////////////////////////////////////////////////////////

std::pair<void*, size_t> System::acquireMemoryMap(string path)
{
    // use the canonical path as a unique identifier for the file
    path = canonicalPath(path);

    // perform the mapping operation in a critical section because
    //  - access to our map dictionary is certainly not thread-safe
    //  - thread-safety of the operating system calls is not so clear
    std::unique_lock<std::mutex> lock(_mapMutex);

    // make a new map only if we don't have one cached
    if (!_maps.count(path))
    {
        MapRecord record;

        // attempt to acquire the map:
        //   - when this fails, the start field remains nullptr and other fields are unspecified
        //   - when this succeeds, all fields will be properly initialized (and thus start will be nonzero)

#if defined(_WIN64)
        // open the file
        record.filehandle =
            CreateFileW(toUTF16(path).get(), GENERIC_READ, FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
        if (record.filehandle != INVALID_HANDLE_VALUE)
        {
            // get the file size
            LARGE_INTEGER filesize;
            if (GetFileSizeEx(record.filehandle, &filesize)) record.length = filesize.QuadPart;

            if (record.length)
            {
                // create the map object (without mapping it into memory)
                record.maphandle = CreateFileMapping(record.filehandle, 0, PAGE_READONLY, 0, 0, 0);
                if (record.maphandle != NULL)
                {
                    // actually map it into memory
                    record.start = MapViewOfFile(record.maphandle, FILE_MAP_READ, 0, 0, 0);
                }
            }

            // if map acquisition failed, close any successful handles
            if (!record.start)
            {
                if (record.maphandle != NULL) CloseHandle(record.maphandle);
                CloseHandle(record.filehandle);
            }
        }
#else
        // open the file
        record.filehandle = open(path.c_str(), O_RDONLY);
        if (record.filehandle != -1)
        {
            // get the file size
            struct stat filesize;
            if (fstat(record.filehandle, &filesize) != -1) record.length = filesize.st_size;

            if (record.length)
            {
                // create the mapping
                auto start = mmap(0, record.length, PROT_READ, MAP_PRIVATE, record.filehandle, 0);
                if (start != MAP_FAILED) record.start = start;
            }

            // if map acquisition failed, close the file
            if (!record.start)
            {
                close(record.filehandle);
            }
        }
#endif

        // if successful, add the map to our dictionary; otherwise return "failure"
        if (record.start)
            _maps.emplace(path, record);
        else
            return std::make_pair(nullptr, 0);
    }

    // retrieve and return the new or existing map
    MapRecord& record = _maps.at(path);
    record.count++;
    return std::make_pair(record.start, record.length);
}

////////////////////////////////////////////////////////////////////

void System::releaseMemoryMap(string path)
{
    if (_maps.count(path))
    {
        // perform the mapping operation in a critical section because
        //  - access to our map dictionary is certainly not thread-safe
        //  - thread-safety of the operating system calls is not so clear
        std::unique_lock<std::mutex> lock(_mapMutex);

        MapRecord& record = _maps.at(path);
        record.count--;

        // actually release the memory map only when count has reached zero
        if (!record.count)
        {

#if defined(_WIN64)
            UnmapViewOfFile(record.start);
            CloseHandle(record.maphandle);
            CloseHandle(record.filehandle);
#else
            munmap(record.start, record.length);
            close(record.filehandle);
#endif

            // remove the map entry from our dictionary
            _maps.erase(path);
        }
    }
}

////////////////////////////////////////////////////////////////////

vector<string> System::stacktrace()
{
    vector<string> result;
    result.emplace_back("Call stack:");

#ifdef _WIN64
#    ifdef _DEBUG
    void* stack[62];
    USHORT frames = CaptureStackBackTrace(0, 62, stack, NULL);
    HANDLE process = GetCurrentProcess();
    SymInitialize(process, NULL, TRUE);
    SYMBOL_INFO* symbol = (SYMBOL_INFO*)calloc(sizeof(SYMBOL_INFO) + 256 * sizeof(char), 1);
    symbol->MaxNameLen = 255;
    symbol->SizeOfStruct = sizeof(SYMBOL_INFO);
    for (USHORT i = 0; i < frames; i++)
    {
        SymFromAddr(process, (DWORD64)(stack[i]), 0, symbol);
        string symbolName(symbol->Name);
        if (symbolName.empty()) symbolName = "unknown symbol";
        result.push_back(std::to_string(frames - i - 1) + " " + symbolName);
    }
    free(symbol);
#    else
    result.push_back("Sorry, stack trace is implemented only in DEBUG mode");
#    endif  // DEBUG
#else
    // Get the strack trace as a list of C-style strings
    const int max_depth = 100;
    void* stack_addrs[max_depth];
    int stack_depth = backtrace(stack_addrs, max_depth);
    char** stack_strings = backtrace_symbols(stack_addrs, stack_depth);

    // Remove contiguous spaces and store into result
    for (int i = 1; i < stack_depth; i++)
    {
        string line(stack_strings[i]);
        auto new_end = std::unique(line.begin(), line.end(), [](char a, char b) { return a == ' ' && b == ' '; });
        line.erase(new_end, line.end());
        result.push_back(line);
    }

    // Free memory // malloc()ed by backtrace_symbols() function
    free(stack_strings);
#endif

    return result;
}

////////////////////////////////////////////////////////////////////

namespace
{

    /*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#    include <Windows.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#    include <sys/param.h>
#    include <sys/types.h>
#    include <unistd.h>
#    if defined(BSD)
#        include <sys/sysctl.h>
#    endif
#endif

    /**
 * Returns the size of physical memory (RAM) in bytes.
 */
    size_t getMemorySize()
    {
#if defined(_WIN32) && (defined(__CYGWIN__) || defined(__CYGWIN32__))
        /* Cygwin under Windows. ------------------------------------ */
        /* New 64-bit MEMORYSTATUSEX isn't available.  Use old 32.bit */
        MEMORYSTATUS status;
        status.dwLength = sizeof(status);
        GlobalMemoryStatus(&status);
        return static_cast<size_t>(status.dwTotalPhys);

#elif defined(_WIN32)
        /* Windows. ------------------------------------------------- */
        /* Use new 64-bit MEMORYSTATUSEX, not old 32-bit MEMORYSTATUS */
        MEMORYSTATUSEX status;
        status.dwLength = sizeof(status);
        GlobalMemoryStatusEx(&status);
        return static_cast<size_t>(status.ullTotalPhys);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
        /* UNIX variants. ------------------------------------------- */
        /* Prefer sysctl() over sysconf() except sysctl() HW_REALMEM and HW_PHYSMEM */

#    if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))
        int mib[2];
        mib[0] = CTL_HW;
#        if defined(HW_MEMSIZE)
        mib[1] = HW_MEMSIZE; /* OSX. --------------------- */
#        elif defined(HW_PHYSMEM64)
        mib[1] = HW_PHYSMEM64; /* NetBSD, OpenBSD. --------- */
#        endif
        int64_t size = 0;    /* 64-bit */
        size_t len = sizeof(size);
        if (sysctl(mib, 2, &size, &len, NULL, 0) == 0) return static_cast<size_t>(size);
        return 0L; /* Failed? */

#    elif defined(_SC_AIX_REALMEM)
        /* AIX. ----------------------------------------------------- */
        return static_cast<size_t>(sysconf(_SC_AIX_REALMEM)) * static_cast<size_t>(1024);

#    elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
        /* FreeBSD, Linux, OpenBSD, and Solaris. -------------------- */
        return static_cast<size_t>(sysconf(_SC_PHYS_PAGES)) * static_cast<size_t>(sysconf(_SC_PAGESIZE));

#    elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
        /* Legacy. -------------------------------------------------- */
        return static_cast<size_t>(sysconf(_SC_PHYS_PAGES)) * static_cast<size_t>(sysconf(_SC_PAGE_SIZE));

#    elif defined(CTL_HW) && (defined(HW_PHYSMEM) || defined(HW_REALMEM))
        /* DragonFly BSD, FreeBSD, NetBSD, OpenBSD, and OSX. -------- */
        int mib[2];
        mib[0] = CTL_HW;
#        if defined(HW_REALMEM)
        mib[1] = HW_REALMEM;   /* FreeBSD. ----------------- */
#        elif defined(HW_PYSMEM)
        mib[1] = HW_PHYSMEM; /* Others. ------------------ */
#        endif
        unsigned int size = 0; /* 32-bit */
        size_t len = sizeof(size);
        if (sysctl(mib, 2, &size, &len, NULL, 0) == 0) return static_cast<size_t>(size);
        return 0L; /* Failed? */
#    endif /* sysctl and sysconf variants */

#else
        return 0L; /* Unknown OS. */
#endif
    }

}  // end anonymous namespace

////////////////////////////////////////////////////////////////////

namespace
{

    /*
     * Author:  David Robert Nadeau
     * Site:    http://NadeauSoftware.com/
     * License: Creative Commons Attribution 3.0 Unported License
     *          http://creativecommons.org/licenses/by/3.0/deed.en_US
     */

#if defined(_WIN32)
#    include <windows.h>
#    include <psapi.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#    include <sys/resource.h>
#    include <unistd.h>
#    if defined(__APPLE__) && defined(__MACH__)
#        include <mach/mach.h>
#    elif (defined(_AIX) || defined(__TOS__AIX__)) \
        || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#        include <fcntl.h>
#        include <procfs.h>
#    elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#        include <stdio.h>
#    endif
#endif

    /**
     * Returns the peak (maximum so far) resident set size (physical
     * memory use) measured in bytes, or zero if the value cannot be
     * determined on this OS.
     */
    size_t getPeakRSS()
    {
#if defined(_WIN32)
        /* Windows -------------------------------------------------- */
        PROCESS_MEMORY_COUNTERS info;
        GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
        return static_cast<size_t>(info.PeakWorkingSetSize);

#elif (defined(_AIX) || defined(__TOS__AIX__)) \
    || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
        /* AIX and Solaris ------------------------------------------ */
        struct psinfo psinfo;
        int fd = -1;
        if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1) return static_cast<size_t>(0); /* Can't open? */
        if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
        {
            close(fd);
            return static_cast<size_t>(0); /* Can't read? */
        }
        close(fd);
        return static_cast<size_t>(psinfo.pr_rssize) * static_cast<size_t>(1024);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
        /* BSD, Linux, and OSX -------------------------------------- */
        struct rusage rusage;
        getrusage(RUSAGE_SELF, &rusage);
#    if defined(__APPLE__) && defined(__MACH__)
        return static_cast<size_t>(rusage.ru_maxrss);
#    else
        return static_cast<size_t>(rusage.ru_maxrss) * static_cast<size_t>(1024);
#    endif

#else
        /* Unknown OS ----------------------------------------------- */
        return static_cast<size_t>(0); /* Unsupported. */
#endif
    }

    /**
     * Returns the current resident set size (physical memory use) measured
     * in bytes, or zero if the value cannot be determined on this OS.
     */
    size_t getCurrentRSS()
    {
#if defined(_WIN32)
        /* Windows -------------------------------------------------- */
        PROCESS_MEMORY_COUNTERS info;
        GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
        return static_cast<size_t>(info.WorkingSetSize);

#elif defined(__APPLE__) && defined(__MACH__)
        /* OSX ------------------------------------------------------ */
        struct mach_task_basic_info info;
        mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
        if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &infoCount)
            != KERN_SUCCESS)
            return static_cast<size_t>(0); /* Can't access? */
        return static_cast<size_t>(info.resident_size);

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
        /* Linux ---------------------------------------------------- */
        long rss = 0L;
        FILE* fp = NULL;
        if ((fp = fopen("/proc/self/statm", "r")) == NULL) return static_cast<size_t>(0); /* Can't open? */
        if (fscanf(fp, "%*s%ld", &rss) != 1)
        {
            fclose(fp);
            return static_cast<size_t>(0); /* Can't read? */
        }
        fclose(fp);
        return static_cast<size_t>(rss) * static_cast<size_t>(sysconf(_SC_PAGESIZE));

#else
        /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
        return static_cast<size_t>(0); /* Unsupported. */
#endif
    }

}  // end anonymous namespace

////////////////////////////////////////////////////////////////////

size_t System::availableMemory()
{
    return getMemorySize();
}

////////////////////////////////////////////////////////////////////

size_t System::peakMemoryUsage()
{
    return getPeakRSS();
}

////////////////////////////////////////////////////////////////////

size_t System::currentMemoryUsage()
{
    return getCurrentRSS();
}

////////////////////////////////////////////////////////////////////
