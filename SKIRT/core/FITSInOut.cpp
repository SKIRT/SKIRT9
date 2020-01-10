/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FITSInOut.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "ProcessManager.hpp"
#include "System.hpp"
#include "fitsio.h"
#include <mutex>

////////////////////////////////////////////////////////////////////

void FITSInOut::read(const SimulationItem* item, string filename, Array& data, int& nx, int& ny, int& nz)
{
    // Determine the path of the input FITS file
    string filepath = item->find<FilePaths>()->input(filename);

    // Read the input data
    FITSInOut::read(filepath, data, nx, ny, nz);

    // Log the file path and the data dimensions
    Log* log = item->find<Log>();
    log->info(item->typeAndName() + " read FITS file " + filepath);
    if (nz > 1) log->info("  Data cube has " + std::to_string(nz) + " frames");
    log->info("  Frame dimensions are " + std::to_string(nx) + " x " + std::to_string(ny));
}

////////////////////////////////////////////////////////////////////

void FITSInOut::write(const SimulationItem* item, string description, string filename, const Array& data,
                      string dataUnits, int nx, int ny, double incx, double incy, double xc, double yc, string xyUnits,
                      const Array& z, string zUnits)
{
    // Only write the FITS file if this process is the root
    if (ProcessManager::isRoot())
    {
        // Determine the path of the output FITS file
        string filepath = item->find<FilePaths>()->output(filename + ".fits");

        // Write the FITS file
        FITSInOut::write(filepath, data, dataUnits, nx, ny, incx, incy, xc, yc, xyUnits, z, zUnits);

        // Log the file path
        item->find<Log>()->info(item->typeAndName() + " wrote " + description + " to FITS file " + filepath);
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // mutex to guard the FITS input/output operations
    std::mutex _mutex;

    // function to report cfitsio errors
    void report_error(string filepath, string action, int status)
    {
        char message[FLEN_STATUS];
        ffgerr(status, message);
        throw FATALERROR("Error while " + action + " FITS file " + filepath + "\n" + string(message));
    }
}

////////////////////////////////////////////////////////////////////

void FITSInOut::read(string filepath, Array& data, int& nx, int& ny, int& nz)
{
    // Acquire a global lock since the cfitsio library is not guaranteed to be reentrant
    std::unique_lock<std::mutex> lock(_mutex);

    // Open the FITS file
    int status = 0;
    fitsfile* fptr;
    ffdopn(&fptr, filepath.c_str(), READONLY, &status);
    if (status) report_error(filepath, "opening", status);

    // Get the dimensions of the data
    int naxis;
    long naxes[3];
    ffgidm(fptr, &naxis, &status);
    ffgisz(fptr, 3, naxes, &status);
    if (status) report_error(filepath, "reading", status);
    nx = naxis > 0 ? naxes[0] : 1;
    ny = naxis > 1 ? naxes[1] : 1;
    nz = naxis > 2 ? naxes[2] : 1;

    // Resize the data container
    size_t nelements = static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz);
    data.resize(nelements);

    // Read the array of pixels from the file
    int dummy;
    ffgpvd(fptr, 0, 1, nelements, 0, &data[0], &dummy, &status);
    if (status) report_error(filepath, "reading", status);

    // Close the file
    ffclos(fptr, &status);
    if (status) report_error(filepath, "reading", status);
}

////////////////////////////////////////////////////////////////////

void FITSInOut::write(string filepath, const Array& data, string dataUnits, int nx, int ny, double incx, double incy,
                      double xc, double yc, string xyUnits, const Array& z, string zUnits)
{
    // Get the z-axis size
    //   0:  a single frame that is not part of a datacube
    //   1:  a datacube with just a single frame, i.e. there is a z-axis with a single grid point
    //  >1:  a datacube with multiple frames, i.e. there is a z-axis with multiple grid points
    int nz = z.size();

    // Verify the data size
    size_t nelements = data.size();
    if (nelements != static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz ? nz : 1))
        throw FATALERROR("Inconsistent data size when creating FITS file " + filepath);
    long naxes[3] = {nx, ny, nz};

    // Acquire a global lock since the cfitsio library is not guaranteed to be reentrant
    // (only when it is built with ./configure --enable-reentrant; make)
    std::unique_lock<std::mutex> lock(_mutex);

    // Generate time stamp
    string stamp = System::timestamp(true);
    stamp.erase(19);  // remove milliseconds

    // Remove any existing file with the same name
    remove(filepath.c_str());

    // Create the fits file
    int status = 0;
    fitsfile* fptr;
    ffdkinit(&fptr, filepath.c_str(), &status);
    if (status) report_error(filepath, "creating", status);

    // Create the primary image (32-bit floating point pixels)
    ffcrim(fptr, FLOAT_IMG, (nz ? 3 : 2), naxes, &status);
    if (status) report_error(filepath, "creating", status);

    // Add the relevant keywords
    ffpkyg(fptr, "BSCALE", 1., 0, "Array value scale", &status);
    ffpkyg(fptr, "BZERO", 0., 0, "Array value offset", &status);
    ffpkys(fptr, "DATE", const_cast<char*>(stamp.c_str()), "Date and time of creation (UTC)", &status);
    ffpkys(fptr, "ORIGIN", const_cast<char*>("SKIRT simulation"), "Astronomical Observatory, Ghent University",
           &status);
    ffpkys(fptr, "BUNIT", const_cast<char*>(dataUnits.c_str()), "Physical unit of the array values", &status);
    ffpkyd(fptr, "CRPIX1", (nx + 1) / 2., 9, "X-axis coordinate system reference pixel", &status);
    ffpkyd(fptr, "CRVAL1", xc, 9, "Coordinate value at X-axis reference pixel", &status);
    ffpkyd(fptr, "CDELT1", incx, 9, "Coordinate increment along X-axis", &status);
    ffpkys(fptr, "CUNIT1", const_cast<char*>(xyUnits.c_str()), "Physical units of the X-axis", &status);
    ffpkys(fptr, "CTYPE1", " ", "Linear X coordinates", &status);
    ffpkyd(fptr, "CRPIX2", (ny + 1) / 2., 9, "Y-axis coordinate system reference pixel", &status);
    ffpkyd(fptr, "CRVAL2", yc, 9, "Coordinate value at Y-axis reference pixel", &status);
    ffpkyd(fptr, "CDELT2", incy, 9, "Coordinate increment along Y-axis", &status);
    ffpkys(fptr, "CUNIT2", const_cast<char*>(xyUnits.c_str()), "Physical units of the Y-axis", &status);
    ffpkys(fptr, "CTYPE2", " ", "Linear Y coordinates", &status);
    if (nz) ffpkys(fptr, "CUNIT3", const_cast<char*>(zUnits.c_str()), "Physical units of the Z-axis", &status);
    if (status) report_error(filepath, "writing", status);

    // Write the array of pixels to the image
    ffpprd(fptr, 0, 1, nelements, const_cast<double*>(&data[0]), &status);
    if (status) report_error(filepath, "writing", status);

    // If the data has 3 dimensions, write a FITS table extension with the values of the third axis
    if (nz)
    {
        // Create the table
        char* ttypev[] = {const_cast<char*>("GRID_POINTS")};
        char* tformv[] = {const_cast<char*>("E16.9")};
        char* tunitv[] = {const_cast<char*>(zUnits.c_str())};
        ffcrtb(fptr, ASCII_TBL, 0, 1, ttypev, tformv, tunitv, "Z-axis coordinate values", &status);
        if (status) report_error(filepath, "writing", status);

        // Write the single column
        void* gridpoints = const_cast<void*>(static_cast<const void*>(begin(z)));
        ffpcl(fptr, TDOUBLE, 1, 1, 1, nz, gridpoints, &status);
        if (status) report_error(filepath, "writing", status);
    }

    // Close the file
    ffclos(fptr, &status);
    if (status) report_error(filepath, "writing", status);
}

////////////////////////////////////////////////////////////////////

void FITSInOut::readColumn(string filepath, Array& data, int n)
{
    // Acquire a global lock since the cfitsio library is not guaranteed to be reentrant
    std::unique_lock<std::mutex> lock(_mutex);

    // Open the FITS file
    int status = 0;
    fitsfile* fptr;
    fftopn(&fptr, filepath.c_str(), READONLY, &status);
    if (status) report_error(filepath, "opening", status);

    // Get the dimensions of the table
    int ncols = 0;
    ffgncl(fptr, &ncols, &status);
    if (status) report_error(filepath, "reading", status);
    long nrows = 0;
    ffgnrw(fptr, &nrows, &status);
    if (status) report_error(filepath, "reading", status);
    if (ncols <= 0 || nrows < n) throw FATALERROR("Not enough table data in FITS file " + filepath);

    // Resize the data container
    data.resize(n);

    // Read the array of values from the first column in the table
    int dummy;
    ffgcvd(fptr, 1, 1, 1, n, 0., &data[0], &dummy, &status);
    if (status) report_error(filepath, "reading", status);

    // Close the file
    ffclos(fptr, &status);
    if (status) report_error(filepath, "reading", status);
}

////////////////////////////////////////////////////////////////////
