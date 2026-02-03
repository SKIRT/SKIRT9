/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FluxRecorder.hpp"
#include "FITSInOut.hpp"
#include "Indices.hpp"
#include "LockFree.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // indices for detector arrays that need calibration
    enum {
        Total = 0,
        Transparent,
        PrimaryDirect,
        PrimaryScattered,
        SecondaryDirect,
        SecondaryScattered,
        SecondaryTransparent,
        TotalQ,
        TotalU,
        TotalV,
        TransparentQ,
        TransparentU,
        TransparentV,
        PrimaryDirectQ,
        PrimaryDirectU,
        PrimaryDirectV,
        PrimaryScatteredQ,
        PrimaryScatteredU,
        PrimaryScatteredV,
        SecondaryDirectQ,
        SecondaryDirectU,
        SecondaryDirectV,
        SecondaryScatteredQ,
        SecondaryScatteredU,
        SecondaryScatteredV,
        SecondaryTransparentQ,
        SecondaryTransparentU,
        SecondaryTransparentV,
        PrimaryScatteredLevel
    };

    // the highest contribution power to track, i.e. the largest k in sum(w_i**k)
    //  - include k=0, which tracks the number of detections
    //  - thus, the number of detector arrays for statistics is this number plus one
    //  - these detector arrays do not need calibration!
    const int maxContributionPower = 4;
}

////////////////////////////////////////////////////////////////////

FluxRecorder::FluxRecorder(const SimulationItem* parentItem) : _parentItem{parentItem} {}

////////////////////////////////////////////////////////////////////

void FluxRecorder::setSimulationInfo(string instrumentName, const WavelengthGrid* lambdagrid, bool hasMedium,
                                     bool hasMediumEmission)
{
    _instrumentName = instrumentName;
    _lambdagrid = lambdagrid;
    _hasMedium = hasMedium;
    _hasMediumEmission = hasMediumEmission;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::setUserFlags(bool recordComponents, int numScatteringLevels, bool recordPolarization,
                                bool recordStatistics)
{
    _recordComponents = recordComponents;
    _numScatteringLevels = numScatteringLevels;
    _recordPolarization = recordPolarization;
    _recordStatistics = recordStatistics;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::setObserverAngles(double inclination, double azimuth, double roll)
{
    _inclination = inclination;
    _azimuth = azimuth;
    _roll = roll;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::setRestFrameDistance(double distance)
{
    _redshift = 0.;
    _angularDiameterDistance = distance;
    _luminosityDistance = distance;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::setObserverFrameRedshift(double redshift, double angularDiameterDistance, double luminosityDistance)
{
    _redshift = redshift;
    _angularDiameterDistance = angularDiameterDistance;
    _luminosityDistance = luminosityDistance;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::includeFluxDensityForDistant()
{
    _includeFluxDensity = true;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::includeSurfaceBrightnessForDistant(int numPixelsX, int numPixelsY, double pixelSizeX,
                                                      double pixelSizeY, double centerX, double centerY)
{
    _includeSurfaceBrightness = true;
    _numPixelsX = numPixelsX;
    _numPixelsY = numPixelsY;
    _pixelSizeX = pixelSizeX;
    _pixelSizeY = pixelSizeY;
    _centerX = centerX;
    _centerY = centerY;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::includeSurfaceBrightnessForLocal(int numPixelsX, int numPixelsY, double solidAnglePerPixel,
                                                    double incrementX, double incrementY, double centerX,
                                                    double centerY, string quantityXY)
{
    _includeSurfaceBrightness = true;
    _local = true;
    _numPixelsX = numPixelsX;
    _numPixelsY = numPixelsY;
    _solidAnglePerPixel = solidAnglePerPixel;
    _incrementX = incrementX;
    _incrementY = incrementY;
    _centerX = centerX;
    _centerY = centerY;
    _quantityXY = quantityXY;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::finalizeConfiguration()
{
    // get a pointer to the medium system, if present
    _ms = _parentItem->find<MediumSystem>(false);

    // get array lengths
    _numPixelsInFrame = _numPixelsX * _numPixelsY;  // convert to size_t before calculating lenIFU
    size_t lenSED = _includeFluxDensity ? _lambdagrid->numBins() : 0;
    size_t lenIFU = _includeSurfaceBrightness ? _numPixelsInFrame * _lambdagrid->numBins() : 0;

    // do not try to record components if there is no medium
    _recordTotalOnly = !_recordComponents || !_hasMedium;

    // allocate the appropriate number of flux detector arrays
    _sed.resize(PrimaryScatteredLevel + _numScatteringLevels);
    _ifu.resize(PrimaryScatteredLevel + _numScatteringLevels);

    // resize the flux detector arrays according to the configuration
    if (_recordTotalOnly)
    {
        _sed[Total].resize(lenSED);
        _ifu[Total].resize(lenIFU);
    }
    else
    {
        _sed[Transparent].resize(lenSED);
        _ifu[Transparent].resize(lenIFU);
        _sed[PrimaryDirect].resize(lenSED);
        _ifu[PrimaryDirect].resize(lenIFU);
        _sed[PrimaryScattered].resize(lenSED);
        _ifu[PrimaryScattered].resize(lenIFU);

        for (int i = 0; i != _numScatteringLevels; ++i)
        {
            _sed[PrimaryScatteredLevel + i].resize(lenSED);
            _ifu[PrimaryScatteredLevel + i].resize(lenIFU);
        }
        if (_hasMediumEmission)
        {
            _sed[SecondaryTransparent].resize(lenSED);
            _ifu[SecondaryTransparent].resize(lenIFU);
            _sed[SecondaryDirect].resize(lenSED);
            _ifu[SecondaryDirect].resize(lenIFU);
            _sed[SecondaryScattered].resize(lenSED);
            _ifu[SecondaryScattered].resize(lenIFU);
        }
    }
    if (_recordPolarization)
    {
        _sed[TotalQ].resize(lenSED);
        _ifu[TotalQ].resize(lenIFU);
        _sed[TotalU].resize(lenSED);
        _ifu[TotalU].resize(lenIFU);
        _sed[TotalV].resize(lenSED);
        _ifu[TotalV].resize(lenIFU);

        if (!_recordTotalOnly)
        {
            _sed[TransparentQ].resize(lenSED);
            _sed[TransparentU].resize(lenSED);
            _sed[TransparentV].resize(lenSED);

            _sed[PrimaryDirectQ].resize(lenSED);
            _sed[PrimaryDirectU].resize(lenSED);
            _sed[PrimaryDirectV].resize(lenSED);

            _sed[PrimaryScatteredQ].resize(lenSED);
            _sed[PrimaryScatteredU].resize(lenSED);
            _sed[PrimaryScatteredV].resize(lenSED);

            if (_hasMediumEmission)
            {
                _sed[SecondaryDirectQ].resize(lenSED);
                _sed[SecondaryDirectU].resize(lenSED);
                _sed[SecondaryDirectV].resize(lenSED);

                _sed[SecondaryScatteredQ].resize(lenSED);
                _sed[SecondaryScatteredU].resize(lenSED);
                _sed[SecondaryScatteredV].resize(lenSED);

                _sed[SecondaryTransparentQ].resize(lenSED);
                _sed[SecondaryTransparentU].resize(lenSED);
                _sed[SecondaryTransparentV].resize(lenSED);
            }
        }
    }

    // allocate and resize the statistics detector arrays
    if (_recordStatistics)
    {
        _wsed.resize(maxContributionPower + 1);
        _wifu.resize(maxContributionPower + 1);
        for (auto& array : _wsed) array.resize(lenSED);
        for (auto& array : _wifu) array.resize(lenIFU);
    }

    // calculate and log allocated memory size
    size_t allocatedSize = 0;
    for (const auto& array : _sed) allocatedSize += array.size();
    for (const auto& array : _ifu) allocatedSize += array.size();
    for (const auto& array : _wsed) allocatedSize += array.size();
    for (const auto& array : _wifu) allocatedSize += array.size();
    _parentItem->find<Log>()->info(_parentItem->typeAndName() + " allocated "
                                   + StringUtils::toMemSizeString(allocatedSize * sizeof(double)) + " of memory");
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::detect(PhotonPacket* pp, int l, double distance)
{
    // abort if we're not recording integrated fluxes and the photon packet arrives outside of the frame
    if (!_includeFluxDensity && l < 0) return;

    // get the photon packet's redshifted wavelength
    double wavelength = pp->wavelength() * (1. + _redshift);

    // get the wavelength bin indices that overlap the photon packet wavelength and perform recording for each
    for (int ell : _lambdagrid->bins(wavelength))
    {
        // get the luminosity contribution from the photon packet,
        // taking into account the transmission for the detector bin at this wavelength
        double L = pp->luminosity() * _lambdagrid->transmission(ell, wavelength);

        // adjust the luminosity for near distance if needed
        if (_local)
        {
            L /= distance * distance;
        }

        // apply the extinction along the path to the recorder
        double Lext = L;
        if (_hasMedium)
        {
            // if this photon packet has already been launched towards an instrument with the same observer type,
            // position and viewing direction, simply recover the stored optical depth from the photon packet;
            // otherwise calculate the optical depth and store it in the photon packet for the next instrument
            double tau;
            if (pp->hasObservedOpticalDepth())
            {
                tau = pp->observedOpticalDepth();
            }
            else
            {
                tau = _ms->getExtinctionOpticalDepth(pp, distance);
                pp->setObservedOpticalDepth(tau);
            }
            Lext *= exp(-tau);
        }

        // get number of scatterings (because we use it a lot)
        int numScatt = pp->numScatt();

        // record in SED arrays
        if (_includeFluxDensity)
        {
            if (_recordTotalOnly)
            {
                LockFree::add(_sed[Total][ell], Lext);
            }
            else
            {
                if (pp->hasPrimaryOrigin())
                {
                    if (numScatt == 0)
                    {
                        LockFree::add(_sed[Transparent][ell], L);
                        LockFree::add(_sed[PrimaryDirect][ell], Lext);
                    }
                    else
                    {
                        LockFree::add(_sed[PrimaryScattered][ell], Lext);
                        if (numScatt <= _numScatteringLevels)
                            LockFree::add(_sed[PrimaryScatteredLevel + numScatt - 1][ell], Lext);
                    }
                }
                else
                {
                    if (numScatt == 0)
                    {
                        LockFree::add(_sed[SecondaryTransparent][ell], L);
                        LockFree::add(_sed[SecondaryDirect][ell], Lext);
                    }
                    else
                    {
                        LockFree::add(_sed[SecondaryScattered][ell], Lext);
                    }
                }
            }
            if (_recordPolarization)
            {
                LockFree::add(_sed[TotalQ][ell], Lext * pp->stokesQ());
                LockFree::add(_sed[TotalU][ell], Lext * pp->stokesU());
                LockFree::add(_sed[TotalV][ell], Lext * pp->stokesV());

                if (!_recordTotalOnly)
                {
                    if (pp->hasPrimaryOrigin())
                    {
                        if (numScatt == 0)
                        {
                            LockFree::add(_sed[TransparentQ][ell], L * pp->stokesQ());
                            LockFree::add(_sed[TransparentU][ell], L * pp->stokesU());
                            LockFree::add(_sed[TransparentV][ell], L * pp->stokesV());
                            LockFree::add(_sed[PrimaryDirectQ][ell], Lext * pp->stokesQ());
                            LockFree::add(_sed[PrimaryDirectU][ell], Lext * pp->stokesU());
                            LockFree::add(_sed[PrimaryDirectV][ell], Lext * pp->stokesV());
                        }
                        else
                        {
                            LockFree::add(_sed[PrimaryScatteredQ][ell], Lext * pp->stokesQ());
                            LockFree::add(_sed[PrimaryScatteredU][ell], Lext * pp->stokesU());
                            LockFree::add(_sed[PrimaryScatteredV][ell], Lext * pp->stokesV());
                        }
                    }
                    else
                    {
                        if (numScatt == 0)
                        {
                            LockFree::add(_sed[SecondaryDirectQ][ell], Lext * pp->stokesQ());
                            LockFree::add(_sed[SecondaryDirectU][ell], Lext * pp->stokesU());
                            LockFree::add(_sed[SecondaryDirectV][ell], Lext * pp->stokesV());
                            LockFree::add(_sed[SecondaryTransparentQ][ell], L * pp->stokesQ());
                            LockFree::add(_sed[SecondaryTransparentU][ell], L * pp->stokesU());
                            LockFree::add(_sed[SecondaryTransparentV][ell], L * pp->stokesV());
                        }
                        else
                        {
                            LockFree::add(_sed[SecondaryScatteredQ][ell], Lext * pp->stokesQ());
                            LockFree::add(_sed[SecondaryScatteredU][ell], Lext * pp->stokesU());
                            LockFree::add(_sed[SecondaryScatteredV][ell], Lext * pp->stokesV());
                        }
                    }
                }
            }
        }

        // record in IFU arrays
        if (_includeSurfaceBrightness && l >= 0)
        {
            size_t lell = l + ell * _numPixelsInFrame;

            if (_recordTotalOnly)
            {
                LockFree::add(_ifu[Total][lell], Lext);
            }
            else
            {
                if (pp->hasPrimaryOrigin())
                {
                    if (numScatt == 0)
                    {
                        LockFree::add(_ifu[Transparent][lell], L);
                        LockFree::add(_ifu[PrimaryDirect][lell], Lext);
                    }
                    else
                    {
                        LockFree::add(_ifu[PrimaryScattered][lell], Lext);
                        if (numScatt <= _numScatteringLevels)
                            LockFree::add(_ifu[PrimaryScatteredLevel + numScatt - 1][lell], Lext);
                    }
                }
                else
                {
                    if (numScatt == 0)
                    {
                        LockFree::add(_ifu[SecondaryTransparent][lell], L);
                        LockFree::add(_ifu[SecondaryDirect][lell], Lext);
                    }
                    else
                    {
                        LockFree::add(_ifu[SecondaryScattered][lell], Lext);
                    }
                }
            }
            if (_recordPolarization)
            {
                LockFree::add(_ifu[TotalQ][lell], Lext * pp->stokesQ());
                LockFree::add(_ifu[TotalU][lell], Lext * pp->stokesU());
                LockFree::add(_ifu[TotalV][lell], Lext * pp->stokesV());
            }
        }

        // record statistics for both SEDs and IFUs
        if (_recordStatistics)
        {
            ContributionList* contributionList = _contributionLists.local();
            if (!contributionList->hasHistoryIndex(pp->historyIndex()))
            {
                recordContributions(contributionList);
                contributionList->reset(pp->historyIndex());
            }
            contributionList->addContribution(ell, l, Lext);
        }
    }
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::flush()
{
    // record the dangling contributions from all threads
    for (ContributionList* contributionList : _contributionLists.all())
    {
        recordContributions(contributionList);
        contributionList->reset();
    }
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::calibrateAndWrite()
{
    // collect recorded data from all processes
    for (auto& array : _sed) ProcessManager::sumToRoot(array);
    for (auto& array : _ifu) ProcessManager::sumToRoot(array);
    for (auto& array : _wsed) ProcessManager::sumToRoot(array);
    for (auto& array : _wifu) ProcessManager::sumToRoot(array);

    // calibrate and write only in the root process
    if (!ProcessManager::isRoot()) return;

    // calculate front factors for converting from recorded quantities to output quantities
    // (for local instruments, the distance correction already happened)
    double fourpid2 = 4. * M_PI * (_local ? 1. : _luminosityDistance * _luminosityDistance);
    double omega = _local ? _solidAnglePerPixel
                          : 4. * atan(0.5 * _pixelSizeX / _angularDiameterDistance)
                                * atan(0.5 * _pixelSizeY / _angularDiameterDistance);

    // convert from recorded quantities to output quantities and from internal units to user-selected output units
    // (for performance reasons, determine the units scaling factor only once for each wavelength)
    Units* units = _parentItem->find<Units>();
    int numWavelengths = _lambdagrid->numBins();
    for (int ell = 0; ell != numWavelengths; ++ell)
    {
        // SEDs
        if (_includeFluxDensity)
        {
            double factor = 1. / fourpid2 / _lambdagrid->effectiveWidth(ell)
                            * units->ofluxdensity(_lambdagrid->wavelength(ell), 1.);
            for (auto& array : _sed)
                if (array.size()) array[ell] *= factor;
        }
        // IFUs
        if (_includeSurfaceBrightness)
        {
            double factor = 1. / fourpid2 / omega / _lambdagrid->effectiveWidth(ell)
                            * units->osurfacebrightness(_lambdagrid->wavelength(ell), 1.);
            size_t begin = ell * _numPixelsInFrame;
            size_t end = begin + _numPixelsInFrame;
            for (auto& array : _ifu)
                if (array.size())
                    for (size_t lell = begin; lell != end; ++lell) array[lell] *= factor;
        }
    }

    // write SEDs to a single text file (with multiple columns)
    if (_includeFluxDensity)
    {
        // Build a list of column names and corresponding pointers to sed arrays (which may be empty)
        vector<string> sedNames;
        vector<Array*> sedArrays;

        // add the total flux; if we didn't record it directly, calculate it now
        sedNames.push_back("total flux");
        Array sedTotal;
        if (_recordTotalOnly)
            sedArrays.push_back(&_sed[Total]);
        else
        {
            sedTotal = _sed[PrimaryDirect] + _sed[PrimaryScattered];
            if (_hasMediumEmission) sedTotal += _sed[SecondaryDirect] + _sed[SecondaryScattered];
            sedArrays.push_back(&sedTotal);
        }

        // add the flux components, if requested
        // we always add all of them, even if some of them are zero
        if (_recordComponents)
        {
            // add transparent flux
            // if we did not actually record components (because there are no media), use the total flux instead
            sedNames.push_back("transparent flux");
            sedArrays.push_back(_recordTotalOnly ? &_sed[Total] : &_sed[Transparent]);

            // add the actual components of the total flux
            sedNames.insert(sedNames.end(), {"direct primary flux", "scattered primary flux", "direct secondary flux",
                                             "scattered secondary flux", "transparent secondary flux"});
            sedArrays.insert(sedArrays.end(), {&_sed[PrimaryDirect], &_sed[PrimaryScattered], &_sed[SecondaryDirect],
                                               &_sed[SecondaryScattered], &_sed[SecondaryTransparent]});
        }

        // add the polarization components, if requested
        if (_recordPolarization)
        {
            sedNames.insert(sedNames.end(), {"total Stokes Q", "total Stokes U", "total Stokes V"});
            sedArrays.insert(sedArrays.end(), {&_sed[TotalQ], &_sed[TotalU], &_sed[TotalV]});

            if (_recordComponents && !_recordTotalOnly)
            {
                sedNames.insert(
                    sedNames.end(),
                    {"transparent Stokes Q", "transparent Stokes U", "transparent Stokes V", "direct primary Stokes Q",
                     "direct primary Stokes U", "direct primary Stokes V", "scattered primary Stokes Q",
                     "scattered primary Stokes U", "scattered primary Stokes V", "direct secondary Stokes Q",
                     "direct secondary Stokes U", "direct secondary Stokes V", "scattered secondary Stokes Q",
                     "scattered secondary Stokes U", "scattered secondary Stokes V", "transparent secondary Stokes Q",
                     "transparent secondary Stokes U", "transparent secondary Stokes V"});
                sedArrays.insert(sedArrays.end(),
                                 {&_sed[TransparentQ], &_sed[TransparentU], &_sed[TransparentV], &_sed[PrimaryDirectQ],
                                  &_sed[PrimaryDirectU], &_sed[PrimaryDirectV], &_sed[PrimaryScatteredQ],
                                  &_sed[PrimaryScatteredU], &_sed[PrimaryScatteredV], &_sed[SecondaryDirectQ],
                                  &_sed[SecondaryDirectU], &_sed[SecondaryDirectV], &_sed[SecondaryScatteredQ],
                                  &_sed[SecondaryScatteredU], &_sed[SecondaryScatteredV], &_sed[SecondaryTransparentQ],
                                  &_sed[SecondaryTransparentU], &_sed[SecondaryTransparentV]});
            }
        }

        // add the scattering levels, if requested, even if they are all zero
        Array empty;
        if (_recordComponents)
            for (int i = 0; i != _numScatteringLevels; ++i)
            {
                sedNames.push_back(std::to_string(i + 1) + "-times scattered primary flux");
                sedArrays.push_back(_recordTotalOnly ? &empty : &_sed[PrimaryScatteredLevel + i]);
            }

        // construct header comment line
        string header = "# SED at ";
        header += "inclination " + StringUtils::toString(units->oposangle(_inclination)) + " " + units->uposangle();
        header += ", azimuth " + StringUtils::toString(units->oposangle(_azimuth)) + " " + units->uposangle();
        if (_recordPolarization)
        {
            header += ", roll " + StringUtils::toString(units->oposangle(_roll)) + " " + units->uposangle();
        }
        if (_redshift)
        {
            header += ", redshift " + StringUtils::toString(_redshift);
            header += ", luminosity distance " + StringUtils::toString(units->odistance(_luminosityDistance)) + " "
                      + units->udistance();
        }
        else
        {
            header +=
                ", distance " + StringUtils::toString(units->odistance(_luminosityDistance)) + " " + units->udistance();
        }

        // open the file and add the column headers
        TextOutFile sedFile(_parentItem, _instrumentName + "_sed", "SED");
        sedFile.writeLine(header);
        sedFile.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
        for (const string& name : sedNames)
        {
            sedFile.addColumn(name + "; " + units->sfluxdensity(), units->ufluxdensity());
        }

        // write the column data
        for (int ell : Indices(numWavelengths, units->rwavelength()))
        {
            vector<double> values({units->owavelength(_lambdagrid->wavelength(ell))});
            for (const Array* array : sedArrays) values.push_back(array->size() ? (*array)[ell] : 0.);
            sedFile.writeRow(values);
        }
        sedFile.close();

        // output statistics to a seperate file
        if (_recordStatistics)
        {
            // open the file and add the column headers
            TextOutFile statFile(_parentItem, _instrumentName + "_sedstats", "SED statistics");
            statFile.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
            for (int k = 0; k <= maxContributionPower; ++k)
            {
                statFile.addColumn("Sum[w_i**" + std::to_string(k) + "]");
            }
            statFile.writeLine("# --> w_i is luminosity contribution (in W) from i_th launched photon");

            // write the column data
            for (int ell : Indices(numWavelengths, units->rwavelength()))
            {
                vector<double> values({units->owavelength(_lambdagrid->wavelength(ell))});
                for (int k = 0; k <= maxContributionPower; ++k) values.push_back(_wsed[k][ell]);
                statFile.writeRow(values);
            }
            statFile.close();
        }
    }

    // write IFUs to FITS files (one file per IFU)
    if (_includeSurfaceBrightness)
    {
        // Build a list of file names and corresponding pointers to ifu arrays (which may be empty)
        vector<string> ifuNames;
        vector<Array*> ifuArrays;

        // add the total flux; if we didn't record it directly, calculate it now
        ifuNames.push_back("total");
        Array ifuTotal;
        if (_recordTotalOnly)
            ifuArrays.push_back(&_ifu[Total]);
        else
        {
            ifuTotal = _ifu[PrimaryDirect] + _ifu[PrimaryScattered];
            if (_hasMediumEmission) ifuTotal += _ifu[SecondaryDirect] + _ifu[SecondaryScattered];
            ifuArrays.push_back(&ifuTotal);
        }

        // add the flux components, if requested
        if (_recordComponents)
        {
            // add the transparent flux only if it may differ from the total flux
            if (!_recordTotalOnly)
            {
                ifuNames.push_back("transparent");
                ifuArrays.push_back(&_ifu[Transparent]);
            }
            // add the actual components of the total flux (empty arrays will be ignored later on)
            ifuNames.insert(ifuNames.end(), {"primarydirect", "primaryscattered", "secondarytransparent",
                                             "secondarydirect", "secondaryscattered"});
            ifuArrays.insert(ifuArrays.end(),
                             {&_ifu[PrimaryDirect], &_ifu[PrimaryScattered], &_ifu[SecondaryTransparent],
                              &_ifu[SecondaryDirect], &_ifu[SecondaryScattered]});
        }

        // add the polarization components, if requested
        if (_recordPolarization)
        {
            ifuNames.insert(ifuNames.end(), {"stokesQ", "stokesU", "stokesV"});
            ifuArrays.insert(ifuArrays.end(), {&_ifu[TotalQ], &_ifu[TotalU], &_ifu[TotalV]});
        }

        // add the scattering levels, if requested
        if (!_recordTotalOnly)
            for (int i = 0; i != _numScatteringLevels; ++i)
            {
                ifuNames.push_back("primaryscattered" + std::to_string(i + 1));
                ifuArrays.push_back(&_ifu[PrimaryScatteredLevel + i]);
            }

        // copy the wavelength grid in output units
        Array wavegrid(numWavelengths);
        for (int ell = 0; ell != numWavelengths; ++ell)
            wavegrid[ell] = units->owavelength(_lambdagrid->wavelength(ell));

        // reverse the ordering of the wavelength grid and frames if necessary
        if (units->rwavelength())
        {
            NR::reverse(wavegrid);

            // flux frames
            for (auto array : ifuArrays)
                if (array->size()) NR::reverse(*array, _numPixelsInFrame);

            // statistics frames
            for (auto& array : _wifu)
                if (array.size()) NR::reverse(array, _numPixelsInFrame);
        }

        // determine spatial axes values and units
        double incx, incy, cx, cy;
        string unitsxy;
        if (_local)
        {
            // for local instruments, use the metadata provided by the instrument
            if (_quantityXY.empty())
            {
                incx = _incrementX;
                incy = _incrementY;
                cx = _centerX;
                cy = _centerY;
                unitsxy = "1";
            }
            else
            {
                incx = units->out(_quantityXY, _incrementX);
                incy = units->out(_quantityXY, _incrementY);
                cx = units->out(_quantityXY, _centerX);
                cy = units->out(_quantityXY, _centerY);
                unitsxy = units->unit(_quantityXY);
            }
        }
        else
        {
            // for distant instruments, convert to angular sizes
            incx = units->oangle(2. * atan(0.5 * _pixelSizeX / _angularDiameterDistance));
            incy = units->oangle(2. * atan(0.5 * _pixelSizeY / _angularDiameterDistance));
            cx = units->oangle(2. * atan(0.5 * _centerX / _angularDiameterDistance));
            cy = units->oangle(2. * atan(0.5 * _centerY / _angularDiameterDistance));
            unitsxy = units->uangle();
        }

        // determine observer info for distant instruments
        std::unique_ptr<FITSInOut::ObserverInfo> obsInfo;
        if (!_local)
        {
            obsInfo = std::make_unique<FITSInOut::ObserverInfo>();
            obsInfo->inclination = _inclination * (180. / M_PI);
            obsInfo->azimuth = _azimuth * (180. / M_PI);
            obsInfo->roll = _roll * (180. / M_PI);
            obsInfo->redshift = _redshift;
            obsInfo->luminosityDistance = units->odistance(_luminosityDistance);
            obsInfo->angularDiameterDistance = units->odistance(_angularDiameterDistance);
            obsInfo->distanceUnits = units->udistance();
        }

        // output the files (ignoring empty arrays)
        int numFiles = ifuNames.size();
        for (int q = 0; q != numFiles; ++q)
            if (ifuArrays[q]->size())
            {
                string filename = _instrumentName + "_" + ifuNames[q];
                string description = ifuNames[q] + " flux";
                FITSInOut::write(_parentItem, description, filename, *(ifuArrays[q]), units->usurfacebrightness(),
                                 _numPixelsX, _numPixelsY, incx, incy, cx, cy, unitsxy, wavegrid, units->uwavelength(),
                                 obsInfo.get());
            }

        // output statistics to additional files
        if (_recordStatistics)
        {
            // the output files have single-precision floating point numbers with range of only about 10^+-38
            // --> scale the values to a range that has a maximum of 10^+-38 to minimize the number of underflows
            const double WMAX = 1e38;
            Array cs(maxContributionPower);
            for (int k = 1; k <= maxContributionPower; ++k)
            {
                cs[k - 1] = pow(WMAX / _wifu[k].max(), 1. / k);  // inverse of WMAX == c**k w[k].max()
            }
            double c = cs.min();
            double cn = 1.;
            for (int k = 0; k <= maxContributionPower; ++k)
            {
                string filename = _instrumentName + "_stats" + std::to_string(k);
                string description = "sum of contributions to the power of " + std::to_string(k);
                _wifu[k] *= cn;
                FITSInOut::write(_parentItem, description, filename, _wifu[k], "", _numPixelsX, _numPixelsY, incx, incy,
                                 cx, cy, unitsxy, wavegrid, units->uwavelength());
                cn *= c;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::recordContributions(ContributionList* contributionList)
{
    // sort the contributions on wavelength and pixel index so that contributions to the same bin are consecutive
    contributionList->sort();
    const vector<Contribution>& contributions = contributionList->contributions();
    size_t numContributions = contributions.size();

    // for SEDs, group contributions on ell index (wavelength bin)
    if (_includeFluxDensity)
    {
        double w = 0;
        for (size_t i = 0; i != numContributions; ++i)
        {
            w += contributions[i].w();
            if (i + 1 == numContributions || contributions[i].ell() != contributions[i + 1].ell())
            {
                int ell = contributions[i].ell();
                double wn = 1.;
                for (int k = 0; k <= maxContributionPower; ++k)
                {
                    LockFree::add(_wsed[k][ell], wn);
                    wn *= w;
                }
                w = 0;
            }
        }
    }

    // for IFUs, group contributions on lell index (wavelength and pixel bins)
    if (_includeSurfaceBrightness)
    {
        double w = 0;
        for (size_t i = 0; i != numContributions; ++i)
        {
            w += contributions[i].w();
            if (i + 1 == numContributions || contributions[i].ell() != contributions[i + 1].ell()
                || contributions[i].l() != contributions[i + 1].l())
            {
                if (contributions[i].l() >= 0)
                {
                    size_t lell = contributions[i].l() + contributions[i].ell() * _numPixelsInFrame;
                    double wn = 1.;
                    for (int k = 0; k <= maxContributionPower; ++k)
                    {
                        LockFree::add(_wifu[k][lell], wn);
                        wn *= w;
                    }
                }
                w = 0;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
