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
#include "TimeGrid.hpp"
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

void FluxRecorder::setSimulationInfo(string instrumentName, bool hasMedium, bool hasMediumEmission)
{
    _instrumentName = instrumentName;
    _hasMedium = hasMedium;
    _hasMediumEmission = hasMediumEmission;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::setWavelengthGrid(const WavelengthGrid* lambdagrid)
{
    _lambdagrid = lambdagrid;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::setTimeGrid(const TimeGrid* timegrid)
{
    _timegrid = timegrid;
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

void FluxRecorder::includeLightCurve()
{
    _includeLightCurve = true;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::includeSpectralTimeMap()
{
    _includeSpectralTimeMap = true;
}

////////////////////////////////////////////////////////////////////

void FluxRecorder::finalizeConfiguration()
{
    // get a pointer to the medium system, if present
    _ms = _parentItem->find<MediumSystem>(false);

    // get array lengths
    _numWavelengths = _lambdagrid->numBins();
    _numPixelsInFrame = _numPixelsX * _numPixelsY;  // convert to size_t before calculating lenIFU
    size_t lenSED = _includeFluxDensity ? _numWavelengths : 0;
    size_t lenIFU = _includeSurfaceBrightness ? _numPixelsInFrame * _numWavelengths : 0;
    size_t lenLC = _includeLightCurve ? _timegrid->numBins() : 0;
    size_t lenSTM = _includeSpectralTimeMap ? _numWavelengths * _timegrid->numBins() : 0;

    // do not try to record components if there is no medium
    _recordTotalOnly = !_recordComponents || !_hasMedium;

    // allocate the appropriate number of flux detector arrays
    _sed.resize(PrimaryScatteredLevel + _numScatteringLevels);
    _ifu.resize(PrimaryScatteredLevel + _numScatteringLevels);
    _lc.resize(PrimaryScatteredLevel + _numScatteringLevels);
    _lcw.resize(PrimaryScatteredLevel + _numScatteringLevels);
    _stm.resize(PrimaryScatteredLevel + _numScatteringLevels);

    // local function to resize flux detector arrays according to the configuration
    // params: vector of flux arrays; desired length of each array; include individual contributions for polarization
    auto resize = [this](vector<Array>& arrays, size_t len, bool polarcomp) {
        if (len)
        {
            if (_recordTotalOnly)
            {
                arrays[Total].resize(len);
            }
            else
            {
                arrays[Transparent].resize(len);
                arrays[PrimaryDirect].resize(len);
                arrays[PrimaryScattered].resize(len);

                for (int i = 0; i != _numScatteringLevels; ++i)
                {
                    arrays[PrimaryScatteredLevel + i].resize(len);
                }
                if (_hasMediumEmission)
                {
                    arrays[SecondaryTransparent].resize(len);
                    arrays[SecondaryDirect].resize(len);
                    arrays[SecondaryScattered].resize(len);
                }
            }
            if (_recordPolarization)
            {
                arrays[TotalQ].resize(len);
                arrays[TotalU].resize(len);
                arrays[TotalV].resize(len);

                if (polarcomp && !_recordTotalOnly)
                {
                    arrays[TransparentQ].resize(len);
                    arrays[TransparentU].resize(len);
                    arrays[TransparentV].resize(len);

                    arrays[PrimaryDirectQ].resize(len);
                    arrays[PrimaryDirectU].resize(len);
                    arrays[PrimaryDirectV].resize(len);

                    arrays[PrimaryScatteredQ].resize(len);
                    arrays[PrimaryScatteredU].resize(len);
                    arrays[PrimaryScatteredV].resize(len);

                    if (_hasMediumEmission)
                    {
                        arrays[SecondaryDirectQ].resize(len);
                        arrays[SecondaryDirectU].resize(len);
                        arrays[SecondaryDirectV].resize(len);

                        arrays[SecondaryScatteredQ].resize(len);
                        arrays[SecondaryScatteredU].resize(len);
                        arrays[SecondaryScatteredV].resize(len);

                        arrays[SecondaryTransparentQ].resize(len);
                        arrays[SecondaryTransparentU].resize(len);
                        arrays[SecondaryTransparentV].resize(len);
                    }
                }
            }
        }
    };

    // resize the flux detector arrays
    resize(_sed, lenSED, true);
    resize(_ifu, lenIFU, false);
    resize(_lc, lenLC, true);
    resize(_lcw, lenLC, true);
    resize(_stm, lenSTM, false);

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
    for (const auto& array : _lc) allocatedSize += array.size();
    for (const auto& array : _lcw) allocatedSize += array.size();
    for (const auto& array : _stm) allocatedSize += array.size();
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
        if (_local) L /= distance * distance;

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

        // local function to record the contribution in the flux detector arrays according to the configuration
        // params: vector of flux arrays; index in array, transparant luminosity, extincted luminosity,
        //         include individual contributions for polarization
        auto record = [this, pp](vector<Array>& arrays, size_t index, double L, double Lext, bool polarcomp) {
            int numScatt = pp->numScatt();

            if (_recordTotalOnly)
            {
                LockFree::add(arrays[Total][index], Lext);
            }
            else
            {
                if (pp->hasPrimaryOrigin())
                {
                    if (numScatt == 0)
                    {
                        LockFree::add(arrays[Transparent][index], L);
                        LockFree::add(arrays[PrimaryDirect][index], Lext);
                    }
                    else
                    {
                        LockFree::add(arrays[PrimaryScattered][index], Lext);
                        if (numScatt <= _numScatteringLevels)
                            LockFree::add(arrays[PrimaryScatteredLevel + numScatt - 1][index], Lext);
                    }
                }
                else
                {
                    if (numScatt == 0)
                    {
                        LockFree::add(arrays[SecondaryTransparent][index], L);
                        LockFree::add(arrays[SecondaryDirect][index], Lext);
                    }
                    else
                    {
                        LockFree::add(arrays[SecondaryScattered][index], Lext);
                    }
                }
            }
            if (_recordPolarization)
            {
                LockFree::add(arrays[TotalQ][index], Lext * pp->stokesQ());
                LockFree::add(arrays[TotalU][index], Lext * pp->stokesU());
                LockFree::add(arrays[TotalV][index], Lext * pp->stokesV());

                if (polarcomp && !_recordTotalOnly)
                {
                    if (pp->hasPrimaryOrigin())
                    {
                        if (numScatt == 0)
                        {
                            LockFree::add(arrays[TransparentQ][index], L * pp->stokesQ());
                            LockFree::add(arrays[TransparentU][index], L * pp->stokesU());
                            LockFree::add(arrays[TransparentV][index], L * pp->stokesV());
                            LockFree::add(arrays[PrimaryDirectQ][index], Lext * pp->stokesQ());
                            LockFree::add(arrays[PrimaryDirectU][index], Lext * pp->stokesU());
                            LockFree::add(arrays[PrimaryDirectV][index], Lext * pp->stokesV());
                        }
                        else
                        {
                            LockFree::add(arrays[PrimaryScatteredQ][index], Lext * pp->stokesQ());
                            LockFree::add(arrays[PrimaryScatteredU][index], Lext * pp->stokesU());
                            LockFree::add(arrays[PrimaryScatteredV][index], Lext * pp->stokesV());
                        }
                    }
                    else
                    {
                        if (numScatt == 0)
                        {
                            LockFree::add(arrays[SecondaryDirectQ][index], Lext * pp->stokesQ());
                            LockFree::add(arrays[SecondaryDirectU][index], Lext * pp->stokesU());
                            LockFree::add(arrays[SecondaryDirectV][index], Lext * pp->stokesV());
                            LockFree::add(arrays[SecondaryTransparentQ][index], L * pp->stokesQ());
                            LockFree::add(arrays[SecondaryTransparentU][index], L * pp->stokesU());
                            LockFree::add(arrays[SecondaryTransparentV][index], L * pp->stokesV());
                        }
                        else
                        {
                            LockFree::add(arrays[SecondaryScatteredQ][index], Lext * pp->stokesQ());
                            LockFree::add(arrays[SecondaryScatteredU][index], Lext * pp->stokesU());
                            LockFree::add(arrays[SecondaryScatteredV][index], Lext * pp->stokesV());
                        }
                    }
                }
            }
        };

        // record in SED arrays
        if (_includeFluxDensity) record(_sed, ell, L, Lext, true);

        // record in IFU arrays
        if (_includeSurfaceBrightness && l >= 0) record(_ifu, l + ell * _numPixelsInFrame, L, Lext, false);

        // if this is a time instrument
        if (_includeLightCurve || _includeSpectralTimeMap)
        {
            // get the time grid bin index corresponding to the distance travelled by the photon packet
            int k = _timegrid->binForDistance(pp->distance());
            if (k >= 0)
            {
                // record in LC arrays
                if (_includeLightCurve)
                {
                    // record both the plain contribution and the contribution multiplied by the wavelength
                    // to allow converting the aggregated value between an amount of energy and a number of photons
                    record(_lc, k, L, Lext, true);
                    record(_lcw, k, L * wavelength, Lext * wavelength, true);
                }

                // record in STM arrays
                if (_includeSpectralTimeMap) record(_stm, ell + k * _numWavelengths, L, Lext, false);
            }
        }

        // record statistics for both SEDs and IFUs (not yet implemented for time instruments)
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
    for (auto& array : _lc) ProcessManager::sumToRoot(array);
    for (auto& array : _lcw) ProcessManager::sumToRoot(array);
    for (auto& array : _stm) ProcessManager::sumToRoot(array);
    for (auto& array : _wsed) ProcessManager::sumToRoot(array);
    for (auto& array : _wifu) ProcessManager::sumToRoot(array);

    // calibrate and write only in the root process
    if (!ProcessManager::isRoot()) return;

    // get units
    Units* units = _parentItem->find<Units>();

    // calculate distance-related front factors for converting from recorded quantities to output quantities
    // (for local instruments, the distance correction already happened)
    double fourpid2 = 4. * M_PI * (_local ? 1. : _luminosityDistance * _luminosityDistance);
    double omega = _local ? _solidAnglePerPixel
                          : 4. * atan(0.5 * _pixelSizeX / _angularDiameterDistance)
                                * atan(0.5 * _pixelSizeY / _angularDiameterDistance);

    // ---------------------- local helper functions ----------------------

    // "global" arrays that can be inserted into the lists constructed by the build functions below
    Array empty, total;

    // local function to build synchronized lists of column names and recorded array pointers
    // params: vector of recorded arrays (in), list of column names (out), list of array pointers (out)
    auto buildCols = [this, &empty, &total](vector<Array>& arrays, vector<string>& colNames,
                                            vector<Array*>& arrayPtrs) {
        // add the total flux; if we didn't record it directly, calculate it now
        colNames.push_back("total flux");
        if (_recordTotalOnly)
            arrayPtrs.push_back(&arrays[Total]);
        else
        {
            total = arrays[PrimaryDirect] + arrays[PrimaryScattered];
            if (_hasMediumEmission) total += arrays[SecondaryDirect] + arrays[SecondaryScattered];
            arrayPtrs.push_back(&total);
        }

        // add the flux components, if requested
        // we always add all of them, even if some of them are zero
        if (_recordComponents)
        {
            // add transparent flux
            // if we did not actually record components (because there are no media), use the total flux instead
            colNames.push_back("transparent flux");
            arrayPtrs.push_back(_recordTotalOnly ? &arrays[Total] : &arrays[Transparent]);

            // add the actual components of the total flux
            colNames.insert(colNames.end(), {"direct primary flux", "scattered primary flux", "direct secondary flux",
                                             "scattered secondary flux", "transparent secondary flux"});
            arrayPtrs.insert(arrayPtrs.end(),
                             {&arrays[PrimaryDirect], &arrays[PrimaryScattered], &arrays[SecondaryDirect],
                              &arrays[SecondaryScattered], &arrays[SecondaryTransparent]});
        }

        // add the polarization components, if requested
        if (_recordPolarization)
        {
            colNames.insert(colNames.end(), {"total Stokes Q", "total Stokes U", "total Stokes V"});
            arrayPtrs.insert(arrayPtrs.end(), {&arrays[TotalQ], &arrays[TotalU], &arrays[TotalV]});

            if (_recordComponents && !_recordTotalOnly)
            {
                colNames.insert(
                    colNames.end(),
                    {"transparent Stokes Q", "transparent Stokes U", "transparent Stokes V", "direct primary Stokes Q",
                     "direct primary Stokes U", "direct primary Stokes V", "scattered primary Stokes Q",
                     "scattered primary Stokes U", "scattered primary Stokes V", "direct secondary Stokes Q",
                     "direct secondary Stokes U", "direct secondary Stokes V", "scattered secondary Stokes Q",
                     "scattered secondary Stokes U", "scattered secondary Stokes V", "transparent secondary Stokes Q",
                     "transparent secondary Stokes U", "transparent secondary Stokes V"});
                arrayPtrs.insert(arrayPtrs.end(),
                                 {&arrays[TransparentQ], &arrays[TransparentU], &arrays[TransparentV],
                                  &arrays[PrimaryDirectQ], &arrays[PrimaryDirectU], &arrays[PrimaryDirectV],
                                  &arrays[PrimaryScatteredQ], &arrays[PrimaryScatteredU], &arrays[PrimaryScatteredV],
                                  &arrays[SecondaryDirectQ], &arrays[SecondaryDirectU], &arrays[SecondaryDirectV],
                                  &arrays[SecondaryScatteredQ], &arrays[SecondaryScatteredU],
                                  &arrays[SecondaryScatteredV], &arrays[SecondaryTransparentQ],
                                  &arrays[SecondaryTransparentU], &arrays[SecondaryTransparentV]});
            }
        }

        // add the scattering levels, if requested, even if they are all zero
        if (_recordComponents)
        {
            for (int i = 0; i != _numScatteringLevels; ++i)
            {
                colNames.push_back(std::to_string(i + 1) + "-times scattered primary flux");
                arrayPtrs.push_back(_recordTotalOnly ? &empty : &arrays[PrimaryScatteredLevel + i]);
            }
        }
    };

    // local function to build synchronized lists of file names and recorded array pointers
    // params: vector of recorded arrays (in), list of file names (out), list of array pointers (out)
    auto buildFiles = [this, &total](vector<Array>& arrays, vector<string>& fileNames, vector<Array*>& arrayPtrs) {
        // add the total flux; if we didn't record it directly, calculate it now
        fileNames.push_back("total");
        if (_recordTotalOnly)
            arrayPtrs.push_back(&arrays[Total]);
        else
        {
            total = arrays[PrimaryDirect] + arrays[PrimaryScattered];
            if (_hasMediumEmission) total += arrays[SecondaryDirect] + arrays[SecondaryScattered];
            arrayPtrs.push_back(&total);
        }

        // add the flux components, if requested
        if (_recordComponents)
        {
            // add the transparent flux only if it may differ from the total flux
            if (!_recordTotalOnly)
            {
                fileNames.push_back("transparent");
                arrayPtrs.push_back(&arrays[Transparent]);
            }
            // add the actual components of the total flux (empty arrays will be ignored later on)
            fileNames.insert(fileNames.end(), {"primarydirect", "primaryscattered", "secondarytransparent",
                                               "secondarydirect", "secondaryscattered"});
            arrayPtrs.insert(arrayPtrs.end(),
                             {&arrays[PrimaryDirect], &arrays[PrimaryScattered], &arrays[SecondaryTransparent],
                              &arrays[SecondaryDirect], &arrays[SecondaryScattered]});
        }

        // add the polarization components, if requested
        if (_recordPolarization)
        {
            fileNames.insert(fileNames.end(), {"stokesQ", "stokesU", "stokesV"});
            arrayPtrs.insert(arrayPtrs.end(), {&arrays[TotalQ], &arrays[TotalU], &arrays[TotalV]});
        }

        // add the scattering levels, if requested
        if (!_recordTotalOnly)
            for (int i = 0; i != _numScatteringLevels; ++i)
            {
                fileNames.push_back("primaryscattered" + std::to_string(i + 1));
                arrayPtrs.push_back(&arrays[PrimaryScatteredLevel + i]);
            }
    };

    // local function to construct header comment line for text column output files
    auto header = [this, units](string name) {
        string header = "# " + name + " at ";
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
        return header;
    };

    // local function that returns unique pointer to new observer info record for distant instruments
    // or a "null" unique pointer for local instruments
    auto observerInfo = [this, units]() {
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
        return obsInfo;
    };

    // ---------------------- SED: flux density ----------------------

    // write SEDs to a single text file with multiple columns
    if (_includeFluxDensity)
    {
        // convert from recorded quantities to output quantities and from internal units to user-selected output units
        // (for performance reasons, determine the units scaling factor only once for each wavelength)
        for (int ell = 0; ell != _numWavelengths; ++ell)
        {
            double factor = 1. / fourpid2 / _lambdagrid->effectiveWidth(ell)
                            * units->ofluxdensity(_lambdagrid->wavelength(ell), 1.);
            for (auto& array : _sed)
                if (array.size()) array[ell] *= factor;
        }

        // build a list of column names and corresponding pointers to sed arrays (which may be empty)
        vector<string> sedNames;
        vector<Array*> sedArrays;
        buildCols(_sed, sedNames, sedArrays);

        // open the file and add the column headers
        TextOutFile sedFile(_parentItem, _instrumentName + "_sed", "SED");
        sedFile.writeLine(header("SED"));
        sedFile.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
        for (const string& name : sedNames)
        {
            sedFile.addColumn(name + "; " + units->sfluxdensity(), units->ufluxdensity());
        }

        // write the column data
        for (int ell : Indices(_numWavelengths, units->rwavelength()))
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
            for (int ell : Indices(_numWavelengths, units->rwavelength()))
            {
                vector<double> values({units->owavelength(_lambdagrid->wavelength(ell))});
                for (int k = 0; k <= maxContributionPower; ++k) values.push_back(_wsed[k][ell]);
                statFile.writeRow(values);
            }
            statFile.close();
        }
    }

    // ---------------------- IFU: surface brightness ----------------------

    // write IFUs to FITS files (one file per IFU)
    if (_includeSurfaceBrightness)
    {
        // convert from recorded quantities to output quantities and from internal units to user-selected output units
        // (for performance reasons, determine the units scaling factor only once for each wavelength)
        for (int ell = 0; ell != _numWavelengths; ++ell)
        {
            double factor = 1. / fourpid2 / omega / _lambdagrid->effectiveWidth(ell)
                            * units->osurfacebrightness(_lambdagrid->wavelength(ell), 1.);
            size_t begin = ell * _numPixelsInFrame;
            size_t end = begin + _numPixelsInFrame;
            for (auto& array : _ifu)
                if (array.size())
                    for (size_t lell = begin; lell != end; ++lell) array[lell] *= factor;
        }

        // copy the wavelength grid in output units
        Array wavegrid(_numWavelengths);
        for (int ell = 0; ell != _numWavelengths; ++ell)
            wavegrid[ell] = units->owavelength(_lambdagrid->wavelength(ell));

        // reverse the ordering of the wavelength grid and frames if necessary
        if (units->rwavelength())
        {
            NR::reverse(wavegrid);

            // flux frames
            for (auto& array : _ifu)
                if (array.size()) NR::reverse(array, _numPixelsInFrame);

            // statistics frames
            for (auto& array : _wifu)
                if (array.size()) NR::reverse(array, _numPixelsInFrame);
        }

        // build a list of file names and corresponding pointers to ifu arrays (which may be empty)
        vector<string> ifuNames;
        vector<Array*> ifuArrays;
        buildFiles(_ifu, ifuNames, ifuArrays);

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

        // output the files (ignoring empty arrays)
        auto info = observerInfo();
        int numFiles = ifuNames.size();
        for (int q = 0; q != numFiles; ++q)
        {
            if (ifuArrays[q]->size())
            {
                string filename = _instrumentName + "_" + ifuNames[q];
                string description = ifuNames[q] + " flux";
                FITSInOut::write(_parentItem, description, filename, *(ifuArrays[q]), units->usurfacebrightness(),
                                 _numPixelsX, _numPixelsY, incx, incy, cx, cy, unitsxy, wavegrid, units->uwavelength(),
                                 info.get());
            }
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

    // ---------------------- LC: light curve ----------------------

    // write LCs to a single text file with multiple columns
    if (_includeLightCurve)
    {
        int numTimes = _timegrid->numBins();
        int numArrays = _lc.size();

        // calibrate and convert values to output units
        for (int k = 0; k != numTimes; ++k)
        {
            double factor = 1. / fourpid2 / _timegrid->width(k);
            for (int i = 0; i != numArrays; ++i)
                if (_lc[i].size()) _lc[i][k] = factor * units->otimefluxdensity(_lc[i][k], _lcw[i][k]);
        }

        // build a list of column names and corresponding pointers to lc arrays (which may be empty)
        vector<string> lcNames;
        vector<Array*> lcArrays;
        buildCols(_lc, lcNames, lcArrays);

        // open the file and add the column headers
        TextOutFile lcFile(_parentItem, _instrumentName + "_lc", "LC");
        lcFile.writeLine(header("LC"));
        lcFile.addColumn("time lag", units->utimelag());
        for (const string& name : lcNames)
        {
            lcFile.addColumn(name + "; " + units->stimefluxdensity(), units->utimefluxdensity());
        }

        // write the column data
        for (int k = 0; k != numTimes; ++k)
        {
            vector<double> values({units->otimelag(_timegrid->time(k))});
            for (const Array* array : lcArrays) values.push_back(array->size() ? (*array)[k] : 0.);
            lcFile.writeRow(values);
        }
        lcFile.close();
    }

    // ---------------------- STM: spectral-time map ----------------------

    // write STMs to FITS files (one file per STM)
    if (_includeSpectralTimeMap)
    {
        int numTimes = _timegrid->numBins();

        // calibrate and convert values to output units
        for (int k = 0; k != numTimes; ++k)
        {
            for (int ell = 0; ell != _numWavelengths; ++ell)
            {
                int kell = ell + k * _numWavelengths;
                double factor = 1. / fourpid2 / _lambdagrid->effectiveWidth(ell) / _timegrid->width(k)
                                * units->ospectraltimefluxdensity(_lambdagrid->wavelength(ell), 1.);
                for (auto& array : _stm)
                    if (array.size()) array[kell] *= factor;
            }
        }

        // copy the wavelength grid in output units
        Array wavegrid(_numWavelengths);
        for (int ell = 0; ell != _numWavelengths; ++ell)
            wavegrid[ell] = units->owavelength(_lambdagrid->wavelength(ell));

        // copy the time grid in output units
        Array timegrid(numTimes);
        for (int k = 0; k != numTimes; ++k) timegrid[k] = units->otimelag(_timegrid->time(k));

        // reverse the ordering of the wavelength grid and flux data if necessary
        if (units->rwavelength())
        {
            // reverse the wavelength grid
            NR::reverse(wavegrid);

            // reverse the flux data
            for (auto& array : _stm)
            {
                if (array.size())
                {
                    for (int k = 0; k < numTimes; ++k)
                    {
                        double* begin = &array[k * _numWavelengths];
                        double* end = begin + _numWavelengths;
                        std::reverse(begin, end);
                    }
                }
            }
        }

        // build a list of file names and corresponding pointers to ifu arrays (which may be empty)
        vector<string> stmNames;
        vector<Array*> stmArrays;
        buildFiles(_stm, stmNames, stmArrays);

        // output the files (ignoring empty arrays)
        auto info = observerInfo();
        int numFiles = stmNames.size();
        for (int q = 0; q != numFiles; ++q)
        {
            if (stmArrays[q]->size())
            {
                string filename = _instrumentName + "_stm_" + stmNames[q];
                string description = stmNames[q] + " flux";
                FITSInOut::writeMap(_parentItem, description, filename, *(stmArrays[q]),
                                    units->uspectraltimefluxdensity(), wavegrid, timegrid, units->uwavelength(),
                                    units->utimelag(), info.get());
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
