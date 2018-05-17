/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FLUXRECORDER_HPP
#define FLUXRECORDER_HPP

#include "Array.hpp"
#include "ThreadLocalMember.hpp"
#include <tuple>
class PhotonPacket;
class WavelengthGrid;

////////////////////////////////////////////////////////////////////

// implementation detail
namespace FluxRecorder_Impl { class ContributionList; }

////////////////////////////////////////////////////////////////////

/** FluxRecorder is a helper class used by instruments to actually record the effects of detected
    photon packets. Each Instrument instance employs a single FluxRecorder instance.

    Features
    --------

    Depending on an instrument's needs, a FluxRecorder instance records an %SED (spectral energy
    density) with a spatially integrated flux density for each wavelength bin, and/or a full IFU
    (integral field unit data cube) with a surface brightness for each pixel in an image frame for
    each wavelength bin. In each case, flux information can be recorded and written to file
    depending on user configuration.

    Specifically, the class can record the individual contributions to the total flux from primary
    and secondary sources (i.e. emission by media), in each case splitting direct and scattered
    radiation. It can also store the flux as it would be seen without any attenuation, i.e. a
    transparent view of the system as if there were no media. Furthermore, the contributions from
    individual scattering levels (for primary radiation only) up to a certain maximum level can be
    recorded separately. When the maximum level is set to zero, this feature is disabled. Finally,
    the class can also record the elements of the Stokes vector for the total flux, \f$(I, Q, U,
    V)\f$, where the intensity \f$I\f$ corresponds to the total flux itself. Note that the Stokes
    \f$Q, U, V\f$ "flux" values can be negative.

    Upon request, the class can also record information intended for calculating statistical
    properties of the results. Let \f$N\f$ denote the number of primary and secondary photon
    packets launched during the peel-off segments of the simulation. We define \f$w_i,
    i=1,\dots,N\f$ as the contribution of the \f$i\f$th photon packet to a particular bin (a
    wavelength bin for an %SED, or a pixel in one of the wavelength frames for an IFU). The value
    of \f$w_i\f$ includes the contributions of all photon packets peeled-off and/or scattered from
    the originally launched packet (called the \em history of that packet). Given this definition,
    the class tracks and outputs the sums \f$\sum_i w_i^k\f$, with \f$k=0,\dots,4\f$ for each bin.
    These sums allow calculating second order statistical properties such as the relative error
    \f$R\f$ and fourth order statistical properties such as the variance of the variance VOV. For
    more information, see, e.g., the user manual for the MCNP code (A General Monte Carlo
    N-Particle Transport Code, Version 5, April 24, 2003, Revised 2/1/2008, Los Alamos National
    Laboratory, USA) or Camps and Baes 2018 (ApJ).

    All of these output possibilities are summarized in the table below. A separate IFU output file
    is written for each line in the table; the first column in the table lists the portion of the
    output filename <tt>prefix_instr_XXX.fits</tt> indicating the IFU type. Note that an output
    file is created only if the corresponding information has been requested \em and it is
    meaningful. For example, if the transparent flux is know to be identical to the total flux
    (because there are no media), the transparent file is not written. Also, if the simulation does
    not include media emission, the secondary flux files are not written.

    IFU file name          | Description | Configured by
    -----------------------|-------------|--------------
    total                  | Total attenuated flux | always on
    transparent            | Transparent flux, as if there were no attenuating media | \em recordComponents = true
    primarydirect          | Direct, unscattered flux from primary sources | \em recordComponents = true
    primaryscattered       | Indirect, scattered flux from primary sources | \em recordComponents = true
    primaryscatteredN      | Flux from primary sources that has scattered N times | .. & \em numScatteringLevels > 0
    secondarydirect        | Direct, unscattered flux emitted by media | \em recordComponents = true
    secondaryscattered     | Indirect, scattered flux emitted by media | \em recordComponents = true
    stokesQ                | Stokes vector element Q for total flux | \em recordPolarization = true
    stokesU                | Stokes vector element U for total flux | \em recordPolarization = true
    stokesV                | Stokes vector element V for total flux | \em recordPolarization = true
    statsN                 | Sum of individual photon contributions to the power of N | \em recordStatistics = true

    The %SED information is written in a maximum of two text column files. The first file, called
    <tt>prefix_instr_sed.txt</tt>, includes a column for each of the requested flux components in
    the order indicated by the table below. The columns for each block are either all present or
    all absent as requested in the configuration, even if the numbers are zero or identical to
    other columns. This allows the column indices to be determined solely based on knowledge of the
    instrument configuration flags, without taking into account whether the simulation actually
    includes the corresponding feature (e.g., media emission or polarization).

    Block                                  | Configured by
    ---------------------------------------|--------------
    Total flux                             | always on
    Transparent flux and 4 flux components | \em recordComponents = true
    Stokes vector elements                 | \em recordPolarization = true
    N-times scattered primary flux         | \em recordComponents = true & \em numScatteringLevels > 0

    The second file, called <tt>prefix_instr_sedstats.txt</tt>, is written only if statistics are
    requested. It includes a column for the wavelength plus a column for each of the individual
    photon contribution sums, for powers from zero to 4.

    Usage
    -----

    A FluxRecorder instance expects a rigourous calling sequence. During setup, the instrument
    configures the FluxRecorder's operation, specifying the wavelength grid, the items to be
    recorded (SED, IFU, individual flux components, statistics) and some additional information on
    the simulation in which the instrument is embedded (e.g., is there any secondary emission). The
    configuration must be completed by calling the finalizeConfiguration() function.

    To record the effects of detecting a photon packet, the instrument invokes the detect()
    function. This function is thread-safe, so it may be (and often is) called from multiple
    parallel execution threads. After the parallel threads have completed the work on a series of
    photon packets, and before the parallel threads are actually destructed, the instrument should
    call the flush() function from a single thread to process any information buffered by the
    detect() function in thread-local storage. Finally, at the end of the simulation, the
    instrument calls the calibrateAndWrite() function to output the recorded information.

    A FluxRecorder instance dynamically adjusts its memory allocation to the configuration and
    simulation characteristics. Detector arrays for individual flux components, polarization, or
    statistics are allocated only when requested in the configuration. Also, for example, if there
    is no secondary emission in the simulation, the corresponding detector arrays are not
    allocated, even if recording of individual components is requested in the configuration.
*/
class FluxRecorder final
{
    //============= Construction - Setup - Destruction =============

public:
    /** The constructor initializes the FluxRecorder to a configuration that records nothing. */
    FluxRecorder();

    /** This function configures information on the simulation in which the recorder is embedded.
        In order of appearance, the arguments specify the name of the associated instrument, the
        wavelength grid of the instrument, whether the simulation includes at least some media, and
        whether the simulation includes emission from those media. */
    void setSimulationInfo(string instrumentName, const WavelengthGrid* lambdagrid,
                           bool hasMedium, bool hasMediumEmission);

    /** This function configures the user requirements for recording, respectively, flux
        components, individual scattering level contributions, polarization, and information for
        calculating statistics. See the documentation in the header of this class for more
        information. */
    void setUserFlags(bool recordComponents, int numScatteringLevels,
                      bool recordPolarization, bool recordStatistics);

    /** This function enables recording of spatially integrated flux densities, i.e. an %SED,
        assuming parallel projection at the specified instrument distance from the model. If both
        includeFluxDensity() and includeSurfaceBrightness() are called, the specified distances
        must be the same. */
    void includeFluxDensity(double distance);

    /** This function enables recording of IFU data cubes, i.e. a surface brightness image frame
        for each wavelength. The recorder assumes parallel projection at the specified instrument
        distance from the model and using the specified frame properties. The distance, the number
        of pixels, and the pixel sizes are used to calibrate the surface brightness; the center
        coordinates are used only for the metadata in the output file. If both includeFluxDensity()
        and includeSurfaceBrightness() are called, the specified distances must be the same. */
    void includeSurfaceBrightness(double distance, int numPixelsX, int numPixelsY,
                                  double pixelSizeX, double pixelSizeY, double centerX, double centerY);

    /** This function completes the configuration of the recorder. It must be called after any of
        the configuration functions, and before the first invocation of the detect() function. */
    void finalizeConfiguration();

    //======================== Other Functions =======================

public:
    /** This function simulates the detection of a photon packet by the recorder as determined by
        its configuration. This function is thread-safe, so it may be (and often is) called from
        multiple parallel execution threads.

        If the luminosity of the photon packet at the site of its last emission or scattering event
        is equal to \f$L_\nu\f$, the fraction that will reach the observer is equal to \f[
        L_\nu^{\text{obs}} = L_\nu\, {\text{e}}^{-\tau_{\ell,{\text{path}}}} \f] with
        \f$\tau_{\ell,{\text{path}}}\f$ the optical depth through the media system towards the
        observer. */
    void detect(const PhotonPacket* pp, int l, double tau);

    /** This function processes and clears any information that may have been buffered by the
        detect() function in thread-local storage. It is not thread-safe. After parallel threads
        have completed the work on a series of photon packets, and before the parallel threads are
        actually destructed, the flush() function should be called from a single thread. */
    void flush();

    /** This function calibrates and outputs the instrument data. The calibration includes
        conversion from specific luminosity to flux density (incorporating distance) and/or to
        surface brightness (i.e. incorporating distance and solid angle per pixel) based on the
        information passed during configuration. The function also converts the resulting values
        from internal units to output units depending on the simulation's choices for flux output
        style.

        For more information on the names and contents of the generated files, see the
        documentation in the header of this class. */
    void calibrateAndWrite();

    //================= Private Types and Functions ===============

private:
    using ContributionList = FluxRecorder_Impl::ContributionList;

    /** This private helper function records the photon history contributions in the specified list
        into the statistics arrays. */
    void recordContributions(ContributionList& contributionList);

    //======================== Data Members ========================

private:
    // information on the simulation being recorded, received from client during configuration
    string _instrumentName;
    const WavelengthGrid* _lambdagrid{nullptr};
    bool _hasMedium{false};
    bool _hasMediumEmission{false};     // relevant only when hasMedium is true

    // recorder configuration, received from client during configuration
    bool _recordComponents{false};
    int _numScatteringLevels{false};    // honored only when recordComponents is true
    bool _recordPolarization{false};
    bool _recordStatistics{false};
    bool _includeFluxDensity{false};
    bool _includeSurfaceBrightness{false};

    // recorder configuration for SEDs and/or IFUs, received from client during configuration
    double _distance{0};

    // recorder configuration for IFUs, received from client during configuration
    int _numPixelsX{0};
    int _numPixelsY{0};
    double _pixelSizeX{0};
    double _pixelSizeY{0};
    double _centerX{0};
    double _centerY{0};

    // cached info, initialized when configuration is finalized
    bool _recordTotalOnly{true};        // becomes false if recordComponents and hasMedium are both true
    size_t _numPixelsInFrame{0};        // number of pixels in a single IFU frame

    // detector arrays that need to be calibrated, initialized when configuration is finalized
    vector<Array> _sed;
    vector<Array> _ifu;

    // detector arrays for statistics that should not be calibrated, initialized when configuration is finalized
    vector<Array> _wsed;
    vector<Array> _wifu;

    // thread-local contribution list
    ThreadLocalMember<ContributionList> _contributionLists;
};

////////////////////////////////////////////////////////////////////

#endif
