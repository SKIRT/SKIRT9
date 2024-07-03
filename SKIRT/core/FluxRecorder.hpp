/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FLUXRECORDER_HPP
#define FLUXRECORDER_HPP

#include "Array.hpp"
#include "ThreadLocalMember.hpp"
#include <tuple>
class MediumSystem;
class PhotonPacket;
class SimulationItem;
class WavelengthGrid;

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
    transparent            | Transparent flux from primary sources | \em recordComponents = true
    primarydirect          | Direct, unscattered flux from primary sources | \em recordComponents = true
    primaryscattered       | Indirect, scattered flux from primary sources | \em recordComponents = true
    primaryscatteredN      | Flux from primary sources that has scattered N times | .. & \em numScatteringLevels > 0
    secondarytransparent   | Transparent flux emitted by media | \em recordComponents = true
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
    Transparent fluxes and flux components | \em recordComponents = true
    Stokes vector elements                 | \em recordPolarization = true
    N-times scattered primary flux         | \em recordComponents = true & \em numScatteringLevels > 0

    The second file, called <tt>prefix_instr_sedstats.txt</tt>, is written only if statistics are
    requested. It includes a column for the wavelength plus a column for each of the individual
    photon contribution sums, for powers from zero to 4.

    Calling sequence
    ----------------

    A FluxRecorder instance expects a rigourous calling sequence. During setup, the instrument
    configures the FluxRecorder's operation, specifying the wavelength grid, the items to be
    recorded (%SED, IFU, individual flux components, statistics), flux calibration settings
    (instrument type, distance, redshift), and some additional information on the simulation in
    which the instrument is embedded (e.g., is there any secondary emission). The configuration
    must be completed by calling the finalizeConfiguration() function.

    To record the effects of detecting a photon packet, the instrument invokes the detect()
    function. This function is thread-safe, so it may be (and often is) called from multiple
    parallel execution threads. After the parallel threads have completed the work on a series of
    photon packets, and before the parallel threads are actually destructed, the instrument should
    call the flush() function from a single thread to process any information buffered by the
    detect() function in thread-local storage. Finally, at the end of the simulation, the
    instrument calls the calibrateAndWrite() function to output the recorded information.

    Flux calibration
    ----------------

    A FluxRecorder instance handles the conversion of the detected photon packet contributions to a
    corresponding observed quantity for each spectral and/or spatial bin. This process is called
    flux calibration and it depends to some extent on whether the instrument associated with the
    flux recorder is "distant" or "local".

    A distant instrument is placed at a large distance from the model under study. The observed
    portion of the sky is approximated by a plane perpendicular to the line of sight rather than a
    part of a sphere. The instrument uses parallel projection and considers all photon packets to
    originate at the same distance from the instrument. Consequently, the incoming photon packet
    contributions can simply be accumulated for each bin and flux calibration can occur on the
    aggregated result. Assuming a non-relativistic distance, the calibration for the flux densities
    \f$F_\lambda\f$ in an %SED is \f[ F_\lambda = \frac{1}{\Delta\lambda} \, \frac{1}{4\pi
    d_\mathrm{m}^2} \, L_\mathrm{bin}\f] where \f$L_\mathrm{bin}\f$ is the aggregated photon packet
    luminosity contribution in a spectral bin, \f$\Delta\lambda\f$ is the width of that bin, and
    \f$d_\mathrm{m}\f$ is the distance from the model to the instrument. The calibration for the
    surface brightness values \f$f_\lambda\f$ in an IFU is \f[ f_\lambda = \frac{1}{\Delta\lambda}
    \, \frac{1}{4\pi d_\mathrm{m}^2} \, \frac{1}{\Omega_\mathrm{distant}}\, L_\mathrm{bin}\f] where
    \f$\Omega_\mathrm{distant}\f$ is the solid angle subtended by the detector pixel corresponding
    to the bin. For a distant instrument using parallel projection, the latter is given by
    \f[\Omega_\mathrm{distant} = 4\arctan\left(\frac{s_\mathrm{x}}{2d_\mathrm{m}}\right)
    \arctan\left(\frac{s_\mathrm{y}}{2d_\mathrm{m}}\right) \f] where \f$s_\mathrm{x}\f$ and
    \f$s_\mathrm{y}\f$ represent the field of view (in model units) of a pixel in each direction.
    For relativistic distances these relations become slightly more involved, as discussed in a
    seperate concept note on cosmological redshift.

    A local instrument is placed near or even inside the model. The mapping between incoming photon
    packet trajectories and detector pixels is determined by the instrument. Although all distances
    are considered to be non-relativistic, photon packets can now originate at various distances
    from the instrument, with a wide dynamic range. The FluxRecorder class does not support
    recording %SEDs for local instruments. The calibration for the surface brightness values
    \f$f_\lambda\f$ in an IFU now becomes \f[ f_\lambda = \frac{1}{\Delta\lambda} \, \frac{1}{4\pi
    d_\mathrm{pp}^2} \, \frac{1}{\Omega_\mathrm{local}}\, L_\mathrm{pp}\f] where
    \f$L_\mathrm{pp}\f$ is the luminosity contribution by a particular photon packet,
    \f$d_\mathrm{pp}\f$ is the distance of the packet's originating position to the instrument, and
    \f$\Omega_\mathrm{local}\f$ is the solid angle subtended by the detector pixel corresponding to
    the bin. The correction for distance now must be performed for each photon packet individually.
    Furthermore, the solid angle subtended by a detector pixel now depends on the projection
    applied by the instrument. The current implementation assumes that the solid angle is identical
    for all pixels in the detector, alhough this is not generally true for all instrument
    projections. In that case, the instrument should specify a representative value and advise the
    user to avoid situations where the actual solid angles deviate much from this value.

    Memory usage
    ------------

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
    /** The constructor initializes the FluxRecorder to a configuration that records nothing. The
        argument specifies a simulation item in the hierarchy of the caller (usually the caller
        itself). It is used to get a human-readable name for the caller and to retrieve a logger
        and a unit system from the simulation hierarchy. */
    FluxRecorder(const SimulationItem* parentItem);

    /** This function configures information on the simulation in which the recorder is embedded.
        In order of appearance, the arguments specify the name of the associated instrument, the
        wavelength grid of the instrument, whether the simulation includes at least some media, and
        whether the simulation includes emission from those media. */
    void setSimulationInfo(string instrumentName, const WavelengthGrid* lambdagrid, bool hasMedium,
                           bool hasMediumEmission);

    /** This function configures the user requirements for recording, respectively, flux
        components, individual scattering level contributions, polarization, and information for
        calculating statistics. See the documentation in the header of this class for more
        information. */
    void setUserFlags(bool recordComponents, int numScatteringLevels, bool recordPolarization, bool recordStatistics);

    /** This function sets the observer angles for a distant instrument associated with this flux
        recorder. These values are listed in the output files as a convenience to the user but are
        not otherwise used. This function should not be called for a local instrument. */
    void setObserverAngles(double inclination, double azimuth, double roll);

    /** This function configures the distance of the recorder in the model's rest frame for a
        distant instrument. The specified distance must be nonzero. The client must call either the
        setRestFrameDistance() or setObserverFrameRedshift() functions, not both. This function
        should not be called for a local instrument. */
    void setRestFrameDistance(double distance);

    /** This function configures the redshift and relativistic distances of the recorder's observer
        frame for a distant instrument. The specified redshift and distances must be nonzero. The
        client must call either the setRestFrameDistance() or setObserverFrameRedshift() functions,
        not both. This function should not be called for a local instrument. */
    void setObserverFrameRedshift(double redshift, double angularDiameterDistance, double luminosityDistance);

    /** This function enables recording of spatially integrated flux densities, i.e. an %SED, for a
        distant instrument. Recording of spatially integrated flux densities for local instruments
        is not supported, so this function should not be called for a local instrument. */
    void includeFluxDensityForDistant();

    /** This function enables recording of IFU data cubes, i.e. a surface brightness image frame
        for each wavelength, for a distant instrument. The number of pixels and the pixel sizes are
        used to calibrate the surface brightness; the center coordinates are used only for the
        metadata in the output file. The coordinates and pixel sizes are converted to angular sizes
        before being written to the metadata of the output file. This function should not be called
        for a local instrument. */
    void includeSurfaceBrightnessForDistant(int numPixelsX, int numPixelsY, double pixelSizeX, double pixelSizeY,
                                            double centerX, double centerY);

    /** This function enables recording of IFU data cubes, i.e. a surface brightness image frame
        for each wavelength, for a local instrument. The solid ange per pixel is used to calibrate
        the surface brightness. The remaining arguments are used only for the metadata in the
        output file. If the quantity string is nonempty, it should specify a valid SKIRT quantity
        name and the increment and center values are converted to the corresponding output units.
        This function should not be called for a distant instrument. */
    void includeSurfaceBrightnessForLocal(int numPixelsX, int numPixelsY, double solidAnglePerPixel, double incrementX,
                                          double incrementY, double centerX, double centerY,
                                          string quantityXY = string());

    /** This function completes the configuration of the recorder. It must be called after any of
        the configuration functions, and before the first invocation of the detect() function. */
    void finalizeConfiguration();

    //======================== Other Functions =======================

public:
    /** This function simulates the detection of a photon packet by the recorder as determined by
        its configuration. This function is thread-safe, so it may be (and often is) called from
        multiple parallel execution threads.

        The calling instrument is responsible for providing the index \em l of the pixel in the
        instrument frame where the photon packet arrives, because this depends on the projection
        being used. In addition, the instrument can specify a \em distance from the photon packet's
        last interaction site to the instrument. For distant instruments with parallel projection,
        this distance should be left at its default value of infinity. For instruments that may be
        placed close by or inside the model, the actual distance should be specified so that the
        flux can be properly calibrated for each individual photon packet.

        All other information is obtained directly or indirectly from the photon packet. If there
        is an obscuring medium, the optical depth from the photon packet's last interaction site to
        the instrument is determined and the corresponding extincton is applied to the packet's
        contribution before detection. */
    void detect(PhotonPacket* pp, int l, double distance = std::numeric_limits<double>::infinity());

    /** This function processes and clears any information that may have been buffered by the
        detect() function in thread-local storage. It is not thread-safe. After parallel threads
        have completed the work on a series of photon packets, and before the parallel threads are
        actually destructed, the flush() function should be called from a single thread. */
    void flush();

    /** This function calibrates and outputs the instrument data. The calibration includes dividing
        the luminosities (W) recorded for each bin by the wavelength bin width to obtain specific
        luminosities (W/m) and further conversion to flux density (incorporating distance) and/or
        to surface brightness (incorporating distance and solid angle per pixel) based on the
        information passed during configuration. The function also converts the resulting values
        from internal units to output units depending on the simulation's choices for flux output
        style.

        For more information on flux calibration and on the names and contents of the generated
        files, see the documentation in the header of this class. */
    void calibrateAndWrite();

    //================= Private Types and Functions ===============

private:
    /** Private data structure to remember a single contribution from a photon packet to a
        statistics bin. */
    class Contribution
    {
    public:
        Contribution(int ell, int l, double w) : _ell(ell), _l(l), _w(w) {}
        bool operator<(const Contribution& c) const { return std::tie(_ell, _l) < std::tie(c._ell, c._l); }
        int ell() const { return _ell; }
        int l() const { return _l; }
        double w() const { return _w; }

    private:
        int _ell{0};   // wavelength index
        int _l{0};     // pixel index (relevant only for IFUs)
        double _w{0};  // contribution
    };

    /** Private data structure to remember a list of contributions for a given photon packet
        history. We assume that all detections for a given history are handled inside the same
        execution thread and that histories (within a particular thread) are handled one after the
        other (i.e. not interleaved). */
    class ContributionList
    {
    public:
        bool hasHistoryIndex(size_t historyIndex) const { return _historyIndex == historyIndex; }
        void addContribution(int ell, int l, double w) { _contributions.emplace_back(ell, l, w); }
        void reset(size_t historyIndex = 0) { _historyIndex = historyIndex, _contributions.clear(); }
        void sort() { std::sort(_contributions.begin(), _contributions.end()); }
        const vector<Contribution>& contributions() const { return _contributions; }

    private:
        size_t _historyIndex{0};
        vector<Contribution> _contributions;
    };

    /** This private helper function records the photon packet history contributions in the
        specified list into the statistics arrays. */
    void recordContributions(ContributionList* contributionList);

    //======================== Data Members ========================

private:
    // information on the simulation being recorded, received from client during configuration
    const SimulationItem* _parentItem{nullptr};
    string _instrumentName;
    const WavelengthGrid* _lambdagrid{nullptr};
    bool _hasMedium{false};
    bool _hasMediumEmission{false};  // relevant only when hasMedium is true

    // recorder configuration, received from client during configuration
    bool _recordComponents{false};
    int _numScatteringLevels{false};  // honored only when recordComponents is true
    bool _recordPolarization{false};
    bool _recordStatistics{false};
    bool _includeFluxDensity{false};
    bool _includeSurfaceBrightness{false};

    // recorder configuration on observer angles, received from client during configuration
    double _inclination{0};
    double _azimuth{0};
    double _roll{0};

    // recorder configuration on distance and redshift, received from client during configuration
    double _redshift{0};
    double _angularDiameterDistance{0};
    double _luminosityDistance{0};

    // recorder configuration for IFUs, received from client during configuration
    bool _local{false};
    int _numPixelsX{0};
    int _numPixelsY{0};
    double _pixelSizeX{0};
    double _pixelSizeY{0};
    double _solidAnglePerPixel{0};
    double _centerX{0};
    double _centerY{0};
    double _incrementX{0};
    double _incrementY{0};
    string _quantityXY;

    // cached info, initialized when configuration is finalized
    MediumSystem* _ms{nullptr};   // pointer to medium system, if present (used only if hasMedium is true)
    bool _recordTotalOnly{true};  // becomes false if recordComponents and hasMedium are both true
    size_t _numPixelsInFrame{0};  // number of pixels in a single IFU frame

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
