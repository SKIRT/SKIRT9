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

/** FluxRecorder is a helper class used by instruments to actually recording the effects of detected photon packets.

    XXX  object represents an advanced instrument that records individual contributions
    to the flux from various sources, including the components of the Stokes vector for polarized
    radiation. For each type of flux contribution being recorded, a FullInstrument object holds two
    array data members: a simple 1D array (the F-array) that stores the integrated flux at every
    wavelength index, and a 3D array (the f-array) corresponding to the surface brightness in every
    pixel, at every wavelength index. The instrument dynamically adapts to the characteristics of
    the simulation it is part of. For example, if there is no dust system, only the direct stellar
    flux is recorded (since there is nothing else). If dust emission is turned off in the dust
    system, the thermal radiation emitted by dust is not recorded (since it will be zero).
    Similarly, if none of the dust mixtures in the dust sytem support polarization, the components
    of the Stokes vector are not recorded. These adjustments save memory and processing time where
    appropriate.

    Specifically, this instrument separately records the contributions to the total flux
    \f$F_\lambda^{\text{tot}}\f$ due to the direct stellar radiation
    \f$F_\lambda^{*,\text{dir}}\f$, the scattered stellar radiation \f$F_\lambda^{*,\text{sca}}\f$,
    the direct dust radiation \f$F_\lambda^{\text{dust},\text{dir}}\f$, and the scattered dust
    radiation \f$F_\lambda^{\text{dust},\text{sca}}\f$, with obviously \f[ F_\lambda^{\text{tot}} =
    F_\lambda^{*,\text{dir}} + F_\lambda^{*,\text{sca}} + F_\lambda^{\text{dust},\text{dir}} +
    F_\lambda^{\text{dust},\text{sca}}. \f] The instrument also stores the flux from the object as
    it would be seen without any dust attenuation, i.e. the transparent view of the system.
    Furthermore, the contributions from individual scattering levels (stellar radiation only)
    up to a certain maximum level \f$N_{\text{max}}\f$ can be recorded separately. By
    default, the maximum level is set to zero, effectively disabling this feature. Finally, the
    instrument also records the values corresponding to each of the elements of the Stokes vector,
    i.e. \f$Q_\lambda^{\text{tot}}\f$, \f$U_\lambda^{\text{tot}}\f$, and
    \f$V_\lambda^{\text{tot}}\f$. Note that these Stokes "flux" values can be negative. */
class FluxRecorder final
{
    //============= Construction - Setup - Destruction =============

public:
    /** The constructor initializes the configuration to record nothing. */
    FluxRecorder();

    /** This function ... */
    void setSimulationInfo(string instrumentName, const WavelengthGrid* lambdagrid,
                           bool hasMedium, bool hasMediumEmission);

    /** This function ... */
    void setUserFlags(bool recordComponents, int numScatteringLevels,
                      bool recordPolarization, bool recordStatistics);

    /** This function ... */
    void includeFluxDensity(double distance);

    /** This function ... */
    void includeSurfaceBrightness(double distance, int numPixelsX, int numPixelsY,
                                  double pixelSizeX, double pixelSizeY, double centerX, double centerY);

    /** This function ... */
    void finalizeConfiguration();

    //======================== Other Functions =======================

public:
    /** This function simulates the detection of a photon package by the instrument. The
        ingredients to be determined are the pixel that the photon package will hit and the
        luminosity that will be collected by the different subdetectors of the instrument (in the
        detection phase, the detectors store luminosities, these will be converted to fluxes and
        surface brightnesses once the simulation is finished). If the total luminosity of the
        photon at the site of its last emission or scattering event is equal to \f$L_\ell\f$, the
        fraction of it that will reach the observer is equal to \f[ L_\ell^{\text{tot}} = L_\ell\,
        {\text{e}}^{-\tau_{\ell,{\text{path}}}} \f] with \f$\tau_{\ell,{\text{path}}}\f$ the
        optical depth through the dust system towards the observer. Depending on the origin of the
        photon package (i.e. emitted by the stars or by the dust) and the number of scattering
        events the package has already experienced, this bit of luminosity must be assigned to the
        correct subdetector. The corresponding flux for the transparent case is simply
        \f$L_\ell^{\text{tra}}\f$. We now only have to add the luminosities to the stored
        luminosity in the correct bin of both the 1D F-vectors and the 3D f-vectors. */
    void detect(const PhotonPacket* pp, int l, double tau);

    /** This function ... */
    void flush();

    /** This function calibrates and outputs the instrument data. The calibration takes care of the
        conversion from bolometric luminosity units to flux density units (for the F-vector) and
        surface brightness units (for the f-vector). Depending on the characteristics of the
        simulation in which the instrument resides (see description in the class header), this
        function creates a number of FITS files named <tt>prefix_instr_total.fits</tt>,
        <tt>prefix_instr_direct.fits</tt>, <tt>prefix_instr_scattered.fits</tt>,
        <tt>prefix_instr_dust.fits</tt>, <tt>prefix_instr_transparent.fits</tt>,
        <tt>prefix_instr_stokesQ.fits</tt>, <tt>prefix_instr_stokesU.fits</tt>, and
        <tt>prefix_instr_stokesV.fits</tt> that contain the frames with the surface brightness in
        every pixel. Note that the Stokes vector components may have negative values. If there is
        only one wavelength, the FITS files are 2D; otherwise the FITS files are 3D. The function
        also creates \f$N_{\text{max}}\f$ separate FITS files named
        <tt>prefix_instr_scatteringlevel1.fits</tt>, etc. that contain the frames with the
        contribution of the individual scattering levels to the surface brightness. Finally, an
        ASCII file <tt>prefix_instrument_sed.dat</tt> is created that contains at least six
        columns, listing the wavelength \f$\lambda\f$ and the total flux, the direct flux, the
        scattered flux, the thermal dust flux, and the transparent flux for that wavelength. These
        six columns are always present, even if the corresponding value is trivially zero. If the
        simulation supports polarization, the subsequent columns list the three Q, U, V components
        of the Stokes vector (which may be negative). The last \f$N_{\text{max}}\f$ columns, if
        any, contain the contribution of the different scattering levels to the total flux. The
        units in which the surface brightness and the flux densities are written depends on the
        simulation's choices for flux output style (Neutral indicates \f$\lambda F_\lambda = \nu
        F_\nu\f$; Wavelength indicates \f$F_\lambda\f$; Frequency indicates \f$F_\nu\f$) and
        actual units. */
    void calibrateAndWrite();

    //================= Private Types and Functions ===============

private:
    using ContributionList = FluxRecorder_Impl::ContributionList;

    /** This function ... */
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
