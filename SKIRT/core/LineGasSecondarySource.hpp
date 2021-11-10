/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINEGASSECONDARYSOURCE_HPP
#define LINEGASSECONDARYSOURCE_HPP

#include "Array.hpp"
#include "SecondarySource.hpp"
class Configuration;
class EmittingGasMix;
class MediumSystem;
class PhotonPacket;
class Random;

//////////////////////////////////////////////////////////////////////

/** LineGasSecondarySource is a helper class that launches secondary emission photon packets from
    emission lines for a given gas medium component. An instance of this class should be
    constructed only if secondary emission from gas is enabled for the simulation.

    For more information on the operation of this class, see the SecondarySourceSystem class. */
class LineGasSecondarySource : public SecondarySource
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a LineGasSecondarySource instance for the medium component with
        the specified index. This medium component must have a material mix that derives from the
        abstract EmittingGasMix class. Before the constructor returns, the newly created object is
        hooked up as a child to the specified parent in the simulation hierarchy (so it will
        automatically be deleted). */
    explicit LineGasSecondarySource(SimulationItem* parent, int h);

    //======================== Other Functions =======================

public:
    /** This function calculates and stores the bolometric luminosities for this source in each
        spatial cell of the simulation, and returns the total bolometric luminosity. */
    double prepareLuminosities() override;

    /** This function prepares the mapping of history indices to spatial cells, given the range of
        history indices allocated to this source. */
    void preparePacketMap(size_t firstIndex, size_t numIndices) override;

    /** This function causes the photon packet \em pp to be launched from one of the cells in the
        spatial grid using the given history index. Because photon packets launched from a given
        spatial cell are usually handled consecutively by the same execution thread, this function
        can and does calculate the line luminosities for a given cell only once. It preserves the
        relevant information from one invocation of the function to the next in a helper object
        allocated with thread-local storage scope. As a result, memory requirements are limited to
        storing the information for only a single cell per execution thread, and the calculation is
        still performed only once per cell.

        Once the line luminosities for the current cell are known, the function randomly selects
        one of the central line wavelengths from the discrete distribution represented by the line
        luminosities and/or from a uniform distribution as governed by the configured wavelength
        bias fraction (the configured bias wavelength distribution is ignored), adjusting the
        launch weight appropriately.

        The selected central line wavelength is adjusted with a random Doppler shift corresponding
        to the thermal motion of the particle emitting that line. The particle velocity is drawn
        from a Maxwell-Boltzman distribution parametrized by the temperature of the medium in the
        emitting cell (obtained from the medium state) and the mass of the emitting particle
        (obtained from the emitting material mix). If the temperature is zero or unavailable, or if
        the particle mass is zero, no thermal Doppler shift is added. The function further obtains
        the bulk velocity of the cell and applies the corresponding macroscopic Doppler shift.
        Other than this, the emission is assumed to be unpolarized and isotropic in the comoving
        frame.

        The function then generates a random position uniformly within the spatial cell.

        Finally, the function actually initializes the photon packet with this information. */
    void launch(PhotonPacket* pp, size_t historyIndex, double L) const override;

    //======================== Data Members ========================

private:
    // initialized by constructor
    int _h{0};  // medium system index of medium component being handled

    // initialized by prepareLuminosities() and preparePacketMap()
    Configuration* _config{nullptr};
    MediumSystem* _ms{nullptr};
    const EmittingGasMix* _mix{nullptr};  // material mix of medium component being handled
    Random* _random{nullptr};

    Array _centers;               // the central line wavelengths published by the material mix
    int _numLines{0};             // the number of lines published by the material mix
    bool _hasTemperature{false};  // true if the component's medium state has a temperature
    Array _masses;                // the particle masses published by the material mix (only if hasTemperature is true)

    Array _Lv;           // the relative bolometric luminosity of each spatial cell (normalized to unity)
    Array _Wv;           // the relative launch weight for each spatial cell (normalized to unity)
    vector<size_t> _Iv;  // first history index allocated to each spatial cell (with extra entry at the end)
};

////////////////////////////////////////////////////////////////

#endif
