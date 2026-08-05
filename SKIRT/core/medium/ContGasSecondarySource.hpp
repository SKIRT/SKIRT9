/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONTGASSECONDARYSOURCE_HPP
#define CONTGASSECONDARYSOURCE_HPP

#include "Array.hpp"
#include "SecondarySource.hpp"
class Configuration;
class EmittingGasMix;
class MediumSystem;
class PhotonPacket;
class Random;

//////////////////////////////////////////////////////////////////////

/** ContGasSecondarySource is a helper class that launches secondary emission photon packets from
    continuum emission for a given gas medium component. An instance of this class should be
    constructed only if secondary emission from gas is enabled for the simulation.

    For more information on the operation of this class, see the SecondarySourceSystem class. */
class ContGasSecondarySource : public SecondarySource
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a ContGasSecondarySource instance for the medium component with
        the specified index. This medium component must have a material mix that derives from the
        abstract EmittingGasMix class. Before the constructor returns, the newly created object is
        hooked up as a child to the specified parent in the simulation hierarchy (so it will
        automatically be deleted). */
    explicit ContGasSecondarySource(SimulationItem* parent, int h);

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
        can and does calculate the emission spectrum for a given cell only once. It preserves the
        relevant information from one invocation of the function to the next in a helper object
        allocated with thread-local storage scope. As a result, memory requirements are limited to
        storing the information for only a single cell per execution thread, and the calculation is
        still performed only once per cell.

        Once the emission spectrum for the current cell is known, the function randomly generates a
        wavelength either from this emission spectrum or from the configured bias wavelength
        distribution, adjusting the launch weight with the proper bias factor. It then generates a
        random position uniformly within the spatial cell. The function obtains the bulk velocity
        of the cell for application of the appropriate macroscopic Doppler shift; other than this
        the emission is assumed to be unpolarized and isotropic in the comoving frame.

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
    Array _Lv;           // the relative bolometric luminosity of each spatial cell (normalized to unity)
    Array _Wv;           // the relative launch weight for each spatial cell (normalized to unity)
    vector<size_t> _Iv;  // first history index allocated to each spatial cell (with extra entry at the end)
};

////////////////////////////////////////////////////////////////

#endif
