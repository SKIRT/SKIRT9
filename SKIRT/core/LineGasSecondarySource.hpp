/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINEGASSECONDARYSOURCE_HPP
#define LINEGASSECONDARYSOURCE_HPP

#include "Array.hpp"
#include "SecondarySource.hpp"
class Configuration;
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
        spatial grid using the given history index.

        TO DO XXX.

        Finally, the function actually initializes the photon packet with this information. */
    void launch(PhotonPacket* pp, size_t historyIndex, double L) const override;

    //======================== Data Members ========================

private:
    // initialized by constructor
    int _h{0};  // medium system index of medium component being handled

    // initialized by prepareLuminosities() and preparePacketMap()
    Configuration* _config{nullptr};
    MediumSystem* _ms{nullptr};
    Random* _random{nullptr};
    double _L{0};        // the total bolometric luminosity of all spatial cells
    Array _Lv;           // the relative bolometric luminosity of each spatial cell (normalized to unity)
    Array _Wv;           // the relative launch weight for each spatial cell (normalized to unity)
    vector<size_t> _Iv;  // first history index allocated to each spatial cell (with extra entry at the end)
};

////////////////////////////////////////////////////////////////

#endif
