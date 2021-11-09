/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTSECONDARYSOURCE_HPP
#define DUSTSECONDARYSOURCE_HPP

#include "Array.hpp"
#include "SecondarySource.hpp"
class Configuration;
class MediumSystem;
class PhotonPacket;
class Random;

//////////////////////////////////////////////////////////////////////

/** DustSecondarySource is a helper class that launches secondary emission photon packets from dust
    media. It handles the aggregated thermal emission from all dust components in the simulation.
    An instance of this class should be constructed only if the simulation configuration includes
    one or more dust medium components and secondary emission from dust is enabled.

    For more information on the operation of this class, see the SecondarySourceSystem class. */
class DustSecondarySource : public SecondarySource
{
    //============= Construction - Setup - Destruction =============

public:
    using SecondarySource::SecondarySource;

    //======================== Other Functions =======================

public:
    /** This function calculates and stores the bolometric dust luminosities in each spatial cell
        of the simulation, and returns the total bolometric dust luminosity. */
    double prepareLuminosities() override;

    /** This function prepares the mapping of history indices to spatial cells, given the range of
        history indices allocated to this source. */
    void preparePacketMap(size_t firstIndex, size_t numIndices) override;

    /** This function causes the photon packet \em pp to be launched from one of the cells in the
        spatial grid using the given history index.

        Before it can randomly emit photon packets from a spatial cell, this function must
        calculate the normalized regular and cumulative dust emission spectrum for the cell from
        the radiation field and the dust properties held by the medium system, using the configured
        emission calculator (LTE or NLTE). Also, it must obtain the average bulk velocity of the
        material in the cell from the medium system. Especially the NLTE emission spectrum
        calculations can be very time-consuming, so it is important to perform them only once for
        each cell. On the other hand, pre-calculating and storing the spectra for all cells in
        advance requires a potentially large amount of memory (proportional to both the number of
        cells and to the number of wavelengths on which the emission is discretized). As described
        in the SecondarySourceSystem class header, photon packets launched from a given spatial
        cell are usually handled consecutively by the same execution thread, allowing this function
        to remember the information calculated for the "current" cell from one invocation to the
        next in a helper object allocated with thread-local storage scope. As a result, memory
        requirements are limited to storing the information for only a single cell per execution
        thread, and the calculation is still performed only once per cell.

        Once the emission spectrum for the current cell is known, the function randomly generates a
        wavelength either from this emission spectrum or from the configured bias wavelength
        distribution, adjusting the launch weight with the proper bias factor. It then generates a
        random position uniformly within the spatial cell.

        If the simulation includes aligned spheroidal grains, the function calculates the
        polarization profile for the emitted photon packet and draws a random direction from the
        corresponding anisotropic phase function. Otherwise, the emission is assumed to be
        unpolarized and isotropic, and a random direction is determined uniformly on the unit
        sphere.

        Finally, the function actually initializes the photon packet with this information. */
    void launch(PhotonPacket* pp, size_t historyIndex, double L) const override;

    //======================== Data Members ========================

private:
    // initialized by prepareLuminosities() and preparePacketMap()
    Configuration* _config{nullptr};
    MediumSystem* _ms{nullptr};
    Random* _random{nullptr};
    Array _Lv;           // the relative bolometric luminosity of each spatial cell (normalized to unity)
    Array _Wv;           // the relative launch weight for each spatial cell (normalized to unity)
    vector<int> _nv;     // the library entry index corresponding to each spatial cell (i.e. map from cells to entries)
    vector<int> _mv;     // the spatial cell indices sorted so that cells belonging to the same entry are consecutive
    vector<size_t> _Iv;  // first history index allocated to each spatial cell (with extra entry at the end)
};

////////////////////////////////////////////////////////////////

#endif
