/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OLIGOWAVELENGTHGRID_HPP
#define OLIGOWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** OligoWavelengthGrid is a subclass of the DisjointWavelengthGrid class representing the
    wavelength grid in an oligochromatic simulation. The class constructs a set of distinct
    nonadjacent wavelength bins with a half bin width given by 1/1000 of the shortest wavelength in
    the list. Refer to the DisjointWavelengthGrid class for more details.

    The list of wavelengths is specified through the constructor, because the oligo wavelength grid
    is created programmatically rather than directly configured by the user. The order of the
    specified wavelengths is not important; they will be sorted anyway. */
class OligoWavelengthGrid : public DisjointWavelengthGrid
{

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically to create an oligochromatic wavelength
        grid. Before the constructor returns, the newly created object is hooked up as a child to
        the specified parent in the simulation hierarchy (so it will automatically be deleted), and
        its setup() function has been called. */
    explicit OligoWavelengthGrid(SimulationItem* parent, const vector<double>& wavelengths);

protected:
    /** This function sets the wavelength grid to the wavelengths specified in the constructor. */
    void setupSelfBefore() override;

    //====================== Data members =====================

private:
    vector<double> _wavelengths;
};

//////////////////////////////////////////////////////////////////////

#endif
