/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CompositeWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include <iostream>

////////////////////////////////////////////////////////////////////

namespace
{
    /*  SegmentedGrid is a helper class for compositing wavelength grids.

        The wavelength grid being handled is represented by a list of bin border points and a
        corresponding list of characteristic wavelengths. A zero characteristic wavelength
        indicates a segment that is not part of the grid, i.e. that lies between two non-adjacent
        bins or beyond the outer grid borders. The two lists always have the same size. The borders
        are listed in strictly increasing order of wavelength (i.e. there are no duplicates) and
        all nonzero characteristic wavelengths lie within their corresponding bin.

        To facilitate the compositing algorithm, the first border point is always zero and the last
        border point is positive infinity. Both these border points have an associated
        characteristic wavelength of zero, indicating that the corresponding segments lie outside
        of the grid.

        The algorithm uses a special (negative) "invalid" value for characteristic wavelengths that
        need to be calculated at a later stage. */
    class SegmentedGrid
    {
    private:
        // initialize the grid to a single segment covering the full possible wavelength range
        vector<double> _borderv = {0., std::numeric_limits<double>::infinity()};
        vector<double> _characv = {0., 0.};

    public:
        // return constant references to our grid data
        const vector<double>& borderv() { return _borderv; }
        const vector<double>& characv() { return _characv; }

        // composite the given wavelength grid into the receiving segmented grid
        void add(const DisjointWavelengthGrid* wlg)
        {
            // loop over the bins of the incoming grid
            int numBins = wlg->numBins();
            for (int ell = 0; ell != numBins; ++ell)
            {
                double chr = wlg->wavelength(ell);
                double lef = wlg->leftBorder(ell);
                double rig = wlg->rightBorder(ell);

                // get the segment indices for the left and right borders
                int indexLef = NR::locate(_borderv, lef);
                int indexRig = NR::locate(_borderv, rig);

                // if both are fully inside the same empty segment, insert the input bin into that segment
                if (indexLef == indexRig && !_characv[indexLef])
                {
                    _borderv.insert(begin(_borderv) + indexLef + 1, rig);
                    _borderv.insert(begin(_borderv) + indexLef + 1, lef);
                    _characv.insert(begin(_characv) + indexLef + 1, 0.);
                    _characv.insert(begin(_characv) + indexLef + 1, chr);
                }
            }
        }

        void print()
        {
            for (size_t k = 0; k != _borderv.size(); ++k)
                std::cout << _borderv[k] * 1e6 << ':' << _characv[k] * 1e6 << "  ";
            std::cout << std::endl;
        }
    };
}

////////////////////////////////////////////////////////////////////

void CompositeWavelengthGrid::setupSelfAfter()
{
    // initialize the result grid to an empty grid
    SegmentedGrid grid;

    // composite the child grids into the result grid, in order of occurrence
    for (auto wlg : _wavelengthGrids)
    {
        grid.add(wlg);
        grid.print();
    }

    // pass the result grid to the base class
    setWavelengthSegments(grid.borderv(), grid.characv());

    // invoke base class at the end because it verifies that the wavelength grid has been initialized
    DisjointWavelengthGrid::setupSelfAfter();
}

////////////////////////////////////////////////////////////////////
