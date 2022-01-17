/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CompositeWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

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
    */
    class SegmentedGrid
    {
    private:
        // log or linear scale, set by constructor
        bool _logScale{false};

        // initialize the grid to a single segment covering the full possible wavelength range
        vector<double> _borderv = {0., std::numeric_limits<double>::infinity()};
        vector<double> _characv = {0., 0.};

        // relative tolerances
        constexpr static double TOLMIN = 1. - 1e-9;
        constexpr static double TOLMAX = 1. + 1e-9;

        // return the characteristic wavelength for given border wavelengths
        double characteristic(double b1, double b2)
        {
            if (_logScale)
                return sqrt(b1 * b2);
            else
                return 0.5 * (b1 + b2);
        }

    public:
        // construct "empty" grid with given scale
        SegmentedGrid(bool logScale) : _logScale{logScale} {}

        // return constant references to our grid data
        const vector<double>& borderv() { return _borderv; }
        const vector<double>& characv() { return _characv; }

        // returns true if the given wavelength is at a segment border within a small tolerance, and false otherwise
        bool isBorder(double wave) const
        {
            int index = NR::locate(_borderv, wave);
            if (wave < _borderv[index] * TOLMAX) return true;
            if (wave > _borderv[index + 1] * TOLMIN) return true;
            return false;
        }

        // returns the index of the segment containing the given wavelength, or if the given wavelength
        // is at a segment border within a small tolerance, the index of that border
        int getIndex(double wave) const
        {
            int index = NR::locate(_borderv, wave);
            if (wave > _borderv[index + 1] * TOLMIN) return index + 1;
            return index;
        }

        // splits the segment containing the given wavelength at that wavelength
        // unless the given wavelength already is at a border
        void splitSegment(double split)
        {
            // only split the segment if the split point is not at a border
            if (!isBorder(split))
            {
                // get the index of the segment to be split
                int index = getIndex(split);

                // determine the characteristic wavelengths for the new segments
                double leftcharac = 0.;
                double rightcharac = 0.;
                double old = _characv[index];
                if (old > 0.)
                {
                    leftcharac = (old < split * TOLMIN) ? old : characteristic(_borderv[index], split);
                    rightcharac = (old > split * TOLMAX) ? old : characteristic(split, _borderv[index + 1]);
                }

                // insert a new segment at the given wavelength
                _characv[index] = leftcharac;
                _borderv.insert(begin(_borderv) + index + 1, split);
                _characv.insert(begin(_characv) + index + 1, rightcharac);
            }
        }

        // replace any segments within the range of the new bin by the new bin;
        // this function assumes that the borders of the new bin are already present
        void replaceSegments(double left, double charac, double right)
        {
            int leftIndex = getIndex(left);
            int rightIndex = getIndex(right);
            if (rightIndex > leftIndex)
            {
                _borderv.erase(begin(_borderv) + leftIndex + 1, begin(_borderv) + rightIndex);
                _characv.erase(begin(_characv) + leftIndex + 1, begin(_characv) + rightIndex);
                _characv[leftIndex] = charac;
            }
        }

        // composite the given wavelength grid into the receiving segmented grid
        void add(const DisjointWavelengthGrid* wlg)
        {
            // loop over the bins of the incoming grid
            int numBins = wlg->numBins();
            for (int ell = 0; ell != numBins; ++ell)
            {
                double charac = wlg->wavelength(ell);
                double left = wlg->leftBorder(ell);
                double right = wlg->rightBorder(ell);

                // if the bin borders are not already present in the grid,
                // insert them by splitting the containing segment
                splitSegment(left);
                splitSegment(right);

                // replace any segments within range by the new bin
                replaceSegments(left, charac, right);
            }
        }
    };
}

////////////////////////////////////////////////////////////////////

void CompositeWavelengthGrid::setupSelfAfter()
{
    // initialize the result grid to an empty grid
    SegmentedGrid grid(log());

    // composite the child grids into the result grid, in order of occurrence
    for (auto wlg : _wavelengthGrids)
    {
        grid.add(wlg);
    }

    // pass the result grid to the base class
    setWavelengthSegments(grid.borderv(), grid.characv());

    // invoke base class at the end because it verifies that the wavelength grid has been initialized
    DisjointWavelengthGrid::setupSelfAfter();
}

////////////////////////////////////////////////////////////////////
