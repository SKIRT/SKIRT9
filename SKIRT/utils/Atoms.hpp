/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ATOMS_HPP
#define ATOMS_HPP

#include "Constants.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include <map>
#include <regex>

////////////////////////////////////////////////////////////////////

namespace
{
    static const std::map<string, short> atomMap = {
        {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},   {"N", 7},  {"O", 8},
        {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16},
        {"Cl", 17}, {"Ar", 18}, {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}};

    static const vector<double> masses = {1.0079,  4.0026, 6.941,   9.01218, 10.81,   12.011,  14.0067, 15.9994,
                                          18.9984, 20.179, 22.9898, 24.305,  26.9815, 28.0855, 30.9738, 32.06,
                                          35.453,  39.948, 39.0983, 40.08,   44.9559, 47.9,    50.9415, 51.996,
                                          54.938,  55.847, 58.9332, 58.7,    63.546,  65.38};
}

/** This class provides utility functions related to the treatment of atomic and ionic species. */
class Atoms
{
public:
    /** This function returns the atomic number of the specified elements. */
    static inline short atomToZ(string element) { return atomMap.at(element); }

    /** This function returns the atomic number of the specified element. */
    static inline double umass(int Z) { return masses[Z - 1]; }

    static inline double mass(short Z) { return masses[Z - 1] * Constants::amu(); }

    static inline std::pair<short, short> parseIon(string ion)
    {
        // read ions
        ion = StringUtils::squeeze(ion);
        // split ion string into element symbol and ionization number
        std::regex pattern("([A-Za-z]+)(\\+?)([0-9]*)");
        std::smatch match;

        if (!std::regex_match(ion, match, pattern)) throw FATALERROR("Could not parse ion format: " + ion);

        int Z, N;
        string element = match[1].str();
        string plus = match[2].str();
        string number = match[3].str();

        Z = atomMap.at(element);

        // enforce format to avoid ambiguities
        if (number.empty())
        {
            if (plus.empty())
                N = Z;  // eg. Fe
            else
                N = Z - 1;  // eg. Fe+
        }
        else
        {
            if (!plus.empty())
                N = Z - StringUtils::toInt(number);  // eg. Fe+14
            else
                throw FATALERROR("Ion format should contain a '+' eg. Fe+14 for ion: " + ion);  // eg. Fe14
        }

        if (N < 0 || N > Z) throw FATALERROR("Invalid electron number for ion: " + ion);

        return std::make_pair(Z, N);
    }
};

#endif
