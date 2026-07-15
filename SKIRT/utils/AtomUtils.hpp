/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ATOMUTILS_HPP
#define ATOMUTILS_HPP

#include "Basics.hpp"

/** This static class provides utility functions related to the treatment of atomic and ionic
    species. */
class AtomUtils final
{
public:
    /** This function returns the atomic number of the specified elements. */
    static short atomToZ(string element);

    /** This function returns the mass of the specified atomic number in SI units. */
    static double mass(short Z);

    static int ionIndex(int Z, int N);

    /** This function returns a pair of the atomic number and number of electrons (Z,N) of the
        specified ion string. The ion string can have the following formats, as listed in the table
        below. Any other format will throw an error. <TABLE>
        <TR><TD><B>Format</B></TD><TD><B>Returns</B></TD></TR>
        <TR><TD><TT>'Z'</TT></TD><TD>\f$(Z,N)=(Z,Z)\f$</TD></TR>
        <TR><TD><TT>'Z+'</TT></TD><TD>\f$(Z,N)=(Z,Z-1)\f$</TD></TR>
        <TR><TD><TT>'Z+I'</TT></TD><TD>\f$(Z,N)=(Z,Z-I)\f$</TD></TR> </TABLE> The number of
        electrons, N, can range from 0 to Z, any other value will result in an error. */
    static std::pair<short, short> parseIon(string ion);
};

#endif
