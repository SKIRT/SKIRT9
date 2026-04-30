/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ATOMS_HPP
#define ATOMS_HPP

#include "Basics.hpp"

/** This static class provides utility functions related to the treatment of atomic and ionic species. */
class Atoms
{
public:
    /** This function returns the atomic number of the specified elements. */
    static short atomToZ(string element);

    /** This function returns the mass of the specified atomic number in kg. */
    static double mass(short Z);

    /** This function returns a pair of the atomic number and number of electrons (Z,N) of
    the specified ion string. The ion string can have the following formats, as listed
    in the table below. Any other format will throw an error. <TABLE>
    <TR><TD><B>Format</B></TD> <TD><B>Returns</B></TD></TR>
    <TR><TD><TT>'Z'</TT></TD> <TD>\f$(Z,N)=(Z,Z)\f$</TD></TR>
    <TR><TD><TT>'Z+'</TT></TD> <TD>\f$(Z,N)=(Z,Z-1)\f$</TD></TR>
    <TR><TD><TT>'Z+N'</TT></TD> <TD>\f$(Z,N)=(Z,N)\f$</TD></TR>
    </TABLE> */
    static std::pair<short, short> parseIon(string ion);
};

#endif
