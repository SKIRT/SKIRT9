/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayAtomicGasMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "Range.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // support the elements with atomic number up to 30
    // 1   2   3   4  5  6  7  8  9  10  11  12  13  14 15 16  17  18 19  20  21  22 23  24  25  26  27  28  29  30
    // H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
    constexpr size_t numAtoms = 30;

    // default abundancies from Table 2 of Anders & Grevesse (1989), the default abundance table in Xspec
    constexpr std::initializer_list<double> defaultAbundancies = {
        1.00E+00, 9.77E-02, 1.45E-11, 1.41E-11, 3.98E-10, 3.63E-04, 1.12E-04, 8.51E-04, 3.63E-08, 1.23E-04,
        2.14E-06, 3.80E-05, 2.95E-06, 3.55E-05, 2.82E-07, 1.62E-05, 3.16E-07, 3.63E-06, 1.32E-07, 2.29E-06,
        1.26E-09, 9.77E-08, 1.00E-08, 4.68E-07, 2.45E-07, 4.68E-05, 8.32E-08, 1.78E-06, 1.62E-08, 3.98E-08};
    static_assert(numAtoms == defaultAbundancies.size(), "Incorrect number of default abundancies");

    // atomic masses (amu) from EADL (Evaluated Atomic Data Library of the Lawrence Livermore National Laboratory)
    constexpr std::initializer_list<double> masses = {
        1.00790e+00, 4.00260e+00, 6.94100e+00, 9.01218e+00, 1.08100e+01, 1.20110e+01, 1.40067e+01, 1.59994e+01,
        1.89984e+01, 2.01790e+01, 2.29898e+01, 2.43050e+01, 2.69815e+01, 2.80855e+01, 3.09738e+01, 3.20600e+01,
        3.54530e+01, 3.99480e+01, 3.90983e+01, 4.00800e+01, 4.49559e+01, 4.79000e+01, 5.09415e+01, 5.19960e+01,
        5.49380e+01, 5.58470e+01, 5.89332e+01, 5.87000e+01, 6.35460e+01, 6.53800e+01};
    static_assert(numAtoms == masses.size(), "Incorrect number of masses");

    // photo-absorption cross sections from Table 1 in Verner, D. A., & Yakovlev, D. G. 1995, A&AS, 109, 125
    // downloaded from www.pa.uky.edu/~verner/photo.html and removed lines for ionized atoms
    struct CrossSectionParams
    {
        short Z;        // atomic number
        short n;        // principal quantum number of the shell
        short l;        // orbital quantum number of the subshell
        double Eth;     // subshell ionization threshold energy (eV)
        double E0;      // fit parameter (eV)
        double sigma0;  // fit parameter (Mb = 10^-22 m^2)
        double ya;      // fit parameter (1)
        double P;       // fit parameter (1)
        double yw;      // fit parameter (1)
    };
    constexpr std::initializer_list<CrossSectionParams> crossSectionParams = {
        {1, 1, 0, 0.1360E+02, 0.4298E+00, 0.5475E+05, 0.3288E+02, 0.2963E+01, 0.0000E+00},
        {2, 1, 0, 0.2459E+02, 0.5996E+01, 0.4470E+04, 0.2199E+01, 0.6098E+01, 0.0000E+00},
        {3, 2, 0, 0.5392E+01, 0.3466E+01, 0.4774E+02, 0.2035E+02, 0.4423E+01, 0.0000E+00},
        {3, 1, 0, 0.6439E+02, 0.2740E+02, 0.1564E+03, 0.3382E+02, 0.1490E+01, 0.0000E+00},
        {4, 2, 0, 0.9323E+01, 0.3427E+01, 0.1423E+04, 0.2708E+01, 0.1064E+02, 0.0000E+00},
        {4, 1, 0, 0.1193E+03, 0.4093E+02, 0.1306E+03, 0.1212E+03, 0.1348E+01, 0.0000E+00},
        {5, 2, 1, 0.8298E+01, 0.6658E+01, 0.3643E+03, 0.9500E+01, 0.5469E+01, 0.5608E+00},
        {5, 2, 0, 0.1405E+02, 0.6495E+01, 0.1525E+05, 0.1352E+01, 0.1218E+02, 0.0000E+00},
        {5, 1, 0, 0.1940E+03, 0.6155E+02, 0.9698E+02, 0.7354E+02, 0.1438E+01, 0.0000E+00},
        {6, 2, 1, 0.1126E+02, 0.9435E+01, 0.1152E+04, 0.5687E+01, 0.6336E+01, 0.4474E+00},
        {6, 2, 0, 0.1939E+02, 0.1026E+02, 0.4564E+04, 0.1568E+01, 0.1085E+02, 0.0000E+00},
        {6, 1, 0, 0.2910E+03, 0.8655E+02, 0.7421E+02, 0.5498E+02, 0.1503E+01, 0.0000E+00},
        {7, 2, 1, 0.1453E+02, 0.1164E+02, 0.1029E+05, 0.2361E+01, 0.8821E+01, 0.4239E+00},
        {7, 2, 0, 0.2541E+02, 0.1482E+02, 0.7722E+03, 0.2306E+01, 0.9139E+01, 0.0000E+00},
        {7, 1, 0, 0.4048E+03, 0.1270E+03, 0.4748E+02, 0.1380E+03, 0.1252E+01, 0.0000E+00},
        {8, 2, 1, 0.1362E+02, 0.1391E+02, 0.1220E+06, 0.1364E+01, 0.1140E+02, 0.4103E+00},
        {8, 2, 0, 0.2848E+02, 0.1994E+02, 0.2415E+03, 0.3241E+01, 0.8037E+01, 0.0000E+00},
        {8, 1, 0, 0.5380E+03, 0.1774E+03, 0.3237E+02, 0.3812E+03, 0.1083E+01, 0.0000E+00},
        {9, 2, 1, 0.1742E+02, 0.1658E+02, 0.2775E+06, 0.1242E+01, 0.1249E+02, 0.3857E+00},
        {9, 2, 0, 0.3786E+02, 0.2568E+02, 0.1097E+03, 0.4297E+01, 0.7303E+01, 0.0000E+00},
        {9, 1, 0, 0.6940E+03, 0.2390E+03, 0.2295E+02, 0.1257E+04, 0.9638E+00, 0.0000E+00},
        {10, 2, 1, 0.2156E+02, 0.2000E+02, 0.1691E+05, 0.2442E+01, 0.1043E+02, 0.3345E+00},
        {10, 2, 0, 0.4847E+02, 0.3204E+02, 0.5615E+02, 0.5808E+01, 0.6678E+01, 0.0000E+00},
        {10, 1, 0, 0.8701E+03, 0.3144E+03, 0.1664E+02, 0.2042E+06, 0.8450E+00, 0.0000E+00},
        {11, 3, 0, 0.5139E+01, 0.5968E+01, 0.1460E+01, 0.2557E+08, 0.3789E+01, 0.0000E+00},
        {11, 2, 1, 0.3814E+02, 0.3655E+02, 0.2486E+03, 0.3222E+03, 0.3570E+01, 0.1465E+00},
        {11, 2, 0, 0.7084E+02, 0.4537E+02, 0.1142E+02, 0.2395E+03, 0.3380E+01, 0.0000E+00},
        {11, 1, 0, 0.1079E+04, 0.4216E+03, 0.1119E+02, 0.5642E+08, 0.7736E+00, 0.0000E+00},
        {12, 3, 0, 0.7646E+01, 0.9393E+01, 0.3034E+01, 0.2625E+08, 0.3923E+01, 0.0000E+00},
        {12, 2, 1, 0.5490E+02, 0.4937E+02, 0.2023E+03, 0.1079E+05, 0.2960E+01, 0.1463E-01},
        {12, 2, 0, 0.9400E+02, 0.4587E+02, 0.1671E+02, 0.2389E+02, 0.4742E+01, 0.0000E+00},
        {12, 1, 0, 0.1311E+04, 0.2711E+03, 0.3561E+02, 0.2374E+02, 0.1952E+01, 0.0000E+00},
        {13, 3, 1, 0.5986E+01, 0.1860E+02, 0.1828E+03, 0.2797E+01, 0.1084E+02, 0.3076E+00},
        {13, 3, 0, 0.1133E+02, 0.1204E+02, 0.5384E+01, 0.4341E+03, 0.4088E+01, 0.0000E+00},
        {13, 2, 1, 0.8040E+02, 0.6445E+02, 0.1735E+03, 0.1131E+05, 0.2762E+01, 0.2337E-01},
        {13, 2, 0, 0.1256E+03, 0.5594E+02, 0.1425E+02, 0.3094E+02, 0.4399E+01, 0.0000E+00},
        {13, 1, 0, 0.1567E+04, 0.3670E+03, 0.2206E+02, 0.4405E+02, 0.1588E+01, 0.0000E+00},
        {14, 3, 1, 0.8152E+01, 0.2212E+02, 0.1845E+03, 0.3849E+01, 0.9721E+01, 0.2921E+00},
        {14, 3, 0, 0.1517E+02, 0.1413E+02, 0.1166E+02, 0.2288E+02, 0.5334E+01, 0.0000E+00},
        {14, 2, 1, 0.1060E+03, 0.7808E+02, 0.1532E+03, 0.5765E+07, 0.2639E+01, 0.2774E-03},
        {14, 2, 0, 0.1560E+03, 0.7017E+02, 0.1166E+02, 0.4742E+02, 0.3933E+01, 0.0000E+00},
        {14, 1, 0, 0.1846E+04, 0.5322E+03, 0.1184E+02, 0.2580E+03, 0.1102E+01, 0.0000E+00},
        {15, 3, 1, 0.1049E+02, 0.2580E+02, 0.9925E+02, 0.6712E+01, 0.8516E+01, 0.2765E+00},
        {15, 3, 0, 0.2017E+02, 0.1658E+02, 0.1125E+02, 0.2613E+02, 0.5205E+01, 0.0000E+00},
        {15, 2, 1, 0.1400E+03, 0.8812E+02, 0.1512E+03, 0.2230E+04, 0.2795E+01, 0.3422E-03},
        {15, 2, 0, 0.1940E+03, 0.8632E+02, 0.9931E+01, 0.6594E+02, 0.3617E+01, 0.0000E+00},
        {15, 1, 0, 0.2154E+04, 0.6472E+03, 0.9167E+01, 0.2562E+03, 0.1063E+01, 0.0000E+00},
        {16, 3, 1, 0.1036E+02, 0.2975E+02, 0.5644E+02, 0.1321E+02, 0.7513E+01, 0.2621E+00},
        {16, 3, 0, 0.2130E+02, 0.1916E+02, 0.1003E+02, 0.3296E+02, 0.5038E+01, 0.0000E+00},
        {16, 2, 1, 0.1700E+03, 0.9152E+02, 0.1883E+03, 0.7193E+02, 0.3633E+01, 0.2485E+00},
        {16, 2, 0, 0.2350E+03, 0.1047E+03, 0.8520E+01, 0.9469E+02, 0.3346E+01, 0.0000E+00},
        {16, 1, 0, 0.2477E+04, 0.8114E+03, 0.6649E+01, 0.3734E+04, 0.8646E+00, 0.0000E+00},
        {17, 3, 1, 0.1297E+02, 0.3398E+02, 0.4539E+02, 0.2232E+02, 0.6896E+01, 0.2479E+00},
        {17, 3, 0, 0.2531E+02, 0.2231E+02, 0.6628E+01, 0.1843E+03, 0.4196E+01, 0.0000E+00},
        {17, 2, 1, 0.2090E+03, 0.8004E+02, 0.3053E+03, 0.3498E+02, 0.4457E+01, 0.2017E+00},
        {17, 2, 0, 0.2780E+03, 0.1092E+03, 0.1059E+02, 0.2491E+02, 0.4205E+01, 0.0000E+00},
        {17, 1, 0, 0.2830E+04, 0.9700E+03, 0.5255E+01, 0.1856E+07, 0.7888E+00, 0.0000E+00},
        {18, 3, 1, 0.1576E+02, 0.3854E+02, 0.4872E+02, 0.2640E+02, 0.6662E+01, 0.2355E+00},
        {18, 3, 0, 0.2892E+02, 0.2525E+02, 0.6394E+01, 0.1700E+03, 0.4223E+01, 0.0000E+00},
        {18, 2, 1, 0.2492E+03, 0.1647E+03, 0.8372E+02, 0.5452E+02, 0.3328E+01, 0.6270E+00},
        {18, 2, 0, 0.3260E+03, 0.1302E+03, 0.9185E+01, 0.2693E+02, 0.4021E+01, 0.0000E+00},
        {18, 1, 0, 0.3203E+04, 0.1135E+04, 0.4280E+01, 0.3285E+08, 0.7631E+00, 0.0000E+00},
        {19, 4, 0, 0.4341E+01, 0.3824E+01, 0.7363E+00, 0.2410E+08, 0.4427E+01, 0.2049E-03},
        {19, 3, 1, 0.2466E+02, 0.4138E+02, 0.2614E+02, 0.2143E+03, 0.5631E+01, 0.2437E+00},
        {19, 3, 0, 0.4080E+02, 0.2910E+02, 0.6377E+01, 0.2229E+04, 0.3587E+01, 0.0000E+00},
        {19, 2, 1, 0.3014E+03, 0.2666E+03, 0.3107E+02, 0.7187E+07, 0.2067E+01, 0.5274E+00},
        {19, 2, 0, 0.3843E+03, 0.1602E+03, 0.6389E+01, 0.1044E+03, 0.3159E+01, 0.0000E+00},
        {19, 1, 0, 0.3614E+04, 0.1171E+04, 0.4540E+01, 0.6165E+04, 0.8392E+00, 0.0000E+00},
        {20, 4, 0, 0.6113E+01, 0.7366E+01, 0.2373E+01, 0.2082E+03, 0.4841E+01, 0.5841E-03},
        {20, 3, 1, 0.3443E+02, 0.4487E+02, 0.9017E+02, 0.1465E+02, 0.7498E+01, 0.2754E+00},
        {20, 3, 0, 0.4830E+02, 0.3012E+02, 0.7227E+01, 0.1736E+03, 0.4165E+01, 0.0000E+00},
        {20, 2, 1, 0.3523E+03, 0.1529E+03, 0.1282E+03, 0.2217E+03, 0.3087E+01, 0.3343E-02},
        {20, 2, 0, 0.4425E+03, 0.1201E+03, 0.1010E+02, 0.2468E+02, 0.4592E+01, 0.0000E+00},
        {20, 1, 0, 0.4043E+04, 0.6947E+03, 0.1586E+02, 0.2563E+02, 0.1966E+01, 0.0000E+00},
        {21, 4, 0, 0.7342E+01, 0.8324E+01, 0.2252E+01, 0.4118E+03, 0.4588E+01, 0.1441E-03},
        {21, 3, 2, 0.8010E+01, 0.9733E+01, 0.2488E+03, 0.2066E+02, 0.1022E+02, 0.3104E+00},
        {21, 3, 1, 0.3360E+02, 0.4936E+02, 0.6176E+02, 0.2186E+02, 0.7081E+01, 0.2671E+00},
        {21, 3, 0, 0.5640E+02, 0.3284E+02, 0.6849E+01, 0.1709E+03, 0.4207E+01, 0.0000E+00},
        {21, 2, 1, 0.4054E+03, 0.1593E+03, 0.1473E+03, 0.6007E+02, 0.3635E+01, 0.3208E-02},
        {21, 2, 0, 0.5032E+03, 0.1180E+03, 0.1065E+02, 0.2172E+02, 0.4911E+01, 0.0000E+00},
        {21, 1, 0, 0.4494E+04, 0.7367E+03, 0.1554E+02, 0.2940E+02, 0.1937E+01, 0.0000E+00},
        {22, 4, 0, 0.6820E+01, 0.9184E+01, 0.2167E+01, 0.4297E+03, 0.4552E+01, 0.3612E-02},
        {22, 3, 2, 0.9940E+01, 0.1102E+02, 0.5478E+03, 0.1610E+02, 0.1096E+02, 0.2972E+00},
        {22, 3, 1, 0.4000E+02, 0.5424E+02, 0.5924E+02, 0.2151E+02, 0.7123E+01, 0.2645E+00},
        {22, 3, 0, 0.6500E+02, 0.3562E+02, 0.6517E+01, 0.1572E+03, 0.4268E+01, 0.0000E+00},
        {22, 2, 1, 0.4640E+03, 0.1814E+03, 0.1291E+03, 0.6629E+02, 0.3531E+01, 0.5519E-04},
        {22, 2, 0, 0.5690E+03, 0.1177E+03, 0.1068E+02, 0.2144E+02, 0.5085E+01, 0.0000E+00},
        {22, 1, 0, 0.4972E+04, 0.6826E+03, 0.2029E+02, 0.2415E+02, 0.2187E+01, 0.0000E+00},
        {23, 4, 0, 0.6740E+01, 0.1002E+02, 0.2059E+01, 0.3914E+03, 0.4565E+01, 0.5057E-03},
        {23, 3, 2, 0.1200E+02, 0.1146E+02, 0.3134E+04, 0.7037E+01, 0.1377E+02, 0.3417E+00},
        {23, 3, 1, 0.4700E+02, 0.5937E+02, 0.7599E+02, 0.1488E+02, 0.7586E+01, 0.2690E+00},
        {23, 3, 0, 0.7700E+02, 0.3668E+02, 0.6690E+01, 0.7759E+02, 0.4706E+01, 0.0000E+00},
        {23, 2, 1, 0.5270E+03, 0.2290E+03, 0.8513E+02, 0.5037E+03, 0.2767E+01, 0.9964E-03},
        {23, 2, 0, 0.6380E+03, 0.9688E+02, 0.1212E+02, 0.1905E+02, 0.5757E+01, 0.0000E+00},
        {23, 1, 0, 0.5475E+04, 0.6550E+03, 0.2420E+02, 0.2343E+02, 0.2326E+01, 0.0000E+00},
        {24, 4, 0, 0.6767E+01, 0.9636E+01, 0.6532E+00, 0.5232E+03, 0.4641E+01, 0.9332E-04},
        {24, 3, 2, 0.8660E+01, 0.7244E+01, 0.1485E+04, 0.9671E+01, 0.1575E+02, 0.7760E+00},
        {24, 3, 1, 0.4900E+02, 0.6567E+02, 0.5313E+02, 0.1981E+02, 0.7258E+01, 0.2601E+00},
        {24, 3, 0, 0.7900E+02, 0.4234E+02, 0.5602E+01, 0.1356E+03, 0.4374E+01, 0.0000E+00},
        {24, 2, 1, 0.5850E+03, 0.1865E+03, 0.1540E+03, 0.5560E+02, 0.3823E+01, 0.5785E-02},
        {24, 2, 0, 0.7030E+03, 0.1588E+03, 0.8555E+01, 0.2258E+02, 0.4789E+01, 0.0000E+00},
        {24, 1, 0, 0.5996E+04, 0.1035E+04, 0.1025E+02, 0.3343E+02, 0.1822E+01, 0.0000E+00},
        {25, 4, 0, 0.7434E+01, 0.1183E+02, 0.1549E+01, 0.2920E+07, 0.4113E+01, 0.3256E-01},
        {25, 3, 2, 0.1430E+02, 0.1311E+02, 0.1668E+05, 0.4497E+01, 0.1646E+02, 0.3881E+00},
        {25, 3, 1, 0.5940E+02, 0.7040E+02, 0.6697E+02, 0.1485E+02, 0.7643E+01, 0.2659E+00},
        {25, 3, 0, 0.9460E+02, 0.4114E+02, 0.6172E+01, 0.5928E+02, 0.4989E+01, 0.0000E+00},
        {25, 2, 1, 0.6554E+03, 0.2737E+03, 0.7498E+02, 0.3952E+03, 0.2793E+01, 0.4661E-01},
        {25, 2, 0, 0.7816E+03, 0.8316E+02, 0.1156E+02, 0.2187E+02, 0.6149E+01, 0.0000E+00},
        {25, 1, 0, 0.6550E+04, 0.8758E+03, 0.1592E+02, 0.3965E+02, 0.1947E+01, 0.0000E+00},
        {26, 4, 0, 0.7902E+01, 0.1277E+02, 0.1468E+01, 0.1116E+06, 0.4112E+01, 0.3238E-01},
        {26, 3, 2, 0.1470E+02, 0.1407E+02, 0.1850E+05, 0.4458E+01, 0.1691E+02, 0.4039E+00},
        {26, 3, 1, 0.6600E+02, 0.7630E+02, 0.6298E+02, 0.1479E+02, 0.7672E+01, 0.2646E+00},
        {26, 3, 0, 0.1040E+03, 0.4334E+02, 0.5921E+01, 0.5293E+02, 0.5129E+01, 0.0000E+00},
        {26, 2, 1, 0.7240E+03, 0.2948E+03, 0.7191E+02, 0.3219E+03, 0.2837E+01, 0.6314E-01},
        {26, 2, 0, 0.8570E+03, 0.5727E+02, 0.1076E+02, 0.2785E+02, 0.6635E+01, 0.0000E+00},
        {26, 1, 0, 0.7124E+04, 0.8044E+03, 0.2055E+02, 0.3633E+02, 0.2118E+01, 0.0000E+00},
        {27, 4, 0, 0.7864E+01, 0.1370E+02, 0.1555E+01, 0.7559E+03, 0.4337E+01, 0.3355E-01},
        {27, 3, 2, 0.1580E+02, 0.1581E+02, 0.4931E+04, 0.6607E+01, 0.1532E+02, 0.3676E+00},
        {27, 3, 1, 0.7300E+02, 0.8256E+02, 0.4587E+02, 0.1973E+02, 0.7331E+01, 0.2573E+00},
        {27, 3, 0, 0.1150E+03, 0.5097E+02, 0.4896E+01, 0.1198E+03, 0.4513E+01, 0.0000E+00},
        {27, 2, 1, 0.8000E+03, 0.2832E+03, 0.9075E+02, 0.7686E+02, 0.3416E+01, 0.4833E-04},
        {27, 2, 0, 0.9400E+03, 0.8888E+02, 0.1030E+02, 0.2797E+02, 0.5913E+01, 0.0000E+00},
        {27, 1, 0, 0.7725E+04, 0.6269E+03, 0.3582E+02, 0.3161E+02, 0.2476E+01, 0.0000E+00},
        {28, 4, 0, 0.7637E+01, 0.1468E+02, 0.1437E+01, 0.7411E+03, 0.4342E+01, 0.3908E-01},
        {28, 3, 2, 0.1700E+02, 0.6063E+01, 0.1186E+04, 0.6823E+01, 0.2223E+02, 0.6227E-02},
        {28, 3, 1, 0.8200E+02, 0.8896E+02, 0.4351E+02, 0.1942E+02, 0.7372E+01, 0.2566E+00},
        {28, 3, 0, 0.1250E+03, 0.5448E+02, 0.4611E+01, 0.1157E+03, 0.4548E+01, 0.0000E+00},
        {28, 2, 1, 0.8760E+03, 0.3043E+03, 0.8611E+02, 0.7868E+02, 0.3408E+01, 0.1680E-04},
        {28, 2, 0, 0.1024E+04, 0.1132E+03, 0.9424E+01, 0.2712E+02, 0.5643E+01, 0.0000E+00},
        {28, 1, 0, 0.8348E+04, 0.7366E+03, 0.2836E+02, 0.3622E+02, 0.2316E+01, 0.0000E+00},
        {29, 4, 0, 0.7726E+01, 0.1436E+02, 0.4681E+00, 0.2383E+04, 0.4224E+01, 0.3736E-01},
        {29, 3, 2, 0.1064E+02, 0.7279E+01, 0.1027E+04, 0.7988E+01, 0.2033E+02, 0.1582E+01},
        {29, 3, 1, 0.8300E+02, 0.9617E+02, 0.4275E+02, 0.1747E+02, 0.7555E+01, 0.2599E+00},
        {29, 3, 0, 0.1288E+03, 0.6198E+02, 0.4164E+01, 0.1158E+03, 0.4488E+01, 0.0000E+00},
        {29, 2, 1, 0.9470E+03, 0.2826E+03, 0.1093E+03, 0.6688E+02, 0.3668E+01, 0.9174E-06},
        {29, 2, 0, 0.1106E+04, 0.2502E+03, 0.5938E+01, 0.2402E+02, 0.4576E+01, 0.0000E+00},
        {29, 1, 0, 0.8988E+04, 0.1788E+04, 0.4870E+01, 0.7645E+02, 0.1451E+01, 0.0000E+00},
        {30, 4, 0, 0.9394E+01, 0.1673E+02, 0.1236E+01, 0.1029E+04, 0.4259E+01, 0.3962E-01},
        {30, 3, 2, 0.1730E+02, 0.1818E+02, 0.1017E+05, 0.5288E+01, 0.1736E+02, 0.4667E+00},
        {30, 3, 1, 0.9700E+02, 0.1025E+03, 0.3930E+02, 0.1876E+02, 0.7456E+01, 0.2559E+00},
        {30, 3, 0, 0.1450E+03, 0.6195E+02, 0.4094E+01, 0.1086E+03, 0.4614E+01, 0.0000E+00},
        {30, 2, 1, 0.1037E+04, 0.3486E+03, 0.7784E+02, 0.8298E+02, 0.3391E+01, 0.1125E-07},
        {30, 2, 0, 0.1203E+04, 0.9755E+02, 0.9077E+01, 0.3219E+02, 0.5888E+01, 0.0000E+00},
        {30, 1, 0, 0.9667E+04, 0.8320E+03, 0.2586E+02, 0.4497E+02, 0.2215E+01, 0.0000E+00}};

    // fluorescence parameters from Perkins et al. 1991, derived from the LLNL Evaluated Electron Data Library (EEDL)
    // there is a separate record for each fluorescence transition
    // fluorescence records MUST be sorted on Z and n in the same way as the cross section parameter records
    struct FluorescenceParams
    {
        short Z;       // atomic number
        short n;       // principal quantum number of the shell (always 1 because we only support K shell)
        double omega;  // fluorescence yield (1)
        double E;      // energy of the emitted photon (eV)
    };
    constexpr std::initializer_list<FluorescenceParams> fluorescenceParams = {
        {5, 1, 6.24610e-04, 183.3},     // Ka1,2 -- approximate data from https://www.globalsino.com/EM/page4672.html
        {6, 1, 5.61488e-04, 282.02},    // Ka2
        {6, 1, 1.12060e-03, 282.03},    // Ka1
        {7, 1, 1.09420e-03, 393.35},    // Ka2
        {7, 1, 2.18181e-03, 393.37},    // Ka1
        {8, 1, 1.90768e-03, 523.09},    // Ka2
        {8, 1, 3.80027e-03, 523.13},    // Ka1
        {9, 1, 3.06841e-03, 671.32},    // Ka2
        {9, 1, 6.10743e-03, 671.39},    // Ka1
        {10, 1, 4.64329e-03, 838.1},    // Ka2
        {10, 1, 9.22967e-03, 838.22},   // Ka1
        {11, 1, 6.68996e-03, 1027.58},  // Ka2
        {11, 1, 1.32959e-02, 1027.78},  // Ka1
        {11, 1, 3.80328e-05, 28.06},    // Kb3
        {11, 1, 7.79865e-05, 28.26},    // Kb1
        {12, 1, 9.27327e-03, 1237.95},  // Ka2
        {12, 1, 1.84189e-02, 1238.26},  // Ka1
        {12, 1, 5.79260e-05, 32.91},    // Kb3
        {12, 1, 1.20420e-04, 33.22},    // Kb1
        {13, 1, 1.23699e-02, 1468.7},   // Ka2
        {13, 1, 2.45528e-02, 1469.17},  // Ka1
        {13, 1, 7.55854e-05, 1545.02},  // Kb3
        {13, 1, 1.50039e-04, 1545.03},  // Kb1
        {14, 1, 1.59791e-02, 1719.83},  // Ka2
        {14, 1, 3.17052e-02, 1720.52},  // Ka1
        {14, 1, 2.72402e-04, 1821.95},  // Kb3
        {14, 1, 5.40444e-04, 1821.98},  // Kb1
        {15, 1, 2.01220e-02, 1991.26},  // Ka2
        {15, 1, 3.98749e-02, 1992.22},  // Ka1
        {15, 1, 6.21469e-04, 2122.02},  // Kb3
        {15, 1, 1.23210e-03, 2122.07},  // Kb1
        {16, 1, 2.47822e-02, 2283.17},  // Ka2
        {16, 1, 4.90644e-02, 2284.5},   // Ka1
        {16, 1, 1.15591e-03, 2445.56},  // Kb3
        {16, 1, 2.28902e-03, 2445.65},  // Kb1
        {17, 1, 2.99473e-02, 2595.42},  // Ka2
        {17, 1, 5.92357e-02, 2597.2},   // Ka1
        {17, 1, 1.90852e-03, 2792.48},  // Kb3
        {17, 1, 3.77764e-03, 2792.61},  // Kb1
        {18, 1, 3.55868e-02, 2928.17},  // Ka2
        {18, 1, 7.03336e-02, 2930.51},  // Ka1
        {18, 1, 2.91258e-03, 3162.98},  // Kb3
        {18, 1, 5.75646e-03, 3163.17},  // Kb1
        {19, 1, 4.18599e-02, 3281.65},  // Ka2
        {19, 1, 8.26518e-02, 3284.66},  // Ka1
        {19, 1, 3.99739e-03, 3559.55},  // Kb3
        {19, 1, 7.91008e-03, 3559.83},  // Kb1
        {20, 1, 4.87196e-02, 3655.96},  // Ka2
        {20, 1, 9.61323e-02, 3659.8},   // Ka1
        {20, 1, 5.17776e-03, 3980.99},  // Kb3
        {20, 1, 1.02489e-02, 3981.4},   // Kb1
        {21, 1, 5.63841e-02, 4052.12},  // Ka2
        {21, 1, 1.11130e-01, 4056.94},  // Ka1
        {21, 1, 6.22301e-03, 4426.2},   // Kb3
        {21, 1, 1.23040e-02, 4426.74},  // Kb1
        {22, 1, 6.45982e-02, 4469.38},  // Ka2
        {22, 1, 1.27140e-01, 4475.36},  // Ka1
        {22, 1, 7.32672e-03, 4895.5},   // Kb3
        {22, 1, 1.44750e-02, 4896.18},  // Kb1
        {23, 1, 7.32827e-02, 4907.72},  // Ka2
        {23, 1, 1.44051e-01, 4915.06},  // Ka1
        {23, 1, 8.48408e-03, 5388.95},  // Kb3
        {23, 1, 1.67441e-02, 5389.81},  // Kb1
        {24, 1, 8.25759e-02, 5367.86},  // Ka2
        {24, 1, 1.62090e-01, 5376.79},  // Ka1
        {24, 1, 9.49389e-03, 5906.35},  // Kb3
        {24, 1, 1.87010e-02, 5907.4},   // Kb1
        {25, 1, 9.17726e-02, 5847.93},  // Ka2
        {25, 1, 1.79949e-01, 5858.68},  // Ka1
        {25, 1, 1.09360e-02, 6448.81},  // Kb3
        {25, 1, 2.15229e-02, 6450.12},  // Kb1
        {26, 1, 1.01391e-01, 6349.85},  // Ka2
        {26, 1, 1.98621e-01, 6362.71},  // Ka1
        {26, 1, 1.22111e-02, 7015.36},  // Kb3
        {26, 1, 2.40042e-02, 7016.95},  // Kb1
        {27, 1, 1.11220e-01, 6873.16},  // Ka2
        {27, 1, 2.17470e-01, 6888.41},  // Ka1
        {27, 1, 1.35020e-02, 7606.53},  // Kb3
        {27, 1, 2.65130e-02, 7608.44},  // Kb1
        {28, 1, 1.21060e-01, 7417.82},  // Ka2
        {28, 1, 2.36419e-01, 7435.78},  // Ka1
        {28, 1, 1.48050e-02, 8222.32},  // Kb3
        {28, 1, 2.90299e-02, 8224.59},  // Kb1
        {29, 1, 1.31119e-01, 7984.67},  // Ka2
        {29, 1, 2.55668e-01, 8005.71},  // Ka1
        {29, 1, 1.58899e-02, 8862.7},   // Kb3
        {29, 1, 3.10908e-02, 8865.34},  // Kb1
        {30, 1, 1.40620e-01, 8571.9},   // Ka2
        {30, 1, 2.73809e-01, 8596.4},   // Ka1
        {30, 1, 1.73890e-02, 9528.66},  // Kb3
        {30, 1, 3.39969e-02, 9531.81},  // Kb1
    };

    // wavelength range over which our cross sections may be nonzero
    constexpr Range nonZeroRange(4e-12, 290e-9);  // 300 keV --> 4.3 eV

    // number of wavelengths per dex in high-resolution grid
    constexpr size_t numWavelengthsPerDex = 2500;
}

namespace
{
    // convert photon energy in eV to and from wavelength in m (same conversion in both directions)
    double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    // return thermal velocity for given gas temperature (in K) and particle mass (in amu)
    double vtherm(double T, double amu) { return sqrt(Constants::k() / Constants::amu() * T / amu); }

    // return cross section in m2 for given energy in eV and cross section parameters,
    // without taking into account thermal dispersion
    double crossSection(double E, const CrossSectionParams& p)
    {
        if (E < p.Eth) return 0.;

        double Q = 5.5 + p.l - 0.5 * p.P;
        double y = E / p.E0;
        double ym1 = y - 1.;
        double F = (ym1 * ym1 + p.yw * p.yw) * std::pow(y, -Q) * std::pow(1. + std::sqrt(y / p.ya), -p.P);
        return 1e-22 * p.sigma0 * F;  // from Mb to m2
    }

    // return cross section in m2 for given energy in eV and cross section parameters,
    // approximating thermal dispersion by replacing the steep threshold transition by a sigmoid
    // error function with given parameters (dispersion and maximum value)
    double crossSection(double E, std::pair<double, double> sigmoid, const CrossSectionParams& p)
    {
        double Es, sigmamax;
        std::tie(Es, sigmamax) = sigmoid;
        if (E <= p.Eth - 2. * Es) return 0.;
        if (E >= p.Eth + 2. * Es) return crossSection(E, p);
        return sigmamax * (0.5 + 0.5 * std::erf((E - p.Eth) / Es));
    }
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();
    auto config = find<Configuration>();

    // ---- abundancies ----

    // verify the number of abundancies; if the list is empty, use our default list
    if (_abundancies.empty())
        _abundancies = defaultAbundancies;
    else if (_abundancies.size() != numAtoms)
        throw FATALERROR("The abundancies list must have exactly " + std::to_string(numAtoms) + " values");

    // ---- thermal dispersion ----

    // copy the atom masses (in amu) into a temporary vector
    vector<double> massv = masses;

    // calculate the parameters for the sigmoid function approximating the convolution with a Gaussian
    // at the threshold energy for each cross section record, and and store the result into a temporary vector;
    // the information includes the thermal energy dispersion at the threshold energy and the intrinsic cross section
    // at the threshold energy plus twice this energy dispersion
    vector<std::pair<double, double>> sigmoidv;
    sigmoidv.reserve(crossSectionParams.size());
    for (const auto& params : crossSectionParams)
    {
        double Es = params.Eth * vtherm(temperature(), massv[params.Z - 1]) / Constants::c();
        double sigmamax = crossSection(params.Eth + 2. * Es, params);
        sigmoidv.emplace_back(Es, sigmamax);
    }

    // calculate and store the thermal velocities corresponding to the fluorescence transitions
    _fluovthermv.reserve(fluorescenceParams.size());
    for (const auto& params : fluorescenceParams)
    {
        _fluovthermv.push_back(vtherm(temperature(), massv[params.Z - 1]));
    }

    // ---- wavelength grid ----

    // construct a wavelength grid for sampling cross sections containing a merged set of grid points
    // in the relevant wavelength range (intersection of simulation range and nonzero range):
    //  - a fine grid in log space that provides sufficient resolution for most applications
    //  - all specific wavelengths mentioned in the configuration of the simulation (grids, normalizations, ...)
    //    ensuring that the cross sections are calculated at exactly these wavelengths
    //  - 7 extra wavelength points around the threshold energies for all transitions,
    //    placed at -2, -4/3, -2/3, 0, 2/3, 4/3, 2 times the thermal energy dispersion

    // we first gather all the wavelength points, in arbitrary order, and then sort them
    vector<double> lambdav;
    lambdav.reserve(5 * numWavelengthsPerDex);

    // get the relevant range (intersection of simulation range and nonzero range)
    Range range = config->simulationWavelengthRange();
    range.intersect(nonZeroRange);

    // add a fine grid in log space;
    // use integer multiples as logarithmic grid points so that the grid is stable for changing wavelength ranges
    constexpr double numPerDex = numWavelengthsPerDex;  // converted to double to avoid casting
    int minLambdaSerial = std::floor(numPerDex * log10(range.min()));
    int maxLambdaSerial = std::ceil(numPerDex * log10(range.max()));
    for (int k = minLambdaSerial; k <= maxLambdaSerial; ++k) lambdav.push_back(pow(10., k / numPerDex));

    // add the wavelengths mentioned in the configuration of the simulation
    for (double lambda : config->simulationWavelengths())
        if (range.contains(lambda)) lambdav.push_back(lambda);

    // add wavelength points around the threshold energies for all transitions
    int index = 0;
    for (const auto& params : crossSectionParams)
    {
        double Es = sigmoidv[index++].first;
        for (double delta : {-2., -4. / 3., -2. / 3., 0., 2. / 3., 4. / 3., 2.})
        {
            double lambda = wavelengthToFromEnergy(params.Eth + delta * Es);
            if (range.contains(lambda)) lambdav.push_back(lambda);
        }
    }

    // add the fluorescence emission wavelengths
    for (const auto& params : fluorescenceParams)
    {
        double lambda = wavelengthToFromEnergy(params.E);
        if (range.contains(lambda)) lambdav.push_back(lambda);
    }

    // add the outer wavelengths of our nonzero range so that there are always at least two points in the grid
    lambdav.push_back(nonZeroRange.min());
    lambdav.push_back(nonZeroRange.max());

    // sort the wavelengths and remove duplicates
    NR::unique(lambdav);
    int numLambda = lambdav.size();

    // derive a wavelength grid that will be used for converting a wavelength to an index in the above array;
    // the grid points are shifted to the left of the actual sample points to approximate rounding
    _lambdav.resize(numLambda);
    _lambdav[0] = lambdav[0];
    for (int ell = 1; ell != numLambda; ++ell)
    {
        _lambdav[ell] = sqrt(lambdav[ell] * lambdav[ell - 1]);
    }

    // ---- extinction ----

    // calculate the extinction cross section at every wavelength; to guarantee that the cross section is zero
    // for wavelengths outside our range, leave the values for the outer wavelength points at zero
    _sigmaextv.resize(numLambda);
    for (int ell = 1; ell < numLambda - 1; ++ell)
    {
        double E = wavelengthToFromEnergy(lambdav[ell]);
        double sigma = 0.;
        int index = 0;
        for (const auto& params : crossSectionParams)
        {
            sigma += crossSection(E, sigmoidv[index++], params) * _abundancies[params.Z - 1];
        }
        _sigmaextv[ell] = sigma;
    }

    // ---- scattering ----

    // calculate and store the fluorescence emission wavelengths
    _fluolambdav.reserve(fluorescenceParams.size());
    for (const auto& params : fluorescenceParams) _fluolambdav.push_back(wavelengthToFromEnergy(params.E));

    // make room for the scattering cross section and the cumulative fluorescence probabilities at every wavelength
    _sigmascav.resize(numLambda);
    _fluocumprobvv.resize(numLambda, 0);

    // provide temporary array for the non-normalized fluorescence contributions (at the current wavelength)
    Array flucontribv(fluorescenceParams.size());

    // calculate the above for every wavelength; as before, leave the values for the outer wavelength points at zero
    for (int ell = 1; ell < numLambda - 1; ++ell)
    {
        double E = wavelengthToFromEnergy(lambdav[ell]);
        double sigma = 0.;

        // interate over both cross section and fluorescence parameter sets in sync
        const auto* flp = fluorescenceParams.begin();
        int index = 0;
        for (const auto& csp : crossSectionParams)
        {
            auto sigmoid = sigmoidv[index++];

            // process all fluorescence parameter sets matching this cross section set
            while (flp != fluorescenceParams.end() && flp->Z == csp.Z && flp->n == csp.n)
            {
                double contribution = crossSection(E, sigmoid, csp) * _abundancies[csp.Z - 1] * flp->omega;
                sigma += contribution;
                flucontribv[flp - fluorescenceParams.begin()] = contribution;
                flp++;
            }
        }

        // store the cross section and determine the normalized cumulative probability distribution
        _sigmascav[ell] = sigma;
        NR::cdf(_fluocumprobvv[ell], flucontribv);
    }
}

////////////////////////////////////////////////////////////////////

int XRayAtomicGasMix::indexForLambda(double lambda) const
{
    return NR::locateClip(_lambdav, lambda);
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType XRayAtomicGasMix::materialType() const
{
    return MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool XRayAtomicGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayAtomicGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity()};
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::sectionAbs(double lambda) const
{
    int index = indexForLambda(lambda);
    return _sigmaextv[index] - _sigmascav[index];
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::sectionSca(double lambda) const
{
    return _sigmascav[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::sectionExt(double lambda) const
{
    return _sigmaextv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionAbs(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionSca(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionExt(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda) const
{
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;
        scatinfo->species = NR::locateClip(_fluocumprobvv[indexForLambda(lambda)], random()->uniform());
        scatinfo->velocity = _fluovthermv[scatinfo->species] * random()->maxwell();
    }
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::peeloffScattering(double& I, double& /*Q*/, double& /*U*/, double& /*V*/, double& lambda,
                                         double w, Direction bfkobs, Direction /*bfky*/, const MaterialState* /*state*/,
                                         const PhotonPacket* pp) const
{
    // draw a random fluorescence channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda);

    // isotropic scattering, so the contribution is trivially 1 (multiplied by the weight for this component)
    I += w;

    // for a random fraction of the events governed by the weight for this component,
    // update the photon packet wavelength to the wavelength of this fluorescence transition,
    // Doppler-shifted out of the atom velocity frame
    if (random()->uniform() <= w)
    {
        lambda = PhotonPacket::shiftedEmissionWavelength(_fluolambdav[scatinfo->species], bfkobs, scatinfo->velocity);
    }
}

////////////////////////////////////////////////////////////////////

void XRayAtomicGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random fluorescence channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = pp->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda);

    // draw a random, isotropic outgoing direction
    Direction bfknew = random()->direction();

    // update the photon packet wavelength to the wavelength of this fluorescence transition,
    // Doppler-shifted out of the atom velocity frame
    lambda = PhotonPacket::shiftedEmissionWavelength(_fluolambdav[scatinfo->species], bfknew, scatinfo->velocity);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

double XRayAtomicGasMix::indicativeTemperature(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    return temperature();
}

////////////////////////////////////////////////////////////////////
