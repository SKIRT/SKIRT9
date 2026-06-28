/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NEBULARLINEEMISSION_HPP
#define NEBULARLINEEMISSION_HPP

#include <algorithm>
#include <cmath>

//////////////////////////////////////////////////////////////////////

/** The NebularLineEmission namespace provides static data tables and functions for computing
    nebular emission from photoionized gas:

    - H recombination lines (Case B): Storey & Hummer (1995) P_B tables with Hui & Gnedin (1997)
      alpha_B. Lines: Lyman-alpha, Balmer (alpha through epsilon), Paschen (alpha, beta),
      Brackett alpha.
    - Metal forbidden lines: pre-tabulated collisional excitation rate coefficients q_col(T, n_e)
      from CHIANTI statistical equilibrium. Lines: [NII] 6548/6583, [OI] 6300/6364,
      [OII] 3727/3729, [OIII] 4363/4959/5007, [SII] 6716/6731.

    Implementation follows McClymont, Smith & Tacchella (2025). */

namespace NebularLineEmission
{
    // ============== Line definitions ==============

    // Line indices for the emission line array
    enum LineIndex {
        // H recombination lines (Storey & Hummer 1995 Case B)
        Lya = 0,   // Lyman-alpha   1215.67 A
        Ha,        // Balmer-alpha   6562.80 A
        Hb,        // Balmer-beta    4861.33 A
        Hg,        // Balmer-gamma   4340.46 A
        Hd,        // Balmer-delta   4101.73 A
        HeBalmer,  // Balmer-epsilon 3970.07 A
        Paa,       // Paschen-alpha  18751 A
        Pab,       // Paschen-beta   12818 A
        Bra,       // Brackett-alpha 40512 A

        // Metal forbidden lines (CHIANTI collisional excitation)
        NII6548,   // [NII] 6548
        NII6583,   // [NII] 6583
        OI6300,    // [OI] 6300
        OI6364,    // [OI] 6364
        OII3729,   // [OII] 3729  (2D 5/2)
        OII3726,   // [OII] 3726  (2D 3/2)
        OIII4363,  // [OIII] 4363
        OIII4959,  // [OIII] 4959
        OIII5007,  // [OIII] 5007
        SII6716,   // [SII] 6716
        SII6731,   // [SII] 6731
    };
    constexpr int numLines = 20;

    // Rest-frame wavelengths [m]
    constexpr double lineWavelengths[numLines] = {
        1215.670e-10,  // Lya
        6562.800e-10,  // Ha
        4861.330e-10,  // Hb
        4340.460e-10,  // Hg
        4101.730e-10,  // Hd
        3970.070e-10,  // He
        1.87510e-6,    // Paa
        1.28180e-6,    // Pab
        4.05120e-6,    // Bra
        6548.050e-10,  // [NII] 6548
        6583.460e-10,  // [NII] 6583
        6300.304e-10,  // [OI] 6300
        6363.776e-10,  // [OI] 6364
        3728.815e-10,  // [OII] 3729
        3726.032e-10,  // [OII] 3726
        4363.210e-10,  // [OIII] 4363
        4958.911e-10,  // [OIII] 4959
        5006.843e-10,  // [OIII] 5007
        6716.440e-10,  // [SII] 6716
        6730.810e-10,  // [SII] 6731
    };

    // Particle mass for Doppler broadening [kg]: proton mass for H, approximate for metals.
    constexpr double Mproton = 1.67262192e-27;
    constexpr double lineMasses[numLines] = {
        Mproton,         // Lya (H)
        Mproton,         // Ha
        Mproton,         // Hb
        Mproton,         // Hg
        Mproton,         // Hd
        Mproton,         // He
        Mproton,         // Paa
        Mproton,         // Pab
        Mproton,         // Bra
        14.0 * Mproton,  // NII
        14.0 * Mproton,  // NII
        16.0 * Mproton,  // OI
        16.0 * Mproton,  // OI
        16.0 * Mproton,  // OII
        16.0 * Mproton,  // OII
        16.0 * Mproton,  // OIII
        16.0 * Mproton,  // OIII
        16.0 * Mproton,  // OIII
        32.0 * Mproton,  // SII
        32.0 * Mproton,  // SII
    };

    // Ion stage index in PhotoIonizationSolver::ionFracs[] for carrier ion of each metal line
    // N: stages 12-19 (NI=12, NII=13, ...)
    // O: stages 20-28 (OI=20, OII=21, OIII=22, ...)
    // S: stages 62-66 (SI=62, SII=63, ...)
    constexpr int lineCarrierIonIndex[numLines] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1,  // H lines: not from ion fractions
        13,                                  // [NII]  -> NII  = ion index 13
        13,                                  // [NII]  -> NII
        20,                                  // [OI]   -> OI   = ion index 20
        20,                                  // [OI]   -> OI
        21,                                  // [OII]  -> OII  = ion index 21
        21,                                  // [OII]  -> OII
        22,                                  // [OIII] -> OIII = ion index 22
        22,                                  // [OIII] -> OIII
        22,                                  // [OIII] -> OIII
        63,                                  // [SII]  -> SII  = ion index 63
        63,                                  // [SII]  -> SII
    };

    // Element index for each metal line (used to get abundance)
    // 0=C, 1=N, 2=O, 3=Ne, 4=Mg, 5=Si, 6=S, 7=Fe
    constexpr int lineElementIndex[numLines] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1,  // H lines
        1,  1,                               // NII -> N
        2,  2,                               // OI -> O
        2,  2,                               // OII -> O
        2,  2,  2,                           // OIII -> O
        6,  6,                               // SII -> S
    };

    // ============== H recombination line data (Storey & Hummer 1995) ==============

    // Temperature grid: 10 points from 500 K to 30000 K (non-uniform in log space)
    constexpr int pbNumT = 10;
    constexpr int pbNumNe = 7;

    // log10(T) values computed at compile time from the actual temperatures
    // T = {500, 1000, 3000, 5000, 7500, 10000, 12500, 15000, 20000, 30000}
    // n_e = {1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8}
    inline double logT_SH95(int i)
    {
        static const double v[pbNumT] = {2.69897000, 3.0,        3.47712125, 3.69897000, 3.87506126,
                                         4.0,        4.09691001, 4.17609126, 4.30103000, 4.47712125};
        return v[i];
    }

    // P_B tables: probability of emitting specific line photon per recombination
    // Indexed as [temperature][ne], 10x7 per line

    // Lyman-alpha
    constexpr double P_Lya[pbNumT][pbNumNe] = {
        {0.80793637, 0.90215197, 0.98157401, 0.99705340, 0.99933853, 0.99934845, 0.99967896},
        {0.78439303, 0.87743841, 0.97423464, 0.99680114, 0.99960243, 0.99967816, 0.99976710},
        {0.74159266, 0.82809599, 0.95734434, 0.99359891, 0.99892162, 0.99927343, 0.99979471},
        {0.71946081, 0.80133228, 0.94581237, 0.99246182, 0.99848958, 0.99940329, 0.99949246},
        {0.70117195, 0.77884530, 0.93441344, 0.99075058, 0.99831301, 0.99960223, 0.99940949},
        {0.68808416, 0.76220739, 0.92503048, 0.98912680, 0.99828523, 0.99940740, 0.99944339},
        {0.67755810, 0.74868826, 0.91762730, 0.98784252, 0.99835238, 0.99903921, 0.99977064},
        {0.66862468, 0.73785038, 0.91033402, 0.98701405, 0.99817046, 0.99892625, 0.99953078},
        {0.65534642, 0.72080463, 0.89823705, 0.98497224, 0.99798854, 0.99881330, 0.99929092},
        {0.63566376, 0.69600394, 0.87952243, 0.98121463, 0.99716593, 0.99896144, 0.99905958}};

    // Balmer-alpha (Halpha 6563 A)
    constexpr double P_Ha[pbNumT][pbNumNe] = {
        {0.60176416, 0.59203990, 0.57777268, 0.55757870, 0.53218274, 0.50309773, 0.47538990},
        {0.57146208, 0.56548156, 0.55537185, 0.54126081, 0.52171947, 0.49779442, 0.47325611},
        {0.51602376, 0.51308374, 0.50870615, 0.50163462, 0.49109775, 0.47695493, 0.46057657},
        {0.48842989, 0.48670009, 0.48384243, 0.47916786, 0.47205133, 0.46211067, 0.45011817},
        {0.46673159, 0.46571739, 0.46369825, 0.46057413, 0.45555148, 0.44857577, 0.43958049},
        {0.45180723, 0.45102825, 0.44947568, 0.44724222, 0.44349888, 0.43815604, 0.43129878},
        {0.44059820, 0.43967487, 0.43880459, 0.43692031, 0.43424518, 0.42968461, 0.42434709},
        {0.43139784, 0.43080340, 0.43015528, 0.42868403, 0.42629594, 0.42286913, 0.41840805},
        {0.41795682, 0.41743320, 0.41691031, 0.41586673, 0.41447832, 0.41195330, 0.40867726},
        {0.39957948, 0.39949886, 0.39933772, 0.39872322, 0.39820159, 0.39687074, 0.39529763}};

    // Balmer-beta (Hbeta 4861 A)
    constexpr double P_Hb[pbNumT][pbNumNe] = {
        {0.11197758, 0.11382091, 0.11630471, 0.11882216, 0.12108799, 0.12288014, 0.12720456},
        {0.11559987, 0.11672083, 0.11811675, 0.11985224, 0.12129050, 0.12256368, 0.12615841},
        {0.11908521, 0.11942732, 0.11988399, 0.12040081, 0.12086866, 0.12131817, 0.12352423},
        {0.11897682, 0.11914759, 0.11943147, 0.11964361, 0.11981875, 0.11999444, 0.12172966},
        {0.11803310, 0.11815162, 0.11823874, 0.11840903, 0.11841351, 0.11845008, 0.11983194},
        {0.11688913, 0.11694270, 0.11695517, 0.11706977, 0.11706124, 0.11708462, 0.11828955},
        {0.11571266, 0.11566493, 0.11578498, 0.11586298, 0.11583936, 0.11573996, 0.11684562},
        {0.11460234, 0.11463318, 0.11469477, 0.11472462, 0.11466497, 0.11459548, 0.11563639},
        {0.11271958, 0.11270919, 0.11273303, 0.11274643, 0.11268794, 0.11258888, 0.11350468},
        {0.10960571, 0.10960828, 0.10961341, 0.10961006, 0.10954052, 0.10948666, 0.11028151}};

    // Balmer-gamma (Hgamma 4340 A)
    constexpr double P_Hg[pbNumT][pbNumNe] = {
        {0.04201613, 0.04300207, 0.04442060, 0.04620360, 0.04833380, 0.05180051, 0.05900800},
        {0.04445549, 0.04508442, 0.04597248, 0.04723395, 0.04878288, 0.05153151, 0.05778401},
        {0.04780603, 0.04803281, 0.04839422, 0.04889559, 0.04959578, 0.05115426, 0.05535776},
        {0.04869447, 0.04880649, 0.04905029, 0.04932562, 0.04972818, 0.05084741, 0.05406425},
        {0.04894908, 0.04902770, 0.04914434, 0.04931948, 0.04956781, 0.05037671, 0.05290131},
        {0.04887845, 0.04892308, 0.04898677, 0.04910065, 0.04927799, 0.04990298, 0.05202193},
        {0.04867211, 0.04867532, 0.04876216, 0.04882675, 0.04896197, 0.04944524, 0.05131098},
        {0.04838950, 0.04842262, 0.04846499, 0.04852574, 0.04859925, 0.04902806, 0.05067135},
        {0.04785069, 0.04784778, 0.04787542, 0.04790009, 0.04793400, 0.04825843, 0.04963792},
        {0.04679314, 0.04678370, 0.04678686, 0.04679788, 0.04678062, 0.04701408, 0.04811984}};

    // Balmer-delta (Hdelta 4102 A)
    constexpr double P_Hd[pbNumT][pbNumNe] = {
        {0.02109056, 0.02162223, 0.02242694, 0.02352853, 0.02523019, 0.02929151, 0.03647097},
        {0.02251388, 0.02286241, 0.02337582, 0.02416699, 0.02544500, 0.02872042, 0.03513255},
        {0.02461022, 0.02474142, 0.02495551, 0.02529094, 0.02594002, 0.02793451, 0.03223331},
        {0.02524504, 0.02531937, 0.02544631, 0.02564597, 0.02606706, 0.02753598, 0.03072347},
        {0.02549994, 0.02554921, 0.02562219, 0.02575119, 0.02603134, 0.02715664, 0.02951235},
        {0.02553874, 0.02557294, 0.02561725, 0.02569328, 0.02590854, 0.02682556, 0.02866888},
        {0.02548653, 0.02548937, 0.02554250, 0.02559362, 0.02575902, 0.02652537, 0.02803061},
        {0.02537701, 0.02539690, 0.02541418, 0.02545984, 0.02558400, 0.02626338, 0.02750274},
        {0.02513928, 0.02513614, 0.02514743, 0.02516997, 0.02525799, 0.02580427, 0.02670334},
        {0.02461963, 0.02461466, 0.02462555, 0.02460818, 0.02465908, 0.02507381, 0.02563413}};

    // Balmer-epsilon (Hepsilon 3970 A)
    constexpr double P_He[pbNumT][pbNumNe] = {
        {0.01236677, 0.01269669, 0.01320054, 0.01396256, 0.01560329, 0.02011311, 0.02620399},
        {0.01324123, 0.01345937, 0.01379294, 0.01434085, 0.01559826, 0.01928819, 0.02486383},
        {0.01455964, 0.01464171, 0.01478131, 0.01503404, 0.01572268, 0.01797897, 0.02178548},
        {0.01497012, 0.01501864, 0.01509961, 0.01525616, 0.01573583, 0.01739606, 0.02017733},
        {0.01514584, 0.01517759, 0.01522259, 0.01532706, 0.01567484, 0.01693436, 0.01892051},
        {0.01518076, 0.01520176, 0.01522822, 0.01530019, 0.01557785, 0.01659052, 0.01806633},
        {0.01516266, 0.01516008, 0.01519246, 0.01523581, 0.01547197, 0.01632220, 0.01743615},
        {0.01509442, 0.01510795, 0.01512412, 0.01516718, 0.01534991, 0.01609385, 0.01694092},
        {0.01495760, 0.01496112, 0.01496463, 0.01498559, 0.01513864, 0.01573279, 0.01619841},
        {0.01465647, 0.01465553, 0.01465163, 0.01465137, 0.01474681, 0.01520321, 0.01525805}};

    // Paschen-alpha (18751 A)
    constexpr double P_Paa[pbNumT][pbNumNe] = {
        {0.32614959, 0.31160767, 0.29095068, 0.26378240, 0.23187496, 0.19849730, 0.16865559},
        {0.28573104, 0.27671586, 0.26290888, 0.24406397, 0.21991833, 0.19281748, 0.16664950},
        {0.21805802, 0.21450247, 0.20864657, 0.19975209, 0.18757159, 0.17207439, 0.15505924},
        {0.18842673, 0.18620048, 0.18256720, 0.17697689, 0.16891281, 0.15828873, 0.14579711},
        {0.16674031, 0.16533443, 0.16293730, 0.15925275, 0.15371473, 0.14629550, 0.13705592},
        {0.15267112, 0.15161900, 0.14988586, 0.14720832, 0.14315584, 0.13757837, 0.13039950},
        {0.14252317, 0.14166471, 0.14047895, 0.13835956, 0.13530560, 0.13082510, 0.12513531},
        {0.13471866, 0.13413161, 0.13311337, 0.13149314, 0.12898747, 0.12542193, 0.12069175},
        {0.12351353, 0.12309690, 0.12241688, 0.12132329, 0.11967100, 0.11724347, 0.11385228},
        {0.10969035, 0.10947782, 0.10914819, 0.10859574, 0.10778889, 0.10657109, 0.10464593}};

    // Paschen-beta (12818 A)
    constexpr double P_Pab[pbNumT][pbNumNe] = {
        {0.07348988, 0.07341128, 0.07269308, 0.07070782, 0.06718798, 0.06260883, 0.05946774},
        {0.07078253, 0.07056248, 0.06978101, 0.06816808, 0.06522349, 0.06126369, 0.05828991},
        {0.06277961, 0.06253897, 0.06202437, 0.06098018, 0.05916578, 0.05653085, 0.05421758},
        {0.05773433, 0.05755307, 0.05717604, 0.05642753, 0.05512529, 0.05315164, 0.05131155},
        {0.05341396, 0.05329313, 0.05301332, 0.05249053, 0.05153516, 0.05006967, 0.04860267},
        {0.05031083, 0.05020266, 0.04998706, 0.04959691, 0.04886864, 0.04771556, 0.04653578},
        {0.04793166, 0.04783464, 0.04770007, 0.04736567, 0.04679777, 0.04583497, 0.04489711},
        {0.04599371, 0.04593356, 0.04581345, 0.04557400, 0.04509822, 0.04428626, 0.04348421},
        {0.04311101, 0.04305376, 0.04296500, 0.04279237, 0.04247102, 0.04187832, 0.04129479},
        {0.03928830, 0.03926736, 0.03921900, 0.03913283, 0.03896478, 0.03862739, 0.03828750}};

    // Brackett-alpha (40512 A)
    constexpr double P_Bra[pbNumT][pbNumNe] = {
        {0.21885188, 0.20451948, 0.18492120, 0.16070099, 0.13417765, 0.10823662, 0.08695037},
        {0.18244521, 0.17373440, 0.16114577, 0.14461364, 0.12492214, 0.10407263, 0.08532126},
        {0.12619750, 0.12290608, 0.11783853, 0.11051795, 0.10092464, 0.08938717, 0.07733933},
        {0.10361101, 0.10162800, 0.09858645, 0.09407032, 0.08785312, 0.08005556, 0.07131783},
        {0.08789567, 0.08670423, 0.08470428, 0.08177624, 0.07762419, 0.07217453, 0.06584045},
        {0.07813115, 0.07725277, 0.07583222, 0.07370886, 0.07067981, 0.06661949, 0.06171954},
        {0.07127436, 0.07057650, 0.06956207, 0.06791772, 0.06561596, 0.06240387, 0.05852338},
        {0.06610734, 0.06560520, 0.06477995, 0.06351096, 0.06166134, 0.05908966, 0.05589889},
        {0.05886696, 0.05851187, 0.05795768, 0.05712213, 0.05588395, 0.05410592, 0.05185238},
        {0.05017424, 0.04999957, 0.04973270, 0.04930727, 0.04868878, 0.04776841, 0.04652094}};

    /** Bilinear interpolation of P_B(T, n_e) from Storey & Hummer (1995) tables.
        Returns probability of emitting the specific line photon per Case B recombination.
        T in K, ne in cm^-3. */
    inline double interpolate_PB(const double (&Ptab)[pbNumT][pbNumNe], double T, double ne)
    {
        // Temperature interpolation indices (log space, non-uniform grid)
        int iT = 8;
        double fracT = 1.0;
        if (T <= 500.)
        {
            iT = 0;
            fracT = 0.;
        }
        else if (T < 30000.)
        {
            double logT = std::log10(T);
            iT = 8;
            while (logT < logT_SH95(iT)) iT--;
            fracT = (logT - logT_SH95(iT)) / (logT_SH95(iT + 1) - logT_SH95(iT));
        }

        // Density interpolation indices (log space, uniform 1 dex spacing from 10^2 to 10^8)
        int iN = 0;
        double fracN = 0.;
        if (ne >= 1e8)
        {
            iN = 5;
            fracN = 1.;
        }
        else if (ne > 1e2)
        {
            double d = std::log10(ne) - 2.0;
            double f = std::floor(d);
            fracN = d - f;
            iN = static_cast<int>(f);
        }

        // Bilinear interpolation
        double fLT = 1. - fracT, fLN = 1. - fracN;
        int jT = iT + 1, jN = iN + 1;
        return fLT * (Ptab[iT][iN] * fLN + Ptab[iT][jN] * fracN) + fracT * (Ptab[jT][iN] * fLN + Ptab[jT][jN] * fracN);
    }

    /** Case B recombination coefficient for HII [cm^3/s]. Hui & Gnedin (1997). */
    inline double alphaBHII(double T)
    {
        double lambda = 315614.0 / T;  // 2 * T_i / T where T_i = 157807 K
        return 2.753e-14 * std::pow(lambda, 1.5) / std::pow(1.0 + std::pow(lambda / 2.74, 0.407), 2.242);
    }

    /** Effective Case B recombination coefficient for a given H line [cm^3/s].
        alpha_eff = P_B(T, ne) * alpha_B(T). */
    inline double alphaEffB(int lineIdx, double T, double ne)
    {
        switch (lineIdx)
        {
            case LineIndex::Lya: return interpolate_PB(P_Lya, T, ne) * alphaBHII(T);
            case LineIndex::Ha: return interpolate_PB(P_Ha, T, ne) * alphaBHII(T);
            case LineIndex::Hb: return interpolate_PB(P_Hb, T, ne) * alphaBHII(T);
            case LineIndex::Hg: return interpolate_PB(P_Hg, T, ne) * alphaBHII(T);
            case LineIndex::Hd: return interpolate_PB(P_Hd, T, ne) * alphaBHII(T);
            case LineIndex::HeBalmer: return interpolate_PB(P_He, T, ne) * alphaBHII(T);
            case LineIndex::Paa: return interpolate_PB(P_Paa, T, ne) * alphaBHII(T);
            case LineIndex::Pab: return interpolate_PB(P_Pab, T, ne) * alphaBHII(T);
            case LineIndex::Bra: return interpolate_PB(P_Bra, T, ne) * alphaBHII(T);
            default: return 0.;
        }
    }

    /** Case B line-emission branching probability P_B(T, ne) per Case B recombination,
        for a given hydrogen line. Dimensionless. Returns 0 for non-H lines. */
    inline double linePBonly(int lineIdx, double T, double ne)
    {
        switch (lineIdx)
        {
            case LineIndex::Lya: return interpolate_PB(P_Lya, T, ne);
            case LineIndex::Ha: return interpolate_PB(P_Ha, T, ne);
            case LineIndex::Hb: return interpolate_PB(P_Hb, T, ne);
            case LineIndex::Hg: return interpolate_PB(P_Hg, T, ne);
            case LineIndex::Hd: return interpolate_PB(P_Hd, T, ne);
            case LineIndex::HeBalmer: return interpolate_PB(P_He, T, ne);
            case LineIndex::Paa: return interpolate_PB(P_Paa, T, ne);
            case LineIndex::Pab: return interpolate_PB(P_Pab, T, ne);
            case LineIndex::Bra: return interpolate_PB(P_Bra, T, ne);
            default: return 0.;
        }
    }

    // ============== Metal forbidden line data (CHIANTI collisional excitation) ==============

    // Grid for q_col tables: 21 temperatures x 19 densities
    constexpr int qcNumT = 21;
    constexpr int qcNumNe = 19;
    constexpr int qcOrderPlusOne = 4;  // Lagrange interpolation order + 1
    constexpr double qcMinLogT = 2.0;
    constexpr double qcMaxLogT = 7.0;
    constexpr double qcMinLogN = 1.0;
    constexpr double qcMaxLogN = 10.0;
    constexpr double qcInvDexT = 4.0;  // 1/0.25 dex
    constexpr double qcInvDexN = 2.0;  // 1/0.5 dex

    // Temperature grid: logT = {2.0, 2.25, 2.5, ..., 7.0}
    inline double qc_logT(int i)
    {
        return 2.0 + 0.25 * i;
    }
    // Density grid: logn = {1.0, 1.5, 2.0, ..., 10.0}
    inline double qc_logne(int i)
    {
        return 1.0 + 0.5 * i;
    }

    /** 4-point Lagrange interpolation. */
    inline double lagrange4(const double x[4], const double y[4], double xval)
    {
        double yval = 0.;
        for (int i = 0; i < 4; ++i)
        {
            double l = 1.;
            for (int j = 0; j < 4; ++j)
            {
                if (i != j) l *= (xval - x[j]) / (x[i] - x[j]);
            }
            yval += y[i] * l;
        }
        return yval;
    }

    // Number of metal lines
    constexpr int numMetalLines = 11;  // NII x2, OI x2, OII x2, OIII x3, SII x2

    // q_col tables: log10(q_col) in cm^3/s, indexed as [ne_index][T_index]
    // Each line has its own qcNumNe x qcNumT table

    // Forward declaration: the actual tables are defined in NebularLineEmission.cpp
    // to avoid a massive header. Only access functions are declared here.

    /** Returns log10(q_col(T, ne)) for the given metal line index (LineIndex::NII6548 through LineIndex::SII6731).
        T in K, ne in cm^-3. Returns the collisional excitation rate coefficient [cm^3/s]. */
    double metalLineQcol(int lineIdx, double T, double ne);

    /** Returns log10(q'_col(T, ne)) for OII-3726 (the doublet companion of OII-3729).
        For other lines, returns 0. */
    double metalLineQpCol(int lineIdx, double T, double ne);

    // ============== Emission computation ==============

    // Physical constants in CGS
    constexpr double h_cgs = 6.62607015e-27;  // Planck constant [erg s]
    constexpr double c_cgs = 2.99792458e10;   // speed of light [cm/s]
    constexpr double kB_cgs = 1.380649e-16;   // Boltzmann constant [erg/K]

    // SI constants
    constexpr double h_SI = 6.62607015e-34;  // [J s]
    constexpr double c_SI = 2.99792458e8;    // [m/s]

    /** Compute H recombination line luminosity [W] in the photon-conserving form
        L = P_B(T, ne) * Gamma_HI * n_HI * V * h*nu_line. Assumes photoionization balance:
        each H-I photoionization eventually produces one Case B recombination, of which a
        fraction P_B emits in the requested line.
        lineIdx: LineIndex::Lya through LineIndex::Bra
        T: temperature [K]
        ne: electron density [cm^-3]
        gammaHI: H I photoionization rate per H atom [s^-1]
        nHI: neutral hydrogen number density [cm^-3]
        V_cm3: cell volume [cm^3]
        Returns luminosity in [W] (SI). */
    inline double hydrogenLineLuminosity(int lineIdx, double T, double ne, double gammaHI, double nHI, double V_cm3)
    {
        double PB = linePBonly(lineIdx, T, ne);
        double lambda_m = lineWavelengths[lineIdx];
        double hnu = h_SI * c_SI / lambda_m;  // photon energy [J]
        // L = h*nu * P_B * Gamma_HI * n_HI * V [W]
        return hnu * PB * gammaHI * nHI * V_cm3;
    }

    /** Compute metal forbidden line luminosity [W] for a given cell.
        lineIdx: LineIndex::NII6548 through LineIndex::SII6731
        T: temperature [K]
        ne: electron density [cm^-3]
        nIon: carrier ion number density [cm^-3]
        V_cm3: cell volume [cm^3]
        Returns luminosity in [W] (SI). */
    inline double metalLineLuminosity(int lineIdx, double T, double ne, double nIon, double V_cm3)
    {
        double qcol = metalLineQcol(lineIdx, T, ne);
        double lambda_m = lineWavelengths[lineIdx];
        double hnu = h_SI * c_SI / lambda_m;
        // L = h*nu * q_col * ne * nIon * V  [total luminosity in W]
        return hnu * qcol * ne * nIon * V_cm3;
    }

}  // namespace NebularLineEmission

//////////////////////////////////////////////////////////////////////

#endif
