/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SkirtUnitDef.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

SkirtUnitDef::SkirtUnitDef()
{
    // get relevant constants in local variables
    double k = Constants::k();
    double pc = Constants::pc();
    double AU = Constants::AU();
    double Msun = Constants::Msun();
    double Lsun = Constants::Lsun();
    double year = Constants::year();

    // *** add units for each physical quantity

    // length
    addUnit("length", "m", 1.);
    addUnit("length", "cm", 1e-2);
    addUnit("length", "mm", 1e-3);
    addUnit("length", "km", 1e3);
    addUnit("length", "AU", AU);
    addUnit("length", "pc", pc);
    addUnit("length", "kpc", 1e3 * pc);
    addUnit("length", "Mpc", 1e6 * pc);

    // distance
    addUnit("distance", "m", 1.);
    addUnit("distance", "cm", 1e-2);
    addUnit("distance", "mm", 1e-3);
    addUnit("distance", "km", 1e3);
    addUnit("distance", "AU", AU);
    addUnit("distance", "pc", pc);
    addUnit("distance", "kpc", 1e3 * pc);
    addUnit("distance", "Mpc", 1e6 * pc);

    // wavelength
    addUnit("wavelength", "m", 1.);
    addUnit("wavelength", "cm", 1e-2);
    addUnit("wavelength", "mm", 1e-3);
    addUnit("wavelength", "micron", 1e-6);
    addUnit("wavelength", "nm", 1e-9);
    addUnit("wavelength", "Angstrom", 1e-10);

    // grainsize
    addUnit("grainsize", "m", 1.);
    addUnit("grainsize", "cm", 1e-2);
    addUnit("grainsize", "mm", 1e-3);
    addUnit("grainsize", "micron", 1e-6);
    addUnit("grainsize", "nm", 1e-9);
    addUnit("grainsize", "Angstrom", 1e-10);

    // cross section
    addUnit("section", "m2", 1.);
    addUnit("section", "cm2", 1e-4);
    addUnit("section", "mm2", 1e-6);

    // volume
    addUnit("volume", "m3", 1.);
    addUnit("volume", "cm3", 1e-6);
    addUnit("volume", "mm3", 1e-9);
    addUnit("volume", "AU3", pow(AU, 3));
    addUnit("volume", "pc3", pow(pc, 3));

    // velocity
    addUnit("velocity", "m/s", 1.);
    addUnit("velocity", "cm/s", 1e-2);
    addUnit("velocity", "mm/s", 1e-3);
    addUnit("velocity", "km/s", 1e3);
    addUnit("velocity", "km/h", 1000. / 3600.);

    // acceleration
    addUnit("acceleration", "m/s2", 1.);
    addUnit("acceleration", "cm/s2", 1e-2);
    addUnit("acceleration", "mm/s2", 1e-3);
    addUnit("acceleration", "km/s2", 1e3);

    // mass
    addUnit("mass", "kg", 1.);
    addUnit("mass", "g", 1e-3);
    addUnit("mass", "Msun", Msun);

    // bulk mass
    addUnit("bulkmass", "kg", 1.);
    addUnit("bulkmass", "g", 1e-3);

    // bulk mass density
    addUnit("bulkmassdensity", "kg/m3", 1.);
    addUnit("bulkmassdensity", "g/cm3", 1e3);

    // mass surface density
    addUnit("masssurfacedensity", "kg/m2", 1.);
    addUnit("masssurfacedensity", "g/cm2", 10.);
    addUnit("masssurfacedensity", "Msun/AU2", Msun / pow(AU, 2));
    addUnit("masssurfacedensity", "Msun/pc2", Msun / pow(pc, 2));

    // mass volume density
    addUnit("massvolumedensity", "kg/m3", 1.);
    addUnit("massvolumedensity", "g/cm3", 1e3);
    addUnit("massvolumedensity", "Msun/AU3", Msun / pow(AU, 3));
    addUnit("massvolumedensity", "Msun/pc3", Msun / pow(pc, 3));

    // mass rate
    addUnit("massrate", "kg/s", 1.);
    addUnit("massrate", "g/s", 1e-3);
    addUnit("massrate", "Msun/yr", Msun / year);

    // number surface density
    addUnit("numbersurfacedensity", "1/m2", 1.);
    addUnit("numbersurfacedensity", "1/cm2", 1e4);
    addUnit("numbersurfacedensity", "1/AU2", 1 / pow(AU, 2));
    addUnit("numbersurfacedensity", "1/pc2", 1 / pow(pc, 2));

    // number volume density
    addUnit("numbervolumedensity", "1/m3", 1.);
    addUnit("numbervolumedensity", "1/cm3", 1e6);
    addUnit("numbervolumedensity", "1/AU3", 1 / pow(AU, 3));
    addUnit("numbervolumedensity", "1/pc3", 1 / pow(pc, 3));

    // mass coefficient
    addUnit("masscoefficient", "m2/kg", 1.);
    addUnit("masscoefficient", "cm2/g", 0.1);

    // time
    addUnit("time", "s", 1.);
    addUnit("time", "yr", year);
    addUnit("time", "Myr", 1e6 * year);
    addUnit("time", "Gyr", 1e9 * year);

    // temperature
    addUnit("temperature", "K", 1.);

    // energy
    addUnit("energy", "J", 1.);
    addUnit("energy", "erg", 1e-7);

    // magnetic field
    addUnit("magneticfield", "T", 1.);
    addUnit("magneticfield", "mT", 1e-3);
    addUnit("magneticfield", "uT", 1e-6);
    addUnit("magneticfield", "nT", 1e-9);
    addUnit("magneticfield", "G", 1e-4);
    addUnit("magneticfield", "mG", 1e-7);
    addUnit("magneticfield", "uG", 1e-10);
    addUnit("magneticfield", "nG", 1e-13);

    // pressure
    addUnit("pressure", "Pa", 1.);
    addUnit("pressure", "N/m2", 1.);
    addUnit("pressure", "J/m3", 1.);
    addUnit("pressure", "bar", 1e5);
    addUnit("pressure", "mbar", 1e2);
    addUnit("pressure", "hPa", 1e2);
    addUnit("pressure", "Ba", 1e-1);
    addUnit("pressure", "erg/cm3", 1e-1);
    addUnit("pressure", "K/m3", k);

    // bolometric luminosity
    addUnit("bolluminosity", "W", 1.);
    addUnit("bolluminosity", "erg/s", 1e-7);
    addUnit("bolluminosity", "Lsun", Lsun);

    // neutral monochromatic luminosity (lambda L_lambda = nu L_nu)
    addUnit("neutralmonluminosity", "W", 1.);
    addUnit("neutralmonluminosity", "erg/s", 1e-7);
    addUnit("neutralmonluminosity", "Lsun", Lsun);

    // neutral flux density (lambda F_lambda = nu F_nu)
    addUnit("neutralfluxdensity", "W/m2", 1.);
    addUnit("neutralfluxdensity", "erg/s/cm2", 1e-3);

    // neutral surface brightness (lambda f_lambda = nu f_nu)
    addUnit("neutralsurfacebrightness", "W/m2/sr", 1.);
    addUnit("neutralsurfacebrightness", "W/m2/arcsec2", 1. / pow(M_PI / (180. * 3600.), 2));
    addUnit("neutralsurfacebrightness", "erg/s/cm2/sr", 1e-3);
    addUnit("neutralsurfacebrightness", "erg/s/cm2/arcsec2", 1e-3 / pow(M_PI / (180. * 3600.), 2));

    // neutral mean intensity or spectral radiance (lambda J_lambda = nu J_nu)
    addUnit("neutralmeanintensity", "W/m2/sr", 1.);
    addUnit("neutralmeanintensity", "W/m2/arcsec2", 1. / pow(M_PI / (180. * 3600.), 2));
    addUnit("neutralmeanintensity", "erg/s/cm2/sr", 1e-3);
    addUnit("neutralmeanintensity", "erg/s/cm2/arcsec2", 1e-3 / pow(M_PI / (180. * 3600.), 2));

    // wavelength monochromatic luminosity (L_lambda)
    addUnit("wavelengthmonluminosity", "W/m", 1.);
    addUnit("wavelengthmonluminosity", "W/micron", 1e6);
    addUnit("wavelengthmonluminosity", "W/Angstrom", 1e10);
    addUnit("wavelengthmonluminosity", "erg/s/cm", 1e-5);
    addUnit("wavelengthmonluminosity", "erg/s/micron", 1e-1);
    addUnit("wavelengthmonluminosity", "erg/s/Angstrom", 1e3);
    addUnit("wavelengthmonluminosity", "Lsun/micron", Lsun * 1e6);

    // wavelength flux density (F_lambda)
    addUnit("wavelengthfluxdensity", "W/m3", 1.);
    addUnit("wavelengthfluxdensity", "W/m2/micron", 1e6);
    addUnit("wavelengthfluxdensity", "W/micron/m2", 1e6);
    addUnit("wavelengthfluxdensity", "W/m2/Angstrom", 1e10);
    addUnit("wavelengthfluxdensity", "W/Angstrom/m2", 1e10);
    addUnit("wavelengthfluxdensity", "erg/s/cm3", 1e-1);
    addUnit("wavelengthfluxdensity", "erg/s/cm2/micron", 1e3);
    addUnit("wavelengthfluxdensity", "erg/s/micron/cm2", 1e3);
    addUnit("wavelengthfluxdensity", "erg/s/cm2/Angstrom", 1e7);
    addUnit("wavelengthfluxdensity", "erg/s/Angstrom/cm2", 1e7);

    // wavelength surface brightness (f_lambda)
    addUnit("wavelengthsurfacebrightness", "W/m3/sr", 1.);
    addUnit("wavelengthsurfacebrightness", "W/m2/micron/sr", 1e6);
    addUnit("wavelengthsurfacebrightness", "W/micron/m2/sr", 1e6);
    addUnit("wavelengthsurfacebrightness", "W/m2/micron/arcsec2", 1e6 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthsurfacebrightness", "W/micron/m2/arcsec2", 1e6 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthsurfacebrightness", "W/m2/Angstrom/sr", 1e10);
    addUnit("wavelengthsurfacebrightness", "W/Angstrom/m2/sr", 1e10);
    addUnit("wavelengthsurfacebrightness", "W/m2/Angstrom/arcsec2", 1e10 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthsurfacebrightness", "W/Angstrom/m2/arcsec2", 1e10 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthsurfacebrightness", "erg/s/cm3/sr", 1e-1);
    addUnit("wavelengthsurfacebrightness", "erg/s/cm2/micron/sr", 1e3);
    addUnit("wavelengthsurfacebrightness", "erg/s/micron/cm2/sr", 1e3);
    addUnit("wavelengthsurfacebrightness", "erg/s/cm2/micron/arcsec2", 1e3 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthsurfacebrightness", "erg/s/micron/cm2/arcsec2", 1e3 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthsurfacebrightness", "erg/s/cm2/Angstrom/sr", 1e7);
    addUnit("wavelengthsurfacebrightness", "erg/s/Angstrom/cm2/sr", 1e7);
    addUnit("wavelengthsurfacebrightness", "erg/s/cm2/Angstrom/arcsec2", 1e7 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthsurfacebrightness", "erg/s/Angstrom/cm2/arcsec2", 1e7 / pow(M_PI / (180. * 3600.), 2));

    // wavelength mean intensity or spectral radiance (J_lambda)
    addUnit("wavelengthmeanintensity", "W/m3/sr", 1.);
    addUnit("wavelengthmeanintensity", "W/m2/micron/sr", 1e6);
    addUnit("wavelengthmeanintensity", "W/micron/m2/sr", 1e6);
    addUnit("wavelengthmeanintensity", "W/m2/micron/arcsec2", 1e6 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthmeanintensity", "W/micron/m2/arcsec2", 1e6 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthmeanintensity", "W/m2/Angstrom/sr", 1e10);
    addUnit("wavelengthmeanintensity", "W/Angstrom/m2/sr", 1e10);
    addUnit("wavelengthmeanintensity", "W/m2/Angstrom/arcsec2", 1e10 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthmeanintensity", "W/Angstrom/m2/arcsec2", 1e10 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthmeanintensity", "erg/s/cm3/sr", 1e-1);
    addUnit("wavelengthmeanintensity", "erg/s/cm2/micron/sr", 1e3);
    addUnit("wavelengthmeanintensity", "erg/s/micron/cm2/sr", 1e3);
    addUnit("wavelengthmeanintensity", "erg/s/cm2/micron/arcsec2", 1e3 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthmeanintensity", "erg/s/micron/cm2/arcsec2", 1e3 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthmeanintensity", "erg/s/cm2/Angstrom/sr", 1e7);
    addUnit("wavelengthmeanintensity", "erg/s/Angstrom/cm2/sr", 1e7);
    addUnit("wavelengthmeanintensity", "erg/s/cm2/Angstrom/arcsec2", 1e7 / pow(M_PI / (180. * 3600.), 2));
    addUnit("wavelengthmeanintensity", "erg/s/Angstrom/cm2/arcsec2", 1e7 / pow(M_PI / (180. * 3600.), 2));

    // frequency monochromatic luminosity (L_nu)
    addUnit("frequencymonluminosity", "W/Hz", 1.);
    addUnit("frequencymonluminosity", "erg/s/Hz", 1e-7);
    addUnit("frequencymonluminosity", "Lsun/Hz", Lsun);

    // frequency flux density (F_nu)
    addUnit("frequencyfluxdensity", "W/m2/Hz", 1.);
    addUnit("frequencyfluxdensity", "W/Hz/m2", 1.);
    addUnit("frequencyfluxdensity", "erg/s/cm2/Hz", 1e-3);
    addUnit("frequencyfluxdensity", "erg/s/Hz/cm2", 1e-3);
    addUnit("frequencyfluxdensity", "Jy", 1e-26);
    addUnit("frequencyfluxdensity", "mJy", 1e-29);
    addUnit("frequencyfluxdensity", "MJy", 1e-20);

    // frequency surface brightness (f_nu)
    addUnit("frequencysurfacebrightness", "W/m2/Hz/sr", 1.);
    addUnit("frequencysurfacebrightness", "W/Hz/m2/sr", 1.);
    addUnit("frequencysurfacebrightness", "W/m2/Hz/arcsec2", 1. / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencysurfacebrightness", "W/Hz/m2/arcsec2", 1. / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencysurfacebrightness", "erg/s/cm2/Hz/sr", 1e-3);
    addUnit("frequencysurfacebrightness", "erg/s/Hz/cm2/sr", 1e-3);
    addUnit("frequencysurfacebrightness", "erg/s/cm2/Hz/arcsec2", 1e-3 / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencysurfacebrightness", "erg/s/Hz/cm2/arcsec2", 1e-3 / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencysurfacebrightness", "Jy/sr", 1e-26);
    addUnit("frequencysurfacebrightness", "Jy/arcsec2", 1e-26 / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencysurfacebrightness", "MJy/sr", 1e-20);
    addUnit("frequencysurfacebrightness", "MJy/arcsec2", 1e-20 / pow(M_PI / (180. * 3600.), 2));

    // frequency mean intensity or spectral radiance (J_nu)
    addUnit("frequencymeanintensity", "W/m2/Hz/sr", 1.);
    addUnit("frequencymeanintensity", "W/Hz/m2/sr", 1.);
    addUnit("frequencymeanintensity", "W/m2/Hz/arcsec2", 1. / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencymeanintensity", "W/Hz/m2/arcsec2", 1. / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencymeanintensity", "erg/s/cm2/Hz/sr", 1e-3);
    addUnit("frequencymeanintensity", "erg/s/Hz/cm2/sr", 1e-3);
    addUnit("frequencymeanintensity", "erg/s/cm2/Hz/arcsec2", 1e-3 / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencymeanintensity", "erg/s/Hz/cm2/arcsec2", 1e-3 / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencymeanintensity", "Jy/sr", 1e-26);
    addUnit("frequencymeanintensity", "Jy/arcsec2", 1e-26 / pow(M_PI / (180. * 3600.), 2));
    addUnit("frequencymeanintensity", "MJy/sr", 1e-20);
    addUnit("frequencymeanintensity", "MJy/arcsec2", 1e-20 / pow(M_PI / (180. * 3600.), 2));

    // angular size (for objects in the sky)
    addUnit("angle", "rad", 1.);
    addUnit("angle", "deg", M_PI / 180.);
    addUnit("angle", "arcsec", M_PI / (180. * 3600.));

    // positioning angle (for instruments)
    addUnit("posangle", "rad", 1.);
    addUnit("posangle", "deg", M_PI / 180.);

    // solid angle
    addUnit("solidangle", "sr", 1.);
    addUnit("solidangle", "arcsec2", pow(M_PI / (180. * 3600.), 2));

    // *** add default units for each physical quantity and for each unit system

    // SI unit system
    addDefaultUnit("SIUnits", "length", "m");
    addDefaultUnit("SIUnits", "distance", "m");
    addDefaultUnit("SIUnits", "wavelength", "m");
    addDefaultUnit("SIUnits", "grainsize", "m");
    addDefaultUnit("SIUnits", "section", "m2");
    addDefaultUnit("SIUnits", "volume", "m3");
    addDefaultUnit("SIUnits", "velocity", "m/s");
    addDefaultUnit("SIUnits", "acceleration", "m/s2");
    addDefaultUnit("SIUnits", "mass", "kg");
    addDefaultUnit("SIUnits", "bulkmass", "kg");
    addDefaultUnit("SIUnits", "bulkmassdensity", "kg/m3");
    addDefaultUnit("SIUnits", "masssurfacedensity", "kg/m2");
    addDefaultUnit("SIUnits", "massvolumedensity", "kg/m3");
    addDefaultUnit("SIUnits", "massrate", "kg/s");
    addDefaultUnit("SIUnits", "numbersurfacedensity", "1/m2");
    addDefaultUnit("SIUnits", "numbervolumedensity", "1/m3");
    addDefaultUnit("SIUnits", "masscoefficient", "m2/kg");
    addDefaultUnit("SIUnits", "time", "s");
    addDefaultUnit("SIUnits", "temperature", "K");
    addDefaultUnit("SIUnits", "energy", "J");
    addDefaultUnit("SIUnits", "magneticfield", "T");
    addDefaultUnit("SIUnits", "pressure", "Pa");
    addDefaultUnit("SIUnits", "bolluminosity", "W");
    addDefaultUnit("SIUnits", "neutralmonluminosity", "W");
    addDefaultUnit("SIUnits", "neutralfluxdensity", "W/m2");
    addDefaultUnit("SIUnits", "neutralsurfacebrightness", "W/m2/sr");
    addDefaultUnit("SIUnits", "neutralmeanintensity", "W/m2/sr");
    addDefaultUnit("SIUnits", "wavelengthmonluminosity", "W/m");
    addDefaultUnit("SIUnits", "wavelengthfluxdensity", "W/m3");
    addDefaultUnit("SIUnits", "wavelengthsurfacebrightness", "W/m3/sr");
    addDefaultUnit("SIUnits", "wavelengthmeanintensity", "W/m3/sr");
    addDefaultUnit("SIUnits", "frequencymonluminosity", "W/Hz");
    addDefaultUnit("SIUnits", "frequencyfluxdensity", "W/m2/Hz");
    addDefaultUnit("SIUnits", "frequencysurfacebrightness", "W/m2/Hz/sr");
    addDefaultUnit("SIUnits", "frequencymeanintensity", "W/m2/Hz/sr");
    addDefaultUnit("SIUnits", "angle", "rad");
    addDefaultUnit("SIUnits", "posangle", "rad");
    addDefaultUnit("SIUnits", "solidangle", "sr");

    // stellar unit system
    addDefaultUnit("StellarUnits", "length", "AU");
    addDefaultUnit("StellarUnits", "distance", "pc");
    addDefaultUnit("StellarUnits", "wavelength", "micron");
    addDefaultUnit("StellarUnits", "grainsize", "micron");
    addDefaultUnit("StellarUnits", "section", "m2");
    addDefaultUnit("StellarUnits", "volume", "AU3");
    addDefaultUnit("StellarUnits", "velocity", "km/s");
    addDefaultUnit("StellarUnits", "acceleration", "m/s2");
    addDefaultUnit("StellarUnits", "mass", "Msun");
    addDefaultUnit("StellarUnits", "bulkmass", "kg");
    addDefaultUnit("StellarUnits", "bulkmassdensity", "kg/m3");
    addDefaultUnit("StellarUnits", "masssurfacedensity", "Msun/AU2");
    addDefaultUnit("StellarUnits", "massvolumedensity", "Msun/AU3");
    addDefaultUnit("StellarUnits", "massrate", "Msun/yr");
    addDefaultUnit("StellarUnits", "numbersurfacedensity", "1/cm3");
    addDefaultUnit("StellarUnits", "numbervolumedensity", "1/cm3");
    addDefaultUnit("StellarUnits", "masscoefficient", "m2/kg");
    addDefaultUnit("StellarUnits", "time", "Gyr");
    addDefaultUnit("StellarUnits", "temperature", "K");
    addDefaultUnit("StellarUnits", "energy", "J");
    addDefaultUnit("StellarUnits", "magneticfield", "uG");
    addDefaultUnit("StellarUnits", "pressure", "K/m3");
    addDefaultUnit("StellarUnits", "bolluminosity", "Lsun");
    addDefaultUnit("StellarUnits", "neutralmonluminosity", "Lsun");
    addDefaultUnit("StellarUnits", "neutralfluxdensity", "W/m2");
    addDefaultUnit("StellarUnits", "neutralsurfacebrightness", "W/m2/arcsec2");
    addDefaultUnit("StellarUnits", "neutralmeanintensity", "W/m2/sr");
    addDefaultUnit("StellarUnits", "wavelengthmonluminosity", "Lsun/micron");
    addDefaultUnit("StellarUnits", "wavelengthfluxdensity", "W/m2/micron");
    addDefaultUnit("StellarUnits", "wavelengthsurfacebrightness", "W/m2/micron/arcsec2");
    addDefaultUnit("StellarUnits", "wavelengthmeanintensity", "W/m2/micron/sr");
    addDefaultUnit("StellarUnits", "frequencymonluminosity", "W/Hz");
    addDefaultUnit("StellarUnits", "frequencyfluxdensity", "Jy");
    addDefaultUnit("StellarUnits", "frequencysurfacebrightness", "MJy/sr");
    addDefaultUnit("StellarUnits", "frequencymeanintensity", "W/m2/Hz/sr");
    addDefaultUnit("StellarUnits", "angle", "arcsec");
    addDefaultUnit("StellarUnits", "posangle", "deg");
    addDefaultUnit("StellarUnits", "solidangle", "arcsec2");

    // extra-galactic unit system
    addDefaultUnit("ExtragalacticUnits", "length", "pc");
    addDefaultUnit("ExtragalacticUnits", "distance", "Mpc");
    addDefaultUnit("ExtragalacticUnits", "wavelength", "micron");
    addDefaultUnit("ExtragalacticUnits", "grainsize", "micron");
    addDefaultUnit("ExtragalacticUnits", "section", "m2");
    addDefaultUnit("ExtragalacticUnits", "volume", "pc3");
    addDefaultUnit("ExtragalacticUnits", "velocity", "km/s");
    addDefaultUnit("ExtragalacticUnits", "acceleration", "m/s2");
    addDefaultUnit("ExtragalacticUnits", "mass", "Msun");
    addDefaultUnit("ExtragalacticUnits", "bulkmass", "kg");
    addDefaultUnit("ExtragalacticUnits", "bulkmassdensity", "kg/m3");
    addDefaultUnit("ExtragalacticUnits", "masssurfacedensity", "Msun/pc2");
    addDefaultUnit("ExtragalacticUnits", "massvolumedensity", "Msun/pc3");
    addDefaultUnit("ExtragalacticUnits", "massrate", "Msun/yr");
    addDefaultUnit("ExtragalacticUnits", "numbersurfacedensity", "1/cm2");
    addDefaultUnit("ExtragalacticUnits", "numbervolumedensity", "1/cm3");
    addDefaultUnit("ExtragalacticUnits", "masscoefficient", "m2/kg");
    addDefaultUnit("ExtragalacticUnits", "time", "Gyr");
    addDefaultUnit("ExtragalacticUnits", "temperature", "K");
    addDefaultUnit("ExtragalacticUnits", "energy", "J");
    addDefaultUnit("ExtragalacticUnits", "magneticfield", "uG");
    addDefaultUnit("ExtragalacticUnits", "pressure", "K/m3");
    addDefaultUnit("ExtragalacticUnits", "bolluminosity", "Lsun");
    addDefaultUnit("ExtragalacticUnits", "neutralmonluminosity", "Lsun");
    addDefaultUnit("ExtragalacticUnits", "neutralfluxdensity", "W/m2");
    addDefaultUnit("ExtragalacticUnits", "neutralsurfacebrightness", "W/m2/arcsec2");
    addDefaultUnit("ExtragalacticUnits", "neutralmeanintensity", "W/m2/sr");
    addDefaultUnit("ExtragalacticUnits", "wavelengthmonluminosity", "Lsun/micron");
    addDefaultUnit("ExtragalacticUnits", "wavelengthfluxdensity", "W/m2/micron");
    addDefaultUnit("ExtragalacticUnits", "wavelengthsurfacebrightness", "W/m2/micron/arcsec2");
    addDefaultUnit("ExtragalacticUnits", "wavelengthmeanintensity", "W/m2/micron/sr");
    addDefaultUnit("ExtragalacticUnits", "frequencymonluminosity", "W/Hz");
    addDefaultUnit("ExtragalacticUnits", "frequencyfluxdensity", "Jy");
    addDefaultUnit("ExtragalacticUnits", "frequencysurfacebrightness", "MJy/sr");
    addDefaultUnit("ExtragalacticUnits", "frequencymeanintensity", "W/m2/Hz/sr");
    addDefaultUnit("ExtragalacticUnits", "angle", "arcsec");
    addDefaultUnit("ExtragalacticUnits", "posangle", "deg");
    addDefaultUnit("ExtragalacticUnits", "solidangle", "arcsec2");
}

////////////////////////////////////////////////////////////////////
