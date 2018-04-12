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
    addUnit("length", "km", 1e3);
    addUnit("length", "AU", AU);
    addUnit("length", "pc", pc);
    addUnit("length", "kpc", 1e3 * pc);
    addUnit("length", "Mpc", 1e6 * pc);

    // distance
    addUnit("distance", "m", 1.);
    addUnit("distance", "cm", 1e-2);
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
    addUnit("wavelength", "A", 1e-10);

    // grainsize
    addUnit("grainsize", "m", 1.);
    addUnit("grainsize", "cm", 1e-2);
    addUnit("grainsize", "mm", 1e-3);
    addUnit("grainsize", "micron", 1e-6);
    addUnit("grainsize", "nm", 1e-9);
    addUnit("grainsize", "A", 1e-10);

    // cross section
    addUnit("section", "m2", 1.);

    // volume
    addUnit("volume", "m3", 1.);
    addUnit("volume", "AU3", pow(AU,3));
    addUnit("volume", "pc3", pow(pc,3));

    // velocity
    addUnit("velocity", "m/s", 1.);
    addUnit("velocity", "km/s", 1e3);

    // mass
    addUnit("mass", "kg", 1.);
    addUnit("mass", "g", 1e-3);
    addUnit("mass", "Msun", Msun);

    // bulk mass
    addUnit("bulkmass", "kg", 1.);

    // bulk mass density
    addUnit("bulkmassdensity", "kg/m3", 1.);
    addUnit("bulkmassdensity", "g/cm3", 1e3);

    // mass surface density
    addUnit("masssurfacedensity", "kg/m2", 1.);
    addUnit("masssurfacedensity", "Msun/AU2", Msun/pow(AU,2));
    addUnit("masssurfacedensity", "Msun/pc2", Msun/pow(pc,2));

    // mass volume density
    addUnit("massvolumedensity", "kg/m3", 1.);
    addUnit("massvolumedensity", "g/cm3", 1e3);
    addUnit("massvolumedensity", "Msun/AU3", Msun/pow(AU,3));
    addUnit("massvolumedensity", "Msun/pc3", Msun/pow(pc,3));

    // opacity
    addUnit("opacity", "m2/kg", 1.);

    // time
    addUnit("time", "s", 1.);
    addUnit("time", "yr", year);
    addUnit("time", "Myr", 1e6*year);
    addUnit("time", "Gyr", 1e9*year);

    // temperature
    addUnit("temperature", "K", 1.);

    // energy
    addUnit("energy", "J", 1.);

    // pressure
    addUnit("pressure", "Pa", 1.);
    addUnit("pressure", "K/m3", k);

    // bolometric luminosity
    addUnit("bolluminosity", "W", 1.);
    addUnit("bolluminosity", "Lsun", Lsun);

    // monochromatic luminosity
    addUnit("monluminosity", "W/m", 1.);
    addUnit("monluminosity", "W/micron", 1e6);
    addUnit("monluminosity", "Lsun/micron", Lsun * 1e6);

    // neutral flux density (lambda F_lambda = nu F_nu)
    addUnit("neutralfluxdensity", "W/m2", 1.);

    // neutral surface brightness (lambda f_lambda = nu f_nu)
    addUnit("neutralsurfacebrightness", "W/m2/sr", 1.);
    addUnit("neutralsurfacebrightness", "W/m2/arcsec2", 1. / pow(M_PI/(180.*3600.),2));

    // wavelength flux density (F_lambda)
    addUnit("wavelengthfluxdensity", "W/m3", 1.);
    addUnit("wavelengthfluxdensity", "W/m2/micron", 1e6);

    // wavelength surface brightness (f_lambda)
    addUnit("wavelengthsurfacebrightness", "W/m3/sr", 1.);
    addUnit("wavelengthsurfacebrightness", "W/m2/micron/sr", 1e6);
    addUnit("wavelengthsurfacebrightness", "W/m2/micron/arcsec2", 1e6 / pow(M_PI/(180.*3600.),2));

    // frequency flux density (F_nu)
    addUnit("frequencyfluxdensity", "W/m2/Hz", 1.);
    addUnit("frequencyfluxdensity", "Jy", 1e-26);
    addUnit("frequencyfluxdensity", "mJy", 1e-29);
    addUnit("frequencyfluxdensity", "MJy", 1e-20);

    // frequency surface brightness (f_nu)
    addUnit("frequencysurfacebrightness", "W/m2/Hz/sr", 1.);
    addUnit("frequencysurfacebrightness", "W/m2/Hz/arcsec2", 1. / pow(M_PI/(180.*3600.),2));
    addUnit("frequencysurfacebrightness", "Jy/sr", 1e-26);
    addUnit("frequencysurfacebrightness", "Jy/arcsec2", 1e-26 / pow(M_PI/(180.*3600.),2));
    addUnit("frequencysurfacebrightness", "MJy/sr", 1e-20);
    addUnit("frequencysurfacebrightness", "MJy/arcsec2", 1e-20 / pow(M_PI/(180.*3600.),2));

    // angular size (for objects in the sky)
    addUnit("angle", "rad", 1.);
    addUnit("angle", "deg", M_PI/180.);
    addUnit("angle", "arcsec", M_PI/(180.*3600.));

    // positioning angle (for instruments)
    addUnit("posangle", "rad", 1.);
    addUnit("posangle", "deg", M_PI/180.);

    // solid angle
    addUnit("solidangle", "sr", 1.);
    addUnit("solidangle", "arcsec2", pow(M_PI/(180.*3600.),2));

    // *** add default units for each physical quantity and for each unit system

    // SI unit system
    addDefaultUnit("SIUnits", "length", "m");
    addDefaultUnit("SIUnits", "distance", "m");
    addDefaultUnit("SIUnits", "wavelength", "m");
    addDefaultUnit("SIUnits", "grainsize", "m");
    addDefaultUnit("SIUnits", "section", "m2");
    addDefaultUnit("SIUnits", "volume", "m3");
    addDefaultUnit("SIUnits", "velocity", "m/s");
    addDefaultUnit("SIUnits", "mass", "kg");
    addDefaultUnit("SIUnits", "bulkmass", "kg");
    addDefaultUnit("SIUnits", "bulkmassdensity", "kg/m3");
    addDefaultUnit("SIUnits", "masssurfacedensity", "kg/m2");
    addDefaultUnit("SIUnits", "massvolumedensity", "kg/m3");
    addDefaultUnit("SIUnits", "opacity", "m2/kg");
    addDefaultUnit("SIUnits", "time", "s");
    addDefaultUnit("SIUnits", "temperature", "K");
    addDefaultUnit("SIUnits", "energy", "J");
    addDefaultUnit("SIUnits", "pressure", "Pa");
    addDefaultUnit("SIUnits", "bolluminosity", "W");
    addDefaultUnit("SIUnits", "monluminosity", "W/m");
    addDefaultUnit("SIUnits", "neutralfluxdensity", "W/m2");
    addDefaultUnit("SIUnits", "neutralsurfacebrightness", "W/m2/sr");
    addDefaultUnit("SIUnits", "wavelengthfluxdensity", "W/m3");
    addDefaultUnit("SIUnits", "wavelengthsurfacebrightness", "W/m3/sr");
    addDefaultUnit("SIUnits", "frequencyfluxdensity", "W/m2/Hz");
    addDefaultUnit("SIUnits", "frequencysurfacebrightness", "W/m2/Hz/sr");
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
    addDefaultUnit("StellarUnits", "mass", "Msun");
    addDefaultUnit("StellarUnits", "bulkmass", "kg");
    addDefaultUnit("StellarUnits", "bulkmassdensity", "kg/m3");
    addDefaultUnit("StellarUnits", "masssurfacedensity", "Msun/AU2");
    addDefaultUnit("StellarUnits", "massvolumedensity", "Msun/AU3");
    addDefaultUnit("StellarUnits", "opacity", "m2/kg");
    addDefaultUnit("StellarUnits", "time", "Gyr");
    addDefaultUnit("StellarUnits", "temperature", "K");
    addDefaultUnit("StellarUnits", "energy", "J");
    addDefaultUnit("StellarUnits", "pressure", "K/m3");
    addDefaultUnit("StellarUnits", "bolluminosity", "Lsun");
    addDefaultUnit("StellarUnits", "monluminosity", "Lsun/micron");
    addDefaultUnit("StellarUnits", "neutralfluxdensity", "W/m2");
    addDefaultUnit("StellarUnits", "neutralsurfacebrightness", "W/m2/arcsec2");
    addDefaultUnit("StellarUnits", "wavelengthfluxdensity", "W/m2/micron");
    addDefaultUnit("StellarUnits", "wavelengthsurfacebrightness", "W/m2/micron/arcsec2");
    addDefaultUnit("StellarUnits", "frequencyfluxdensity", "Jy");
    addDefaultUnit("StellarUnits", "frequencysurfacebrightness", "MJy/sr");
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
    addDefaultUnit("ExtragalacticUnits", "mass", "Msun");
    addDefaultUnit("ExtragalacticUnits", "bulkmass", "kg");
    addDefaultUnit("ExtragalacticUnits", "bulkmassdensity", "kg/m3");
    addDefaultUnit("ExtragalacticUnits", "masssurfacedensity", "Msun/pc2");
    addDefaultUnit("ExtragalacticUnits", "massvolumedensity", "Msun/pc3");
    addDefaultUnit("ExtragalacticUnits", "opacity", "m2/kg");
    addDefaultUnit("ExtragalacticUnits", "time", "Gyr");
    addDefaultUnit("ExtragalacticUnits", "temperature", "K");
    addDefaultUnit("ExtragalacticUnits", "energy", "J");
    addDefaultUnit("ExtragalacticUnits", "pressure", "K/m3");
    addDefaultUnit("ExtragalacticUnits", "bolluminosity", "Lsun");
    addDefaultUnit("ExtragalacticUnits", "monluminosity", "Lsun/micron");
    addDefaultUnit("ExtragalacticUnits", "neutralfluxdensity", "W/m2");
    addDefaultUnit("ExtragalacticUnits", "neutralsurfacebrightness", "W/m2/arcsec2");
    addDefaultUnit("ExtragalacticUnits", "wavelengthfluxdensity", "W/m2/micron");
    addDefaultUnit("ExtragalacticUnits", "wavelengthsurfacebrightness", "W/m2/micron/arcsec2");
    addDefaultUnit("ExtragalacticUnits", "frequencyfluxdensity", "Jy");
    addDefaultUnit("ExtragalacticUnits", "frequencysurfacebrightness", "MJy/sr");
    addDefaultUnit("ExtragalacticUnits", "angle", "arcsec");
    addDefaultUnit("ExtragalacticUnits", "posangle", "deg");
    addDefaultUnit("ExtragalacticUnits", "solidangle", "arcsec2");
}

////////////////////////////////////////////////////////////////////
