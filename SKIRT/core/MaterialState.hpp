/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALSTATE_HPP
#define MATERIALSTATE_HPP

#include "MediumState.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the MaterialState class provides access to the medium state for a particular
    spatial cell and medium component, as determined during construction. This includes the common
    state variables for the cell as well as the specific state variables for the medium component.
    */
class MaterialState
{
    //============= Construction - Setup - Destruction =============

public:
    /** The constructor arguments specify the spatial cell and medium component to be represented,
        as well as a reference to the medium state object that holds the actual data. The constant
        reference is cast internally to a writable reference. This is ugly, but it avoids providing
        two versions of this class (read-only and writable). */
    MaterialState(const MediumState& ms, int m, int h) : _ms{const_cast<MediumState&>(ms)}, _m{m}, _h{h} {}

    //============= Setting =============

public:
    /** This function sets the number density \f$n\f$ of the medium component in the spatial cell.
        */
    void setNumberDensity(double value) { _ms.setNumberDensity(_m, _h, value); }

    /** This function sets the metallicity \f$Z\f$ of the medium component in the spatial cell. */
    void setMetallicity(double value) { _ms.setMetallicity(_m, _h, value); }

    /** This function sets the temperature \f$T\f$ of the medium component in the spatial cell. */
    void setTemperature(double value) { _ms.setTemperature(_m, _h, value); }

    /** This function sets the custom state variable with index \f$i\f$ of the medium component in
        the spatial cell. */
    void setCustom(int i, double value) { _ms.setCustom(_m, _h, i, value); }

    //============= Querying =============

public:
    /** This function returns the index \f$m\f$ of the spatial cell being represented. */
    int cellIndex() const { return _m; }

    /** This function returns the index \f$h\f$ of the medium component being represented. */
    int mediumIndex() const { return _h; }

    /** This function returns the volume \f$V\f$ of the spatial cell. */
    double volume() const { return _ms.volume(_m); }

    /** This function returns the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in
        the spatial cell. */
    Vec bulkVelocity() const { return _ms.bulkVelocity(_m); }

    /** This function returns the magnetic field \f${\boldsymbol{B}}\f$ in the spatial cell. */
    Vec magneticField() const { return _ms.magneticField(_m); }

    /** This function returns the number density \f$n\f$ of the medium component in the spatial
        cell. */
    double numberDensity() const { return _ms.numberDensity(_m, _h); }

    /** This function returns the metallicity \f$Z\f$ of the medium component in the spatial cell.
        */
    double metallicity() const { return _ms.metallicity(_m, _h); }

    /** This function returns the temperature \f$T\f$ of the medium component in the spatial cell.
        */
    double temperature() const { return _ms.temperature(_m, _h); }

    /** This function returns the custom state variable with index \f$i\f$ of the medium component
        in the spatial cell. */
    double custom(int i) const { return _ms.custom(_m, _h, i); }

    //======================== Data Members ========================

private:
    MediumState& _ms;
    int _m{-1};
    int _h{-1};
};

////////////////////////////////////////////////////////////////

#endif
