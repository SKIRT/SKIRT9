/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef UPDATE_STATUS_HPP
#define UPDATE_STATUS_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** UpdateStatus is a helper class for representing a tri-state enumeration value that indicates
    the update status of an entity (which can be a single value or an aggregate of values): \em
    NotUpdated means the entity has not been changed and thus has converged; \em UpdatedConverged
    means the entity has been changed but can be considered to have converged; \em
    UpdatedNotConverged means the entity has been changed and has not yet converged.

    This tri-state is implemented as a class rather than a straight enumeration to facilitate
    streamlined operations and queries on its value, as can be seen from the various functions
    defined for this class. UpdateStatus instances can be copied or moved at will. The tri-state
    value is stored in a single byte to limit the size of arrays of UpdateStatus instances. */
class UpdateStatus final
{
    //======================== Data Members ========================

private:
    // bit_0 = 1 indicates updated; bit_1 = 1 indicates not converged
    // this allows status values to be ORed when aggregating them
    const static uint8_t NotUpdated = 0;
    const static uint8_t UpdatedConverged = 1;
    const static uint8_t UpdatedNotConverged = 3;
    uint8_t _status{NotUpdated};

    //============= Construction =============

public:
    /** This is the default constructor; the status is initialized to \em NotUpdated. */
    UpdateStatus() {}

    //============= Setters =============

    /** If the current status is \em NotUpdated, this function changes it to \em UpdatedConverged.
        Otherwise, the function has no effect. Specifically, an \em UpdatedNotConverged status is
        never "downgraded" to \em UpdatedConverged. */
    void updateConverged() { _status |= UpdatedConverged; }

    /** This function sets the status to \em UpdatedNotConverged regardless of its current value.
        */
    void updateNotConverged() { _status = UpdatedNotConverged; }

    /** This function updates the receiving status according the specified other status, which is
        useful when determining the aggregated status for a set of entities. Specifically, if the
        other status is \em UpdatedNotConverged, the effect is that of the updateNotConverged()
        function; if the other status is \em UpdatedConverged, the effect is that of the
        updateConverged() function; and if the other status is \em NotUpdated, there is no effect.
        */
    void update(UpdateStatus other) { _status |= other._status; }

    //============= Queries =============

    /** This function returns true if the status is \em UpdatedConverged or \em
        UpdatedNotConverged, and false if the status is \em NotUpdated. */
    bool isUpdated() const { return _status != NotUpdated; }

    /** This function returns true if the status is \em NotUpdated or \em UpdatedConverged, and
        false if the status is \em UpdatedNotConverged. */
    bool isConverged() const { return _status != UpdatedNotConverged; }
};

//////////////////////////////////////////////////////////////////////

#endif
