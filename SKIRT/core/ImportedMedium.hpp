/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDMEDIUM_HPP
#define IMPORTEDMEDIUM_HPP

#include "Medium.hpp"
#include "MaterialMix.hpp"
#include "SiteListInterface.hpp"
class Snapshot;

////////////////////////////////////////////////////////////////////

/** ImportedMedium is an abstract class representing transfer media for which the spatial material
    density distribution (and, optionally, related properties such as the bulk velocity) is
    imported from an input file. The input data is usually derived from a hydrodynamical simulation
    snapshot. Various types of snapshots are supported by subclasses of this class. Refer to the
    subclass documentation for information on the file format. */
class ImportedMedium : public Medium, public SiteListInterface
{
    ITEM_ABSTRACT(ImportedMedium, Medium, "a transfer medium imported from snapshot data")
        ATTRIBUTE_TYPE_INSERT(ImportedMedium, "Dimension3,SiteListInterface")

    PROPERTY_STRING(filename, "the name of the file to be imported")

    ATTRIBUTE_SUB_PROPERTIES_HERE(ImportedMedium)

    PROPERTY_DOUBLE(massFraction, "the fraction of the mass to be included (or one to include all)")
        ATTRIBUTE_MIN_VALUE(massFraction, "[0")
        ATTRIBUTE_MAX_VALUE(massFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(massFraction, "1")

    PROPERTY_BOOL(importMetallicity, "import a metallicity column")
        ATTRIBUTE_DEFAULT_VALUE(importMetallicity, "false")

    PROPERTY_BOOL(importTemperature, "import a temperature column")
        ATTRIBUTE_DEFAULT_VALUE(importTemperature, "false")

    PROPERTY_DOUBLE(maxTemperature, "the maximum temperature for included mass (or zero to include all)")
        ATTRIBUTE_QUANTITY(maxTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(maxTemperature, "[0 K")
        ATTRIBUTE_MAX_VALUE(maxTemperature, "1000000 K]")
        ATTRIBUTE_DEFAULT_VALUE(maxTemperature, "0 K")
        ATTRIBUTE_RELEVANT_IF(maxTemperature, "importTemperature")

    PROPERTY_BOOL(importVelocity, "import velocity components (3 columns)")
        ATTRIBUTE_DEFAULT_VALUE(importVelocity, "false")
        ATTRIBUTE_DISPLAYED_IF(importVelocity, "(Panchromatic&Level2)|Level3")

    PROPERTY_ITEM(materialMix, MaterialMix, "the material type and properties throughout the medium")
        ATTRIBUTE_DEFAULT_VALUE(materialMix, "MeanInterstellarDustMix")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function imports the snapshot data from the input file through a Snapshot object of
        the appropriate type. Specifically, it first calls the createSnapshot() function, which
        must be implemented in a subclass, to construct and open a Snapshot object of the
        appropriate type. It then passes the user-configurable options of this class to the
        Snapshot object and tells it to import the data. */
    void setupSelfAfter() override;

    /** This function constructs a new Snapshot object of the type appropriate for the subclass,
        calls its open() function, configures the mass or density column, and returns a pointer to
        the object. Ownership of the Snapshot object is transferred to the caller. */
    virtual Snapshot* createAndOpenSnapshot() = 0;

    /** The destructor deletes the snapshot object, if present. */
    ~ImportedMedium();

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the medium, which is always 3 for an imported
        medium. */
    int dimension() const override;

    /** This function returns the MaterialMix object defining the material properties for the
        medium at the specified position. In the current implementation, the same object is
        returned regardless of position. This may change in the future. */
    const MaterialMix* mix(Position bfr = Position()) const override;

    /** This function returns the bulk velocity of the medium at the specified position. If the \em
        importVelocity flag is enabled, it simply calls the corresponding function in the snapshot
        object; otherwise is returns zero velocity. */
    Vec bulkVelocity(Position bfr) const override;

    /** This function returns the number density of the medium at the specified position. */
    double numberDensity(Position bfr) const override;

    /** This function returns the total number of material entities in the medium. */
    double number() const override;

    /** This function returns the mass density of the medium at the specified position. */
    double massDensity(Position bfr) const override;

    /** This function returns the total mass in the medium. */
    double mass() const override;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$
        along the full X axis of the model coordinate system. */
    double opticalDepthX(double lambda) const override;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$
        along the full Y axis of the model coordinate system. */
    double opticalDepthY(double lambda) const override;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$
        along the full Z axis of the model coordinate system. */
    double opticalDepthZ(double lambda) const override;

    /** This function generates a random position sampled from the medium's spatial density
        distribution. It simply calls the corresponding function in the snapshot object. In the
        current implementation, the conversion from number to mass is the same throughout the
        medium's spatial domain, meaning that there is no difference between sampling from the
        number density or the mass density. This may change in the future. */
    Position generatePosition() const override;

    /** This function returns the number of entities (particles or cells) used for defining the
        density distribution represented by the snapshot. The function is part of the
        SiteListInterface. It simply calls the corresponding function in the snapshot object. */
    int numSites() const override;

    /** This function returns the coordinates of the entity (particle or cell) with the specified
        zero-based index. If the index is out of range, the behavior is undefined. The function is
        part of the SiteListInterface. It simply calls the corresponding function in the snapshot
        object. */
    Position sitePosition(int index) const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Snapshot* _snapshot{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
