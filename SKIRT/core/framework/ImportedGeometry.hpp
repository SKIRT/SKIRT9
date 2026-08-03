/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDGEOMETRY_HPP
#define IMPORTEDGEOMETRY_HPP

#include "GenGeometry.hpp"
#include "SiteListInterface.hpp"
class Snapshot;

////////////////////////////////////////////////////////////////////

/** ImportedGeometry is an abstract class for describing 3D geometries with a spatial density
    distribution imported from an input file. The input data is usually derived from a
    hydrodynamical simulation snapshot. Various types of snapshots are supported by subclasses of
    this class. Refer to the subclass documentation for information on the file format. */
class ImportedGeometry : public GenGeometry, public SiteListInterface
{
    ITEM_ABSTRACT(ImportedGeometry, GenGeometry, "a geometry imported from snapshot data")
        ATTRIBUTE_TYPE_INSERT(ImportedGeometry, "SiteListInterface")

        PROPERTY_STRING(filename, "the name of the file to be imported")

        ATTRIBUTE_SUB_PROPERTIES_HERE(ImportedGeometry)

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

        PROPERTY_STRING(useColumns, "a list of names corresponding to columns in the file to be imported")
        ATTRIBUTE_DEFAULT_VALUE(useColumns, "")
        ATTRIBUTE_REQUIRED_IF(useColumns, "false")
        ATTRIBUTE_DISPLAYED_IF(useColumns, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function imports the snapshot data from the input file through a Snapshot object of
        the appropriate type. Specifically, it first calls the createSnapshot() function, which
        must be implemented in a subclass, to construct and open a Snapshot object of the appropriate type.
        It then passes the user-configurable options of this class to the Snapshot object and tells
        it to import the data. */
    void setupSelfAfter() override;

    /** This function constructs a new Snapshot object of the type appropriate for the subclass,
        calls its open() function, configures the mass or density column, and returns a pointer to
        the object. Ownership of the Snapshot object is transferred to the caller. */
    virtual Snapshot* createAndOpenSnapshot() = 0;

    /** The destructor deletes the snapshot object, if present. */
    ~ImportedGeometry();

    //======================== Other Functions =======================

public:
    /** This function returns the density for this geometry at the given position. It simply calls
        the corresponding function in the snapshot object. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry's spatial density distribution.
        It simply calls the corresponding function in the snapshot object. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density of the geometry, defined as the
        integration of the density along the entire X-axis. It simply calls the corresponding
        function in the snapshot object. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density of the geometry, defined as the
        integration of the density along the entire Y-axis. It simply calls the corresponding
        function in the snapshot object. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density of the geometry, defined as the
        integration of the density along the entire Z-axis. It simply calls the corresponding
        function in the snapshot object. */
    double SigmaZ() const override;

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
    double _norm{0.};
};

////////////////////////////////////////////////////////////////////

#endif
