/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDMEDIUM_HPP
#define IMPORTEDMEDIUM_HPP

#include "MaterialMix.hpp"
#include "MaterialMixFamily.hpp"
#include "Medium.hpp"
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
        ATTRIBUTE_RELEVANT_IF(importVelocity, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(importVelocity, "Level2")
        ATTRIBUTE_INSERT(importVelocity, "importVelocity:MediumVelocity")

        PROPERTY_BOOL(importMagneticField, "import magnetic field components (3 columns)")
        ATTRIBUTE_DEFAULT_VALUE(importMagneticField, "false")
        ATTRIBUTE_DISPLAYED_IF(importMagneticField, "Level3")
        ATTRIBUTE_INSERT(importMagneticField, "importMagneticField:MagneticField")

        PROPERTY_BOOL(importVariableMixParams, "import parameter(s) to select a spatially varying material mix")
        ATTRIBUTE_DEFAULT_VALUE(importVariableMixParams, "false")
        ATTRIBUTE_DISPLAYED_IF(importVariableMixParams, "(!NonIdentitySpatialCellLibrary)&Level2")

        PROPERTY_STRING(useColumns, "a list of names corresponding to columns in the file to be imported")
        ATTRIBUTE_DEFAULT_VALUE(useColumns, "")
        ATTRIBUTE_REQUIRED_IF(useColumns, "false")
        ATTRIBUTE_DISPLAYED_IF(useColumns, "Level3")

        PROPERTY_ITEM(materialMix, MaterialMix, "the material type and properties throughout the medium")
        ATTRIBUTE_DEFAULT_VALUE(materialMix, "MeanInterstellarDustMix")
        ATTRIBUTE_RELEVANT_IF(materialMix, "!importVariableMixParams")

        PROPERTY_ITEM(materialMixFamily, MaterialMixFamily,
                      "the material mix family addressed by the imported parameter(s)")
        ATTRIBUTE_DEFAULT_VALUE(materialMixFamily, "SelectDustMixFamily")
        ATTRIBUTE_RELEVANT_IF(materialMixFamily, "importVariableMixParams")

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
        medium at the specified position.

        If the \em importVariableMixParams flag is enabled, the appropriate material mix is
        selected from the configured material mix family based on the value of the imported
        parameters for the specified position. If the flag is disabled, the fixed configured
        material mix is returned regardless of position. */
    const MaterialMix* mix(Position bfr) const override;

    /** This function returns a default MaterialMix object representative of the material
        properties of the medium.

        If the \em importVariableMixParams flag is enabled, a default material mix is selected from
        the configured material mix family using the appropriate number of parameters with a value
        of zero. If the flag is disabled, the fixed configured material mix is returned. */
    const MaterialMix* mix() const override;

    /** This function returns the configured value of the \em importVariableMixParams flag. If
        true, this medium may return a different MaterialMix object depending on the specified
        position. If false, the same object is always returned. */
    bool hasVariableMix() const override;

    /** This function returns true if the \em importVelocity flag is enabled for the medium. */
    bool hasVelocity() const override;

    /** This function returns the bulk velocity of the medium at the specified position. If the \em
        importVelocity flag is enabled, it simply calls the corresponding function in the snapshot
        object; otherwise it returns zero velocity. */
    Vec bulkVelocity(Position bfr) const override;

    /** This function returns true if the \em importMagneticField flag is enabled for the medium. */
    bool hasMagneticField() const override;

    /** This function returns the magnetic field vector of the medium at the specified position. If
        the \em importMagneticField flag is enabled, it simply calls the corresponding function in
        the snapshot object; otherwise it returns a zero magnetic field. */
    Vec magneticField(Position bfr) const override;

    /** This function returns true if the medium has a gas or electron material mix and the \em
        importMetallicity flag is enabled for the medium. */
    bool hasMetallicity() const override;

    /** This function returns the metallicity of the medium at the specified position, if defined
        in the input model. Specifically, if the medium has a gas or electron material mix and the
        \em importMetallicity flag is enabled, the function returns the imported metallicity at the
        specified position. Otherwise, the function returns zero. */
    double metallicity(Position bfr) const override;

    /** This function returns true if the medium has a gas or electron material mix and the \em
        importTemperature flag is enabled for the medium. */
    bool hasTemperature() const override;

    /** This function returns the temperature of the medium at the specified position, if defined
        in the input model. Specifically, if the medium has a gas or electron material mix and the
        \em importTemperature flag is enabled, the function returns the imported temperature at the
        specified position. Otherwise, the function returns zero. */
    double temperature(Position bfr) const override;

    /** This function returns true if custom input model parameters are available for this medium.
        See the parameters() function for more information. */
    bool hasParameters() const override;

    /** If custom input model parameters are available for this medium, this function stores the
        parameter values at the specified position into the given array. If the position is outside
        the domain, the parameter values default to zero. If no custom input model parameters are
        available for this medium, the array is resized to zero length.

        For an imported medium that is configured with a material mix, custom input model
        parameters are imported from the snapshot as requested by the MaterialMix::parameterInfo()
        function. The number and order of parameters returned by this function then reflects the
        number and order of parameter descriptions returned by the MaterialMix::parameterInfo()
        function. If the medium is configured with a material mix \em family, any imported
        parameters are used to select a member of the material mix family at each location, and
        this function returns an empty array. In other words, requests by material mix family
        members for custom input model parameters are ignored. */
    void parameters(Position bfr, Array& params) const override;

    /** This function returns the number density of the medium at the specified position. */
    double numberDensity(Position bfr) const override;

    /** This function returns the total number of material entities in the medium. The function
        uses the default material mix (the one at the origin) throughout the complete spatial
        domain; if the \em importVariableMixParams flag is enabled, this is an approximation. */
    double number() const override;

    /** This function returns the mass density of the medium at the specified position. */
    double massDensity(Position bfr) const override;

    /** This function returns the total mass in the medium. The function uses the default material
        mix (the one at the origin) throughout the complete spatial domain; if the \em
        importVariableMixParams flag is enabled, this is an approximation. */
    double mass() const override;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$ along the
        full X axis of the model coordinate system. The function uses the default material mix (the
        one at the origin) throughout the complete spatial domain; if the \em
        importVariableMixParams flag is enabled, this is an approximation. */
    double opticalDepthX(double lambda) const override;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$ along the
        full Y axis of the model coordinate system. The function uses the default material mix (the
        one at the origin) throughout the complete spatial domain; if the \em
        importVariableMixParams flag is enabled, this is an approximation. */
    double opticalDepthY(double lambda) const override;

    /** This function returns the optical depth of the medium at wavelength \f$\lambda\f$ along the
        full Z axis of the model coordinate system. The function uses the default material mix (the
        one at the origin) throughout the complete spatial domain; if the \em
        importVariableMixParams flag is enabled, this is an approximation. */
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

    /** This function returns (a pointer to) the snapshot object associated with this imported
        medium. It is intended to provide InputModelProbe instances with direct access to the
        snapshot for probing imported information that is not otherwise made available to the
        simulation. To preserve proper data encapsulation, this function should \em not be called
        from anywhere else in the simulation machinery. */
    const Snapshot* snapshot() const;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Snapshot* _snapshot{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
