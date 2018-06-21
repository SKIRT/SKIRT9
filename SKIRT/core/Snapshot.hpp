/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SNAPSHOT_HPP
#define SNAPSHOT_HPP

#include "Position.hpp"
#include "SnapshotParameter.hpp"
class SimulationItem;
class TextInFile;

////////////////////////////////////////////////////////////////////

/** Snapshot is an abstract base class for representing snapshot data that is imported from a
    column text file and is used to define the properties of primary sources or transfer media.
    This base class offers facilities to handle the input options, such as determining which
    columns are to be imported from the text file, and defines a common interface for all snapshot
    types, for example to obtain the density at a given position. The actual implementation must be
    provided in a derived class for each snapshot type.

    Some snapshot types use smoothed particles and some use spatial cells as their basic
    constituents. To avoid mentioning both, we use the generic term \em entity for referring to
    either particles or cells.

    A Snapshot instance must be initialized using the following strict calling sequence; calling
    any of the getters before the initialization sequence has completed results in undefined
    behavior:

      - open the input file by calling the open() function
      - configure the input file columns by calling the applicable importXXX() functions in the
        order the columns are expected to appear in the text file;
      - configure other options by calling the applicable setXXX() functions;
      - read the file calling the readAndClose() function.

    The functions for extracting information from a snapshot once it is loaded are thread-safe.
    However, the functions for configuring options and reading the snapshot file are \em not
    thread-safe; they must be called from a single thread (for each snapshot object). */
class Snapshot
{
    //================= Construction - Destruction =================

public:
    /** The default constructor initializes the snapshot in an invalid state; see the description
        of the required calling sequence in the class header. */
    Snapshot();

    /** The destructor deletes the input file object, if it would still be present. The fact that
        the destructor is virtual ensures that derived classes can be propertly destructed through
        a pointer to this base class. */
    virtual ~Snapshot();

    //========== Managing the input file ==========

public:
    /** This function creates an input file object corresponding to the specified file and opens it
        for reading; if the file can't be opened, a FatalError is thrown. It must be called \em
        before invoking any of the configuration functions(). This function takes several
        arguments: (1) \em item specifies a simulation item in the hierarchy of the caller (usually
        the caller itself) used to retrieve the input file path and an appropriate logger; (2) \em
        filename specifies the name of the file, including filename extension but excluding path
        and simulation prefix; (3) \em description describes the contents of the file for use in
        the log message issued after the file is successfully opened.

        The function is virtual so that subclasses can override it, for example to perform some
        initial configuration after opening but before client configuration occurs. */
    virtual void open(const SimulationItem* item, string filename, string description);

    /** This function reads the snapshot data from the input file, honoring the options set through
        the configuration functions, stores the data for later use, and closes the file.

        The implementation in this base class simply closes the file and deletes the corresponding
        file object. Subclasses must override the function to actually read the data and then call
        the implementation in this base class. */
    virtual void readAndClose();

protected:
    /** This function returns a pointer to the input file object. It is intended for use in
        subclasses. */
    TextInFile* infile() { return _infile; }

    //========== Configuring options ==========

public:
    /** This function configures the snapshot to import a spatial position with three components
        (x,y,z). */
    void importPosition();

    /** This function configures the snapshot to import a spatial radial size. */
    void importSize();

    /** This function configures the snapshot to import a mass density per unit of volume. The
        importMass() and importDensity() options are mutually exclusive; calling both functions for
        the same snapshot results in undefined behavior. */
    void importDensity();

    /** This function configures the snapshot to import a mass (i.e. density integrated over
        volume). The importMass() and importDensity() options are mutually exclusive; calling both
        functions for the same snapshot results in undefined behavior. */
    void importMass();

    /** This function configures the snapshot to import a metallicity fraction. */
    void importMetallicity();

    /** This function configures the snapshot to import a temperature. */
    void importTemperature();

    /** This function configures the snapshot to import a velocity with three components (x,y,z).
        */
    void importVelocity();

    /** This function configures the snapshot to import alignment characteristics including
        orientation and fraction of alignment. TO DO: provide details and implement. */
    void importAlignment();

    /** This function configures the snapshot to import a sequence of parameters as described by
        the specified list of SnapshotParameter metadata objects. */
    void importParameters(const vector<SnapshotParameter>& parameters);

    /** This function sets the policy for calculating the mass or mass density at a given spatial
        position in the snapshot. If the policy is not set during configuration, or if neither of
        the importMass() and importDensity() functions was invoked, attempting to obtain mass or
        mass density values results in undefined behavior.

        Otherwise, if the \em maxTemperature argument is positive, and the temperature field is
        being imported and its value for a given entity is positive as well, and the imported
        temperature is above the maximum temperature, then the density corresponding to this entity
        is considered to be zero.

        Otherwise, the intrinsic mass or mass density is obtained from the mass or density fields
        using a mechanism that depends on the snapshot type and possibly on the configuration.
        Then, this value is multiplied by the metallicity fraction if this field has been imported,
        and then it is multiplied by the value of the \em multiplier argument to obtain the final
        result. */
    void setMassDensityPolicy(double multiplier = 1., double maxTemperature = 0.);

protected:
    /** This function returns the column index of the first position field, or -1 if this is not
        being imported, for use by subclasses. */
    int positionIndex() { return _positionIndex; }

    /** This function returns the column index of the size field, or -1 if this is not being
        imported, for use by subclasses. */
    int sizeIndex() { return _sizeIndex; }

    /** This function returns the column index of the density field, or -1 if this is not being
        imported, for use by subclasses. */
    int densityIndex() { return _densityIndex; }

    /** This function returns the column index of the mass field, or -1 if this is not being
        imported, for use by subclasses. */
    int massIndex() { return _massIndex; }

    /** This function returns the column index of the metallicity field, or -1 if this is not being
        imported, for use by subclasses. */
    int metallicityIndex() { return _metallicityIndex; }

    /** This function returns the column index of the temperature field, or -1 if this is not being
        imported, for use by subclasses. */
    int temperatureIndex() { return _temperatureIndex; }

    /** This function returns the column index of the first velocity field, or -1 if this is not
        being imported, for use by subclasses. */
    int velocityIndex() { return _velocityIndex; }

    /** This function returns the column index of the first field in the parameter list, or -1 if
        this is not being imported, for use by subclasses. */
    int parametersIndex() { return _parametersIndex; }

    /** This function returns the mass or mass density multiplier configured by the user, or zero
        if the user did not configure the mass or mass density policy, for use by subclasses. */
    double multiplier() { return _multiplier; }

    /** This function returns the maximum temperature configured by the user for an entity to have
        mass, or zero if the user did not configure the mass or mass density policy, for use by
        subclasses. */
    double maxTemperature() { return _maxTemperature; }

    //======================== Data Members ========================

private:
    // data members initialized during configuration
    TextInFile* _infile{nullptr};

    // column indices
    int _nextIndex{0};
    int _positionIndex{-1};
    int _sizeIndex{-1};
    int _densityIndex{-1};
    int _massIndex{-1};
    int _metallicityIndex{-1};
    int _temperatureIndex{-1};
    int _velocityIndex{-1};
    int _parametersIndex{-1};
    int _numParameters{0};

    // mass and mass density policy
    double _multiplier{0.};
    double _maxTemperature{0.};
};

////////////////////////////////////////////////////////////////////

#endif
