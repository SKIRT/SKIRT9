/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SNAPSHOT_HPP
#define SNAPSHOT_HPP

#include "Position.hpp"
#include "Array.hpp"
#include "Box.hpp"
#include "SnapshotParameter.hpp"
class Log;
class Random;
class SimulationItem;
class TextInFile;
class Units;

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

    //========== Input file and context ==========

public:
    /** This function creates an input file object corresponding to the specified file and opens it
        for reading; if the file can't be opened, a FatalError is thrown. It must be called \em
        before invoking any of the configuration functions(). This function takes several
        arguments: (1) \em item specifies a simulation item in the hierarchy of the caller (usually
        the caller itself) used to retrieve context such as an appropriate logger; (2) \em
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
    /** This function retrieves the relevant simulation hierarchy context (such as an appropriate
        logger) from the specified simulation item (usually the caller itself). It is intended for
        use in subclasses that support use cases where the open() function is never invoked. */
    void setContext(const SimulationItem* item);

    /** This function returns a pointer to the input file object. It is intended for use in
        subclasses. */
    TextInFile* infile() { return _infile; }

    /** This function returns a pointer to an appropriate log object. It is intended for use in
        subclasses. */
    Log* log() const { return _log; }

    /** This function returns a pointer to an appropriate units object. It is intended for use in
        subclasses. */
    Units* units() const { return _units; }

    /** This function returns a pointer to an appropriate random generator. It is intended for use
        in subclasses. */
    Random* random() const { return _random; }

    //========== Configuration ==========

public:
    /** This function configures the snapshot to import a spatial position with three components
        (x,y,z). The default unit is pc. */
    void importPosition();

    /** This function configures the snapshot to import a spatial radial size. The default unit is
        pc. */
    void importSize();

    /** This function configures the snapshot to import a mass density per unit of volume. The
        default unit is Msun/pc3. The importMassDensity(), importMass(), importNumberDensity(), and
        importNumber() options are mutually exclusive; calling more than one of these functions for
        the same snapshot results in undefined behavior. */
    void importMassDensity();

    /** This function configures the snapshot to import a mass (i.e. mass density integrated over
        volume). The default unit is Msun. The importMassDensity(), importMass(),
        importNumberDensity(), and importNumber() options are mutually exclusive; calling more than
        one of these functions for the same snapshot results in undefined behavior. */
    void importMass();

    /** This function configures the snapshot to import a number density per unit of volume. The
        default unit is 1/cm3. The importMassDensity(), importMass(), importNumberDensity(), and
        importNumber() options are mutually exclusive; calling more than one of these functions for
        the same snapshot results in undefined behavior. */
    void importNumberDensity();

    /** This function configures the snapshot to import a number (i.e. number density integrated
        over volume). The default unit is 1. The importMassDensity(), importMass(),
        importNumberDensity(), and importNumber() options are mutually exclusive; calling more than
        one of these functions for the same snapshot results in undefined behavior. */
    void importNumber();

    /** This function configures the snapshot to import a (dimensionless) metallicity fraction. */
    void importMetallicity();

    /** This function configures the snapshot to import a temperature. The default unit is K. */
    void importTemperature();

    /** This function configures the snapshot to import a velocity with three components (x,y,z).
        The default unit is km/s. */
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

    /** This function returns true if the snapshot holds number (density) values, and false if it
        holds mass (density) values (or if no mass or density column is being imported). */
    bool holdsNumber() const { return _holdsNumber; }

protected:
    /** This function returns the column index of the first position field, or -1 if this is not
        being imported, for use by subclasses. */
    int positionIndex() const { return _positionIndex; }

    /** This function returns the column index of the size field, or -1 if this is not being
        imported, for use by subclasses. */
    int sizeIndex() const { return _sizeIndex; }

    /** This function returns the column index of the density field, or -1 if this is not being
        imported, for use by subclasses. */
    int densityIndex() const { return _densityIndex; }

    /** This function returns the column index of the mass field, or -1 if this is not being
        imported, for use by subclasses. */
    int massIndex() const { return _massIndex; }

    /** This function returns the column index of the metallicity field, or -1 if this is not being
        imported, for use by subclasses. */
    int metallicityIndex() const { return _metallicityIndex; }

    /** This function returns the column index of the temperature field, or -1 if this is not being
        imported, for use by subclasses. */
    int temperatureIndex() const { return _temperatureIndex; }

    /** This function returns the column index of the first velocity field, or -1 if this is not
        being imported, for use by subclasses. */
    int velocityIndex() const { return _velocityIndex; }

    /** This function returns the column index of the first field in the parameter list, or -1 if
        this is not being imported, for use by subclasses. */
    int parametersIndex() const { return _parametersIndex; }

    /** This function returns the number of parameters being imported (which may be zero), for use
        by subclasses. */
    int numParameters() const { return _numParameters; }

    /** This function returns the mass or mass density multiplier configured by the user, or zero
        if the user did not configure the mass or mass density policy, for use by subclasses. */
    double multiplier() const { return _multiplier; }

    /** This function returns the maximum temperature configured by the user for an entity to have
        mass, or zero if the user did not configure the mass or mass density policy, for use by
        subclasses. */
    double maxTemperature() const { return _maxTemperature; }

    /** This function returns true if the user configured the mass or mass density policy and false
        otherwise, for use by subclasses. */
    bool hasMassDensityPolicy() const { return _hasDensityPolicy; }

    /** This function returns true if the user configured the mass or mass density policy with a
        nonzero temperature and the temperature field is being imported. Returns false otherwise.
        For use by subclasses. */
    bool hasTemperatureCutoff() const { return _hasDensityPolicy && _maxTemperature>0 && _temperatureIndex>=0; }

    //============== Interrogation (to be implemented in subclass) =============

public:
    /** This function returns the extent of the complete spatial domain of the snapshot as a box
        lined up with the coordinate axes. */
    virtual Box extent() const = 0;

    /** This function returns the number of entities \f$N_\mathrm{ent}\f$ in the snapshot. */
    virtual int numEntities() const = 0;

    /** This function returns a characteristic position for the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. Such a position is always available for all snapshot types, regardless
        of whether a position is explicitly being imported or not. If the index is out of range,
        the behavior is undefined. */
    virtual Position position(int m) const = 0;

    /** This function returns the velocity of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the velocity is not being imported, or the index is out of range,
        the behavior is undefined. */
    virtual Vec velocity(int m) const = 0;

    /** This function returns the velocity of the entity nearest to (or at) the specified point
        \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero velocity. If
        the velocity is not being imported, the behavior is undefined. */
    virtual Vec velocity(Position bfr) const = 0;

    /** This function stores the parameters of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$ into the given array. If parameters are not being imported, or the
        index is out of range, the behavior is undefined. */
    virtual void parameters(int m, Array& params) const = 0;

    /** This function returns the mass density represented by the snapshot at a given point
        \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero. If no
        density policy has been set or no mass/density information is being imported, the behavior
        is undefined. */
    virtual double density(Position bfr) const = 0;

    /** This function returns the total mass represented by the snapshot, which is equivalent to
        the mass density integrated over the complete spatial domain. If no density policy has been
        set or no mass/density information is being imported, the behavior is undefined. */
    virtual double mass() const = 0;

    /** This function returns a random position within the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$, drawn from an appropriate probability distribution depending on the
        snapshot type (e.g. uniform for cells, and some smoothing kernel for particles). If the
        index is out of range, the behavior is undefined. */
    virtual Position generatePosition(int m) const = 0;

    /** This function returns a random position within the spatial domain of the snapshot, drawn
        from the mass density distribution represented by the snapshot. If no density policy has
        been set or no mass/density information is being imported, the behavior is undefined. */
    virtual Position generatePosition() const = 0;

    //============== Interrogation implemented here =============

    /** This function returns the volume of the complete domain of the snapshot, taken to be a box
        lined up with the coordinate axes. */
    double volume() const;

    /** This function returns the X-axis surface density of the density distribution represented by
        the snapshot, defined as the integration of the density along the entire X-axis, \f[
        \Sigma_X = \int_{x_\text{min}}^{x_\text{max}} \rho(x,0,0)\, {\text{d}}x.\f] This integral
        is calculated numerically using 10000 samples along the X-axis. If no density policy has
        been set or no mass/density information is being imported, the behavior is undefined. */
    double SigmaX() const;

    /** This function returns the Y-axis surface density of the density distribution represented by
        the snapshot, defined as the integration of the density along the entire Y-axis, \f[
        \Sigma_Y = \int_{y_\text{min}}^{y_\text{max}} \rho(0,y,0)\, {\text{d}}y.\f] This integral
        is calculated numerically using 10000 samples along the Y-axis. If no density policy has
        been set or no mass/density information is being imported, the behavior is undefined. */
    double SigmaY() const;

    /** This function returns the Z-axis surface density of the density distribution represented by
        the snapshot, defined as the integration of the density along the entire Z-axis, \f[
        \Sigma_Z = \int_{z_\text{min}}^{z_\text{max}} \rho(0,0,z)\, {\text{d}}z.\f] This integral
        is calculated numerically using 10000 samples along the Z-axis. If no density policy has
        been set or no mass/density information is being imported, the behavior is undefined. */
    double SigmaZ() const;

    //======================== Data Members ========================

private:
    // data members initialized during configuration
    TextInFile* _infile{nullptr};
    Log* _log{nullptr};
    Units* _units{nullptr};
    Random* _random{nullptr};

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
    bool _hasDensityPolicy{false};
    bool _holdsNumber{false};   // true if snapshot holds number (density); false if it holds mass (density)
};

////////////////////////////////////////////////////////////////////

#endif
