/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SNAPSHOT_HPP
#define SNAPSHOT_HPP

#include "Array.hpp"
#include "Box.hpp"
#include "Direction.hpp"
#include "Position.hpp"
#include "SnapshotParameter.hpp"
class EntityCollection;
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
        the log message issued after the file is successfully opened. */
    void open(const SimulationItem* item, string filename, string description);

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
    /** This function specifies a mapping (defined by the \em columns argument) between the
        "physical" columns in the file being imported (defined by the column information in the
        file header) and the required "logical" columns (defined by calling the other configuration
        functions). For information on the syntax and semantics of the \em columns string contents,
        refer to the description of the TextInFile::useColumns() function.

        This function can be called with a non-empty \em columns string at most once for each
        snapshot, and such invocation should occur \em before the first invocation of any other
        configuration function. Calling this function with an empty \em columns string is
        equivalent to not calling it at all. */
    void useColumns(string columns);

    /** This function configures the snapshot to import a spatial position with three components
        (x,y,z). The default unit is pc. */
    void importPosition();

    /** This function configures the snapshot to import a spatial radial size. The default unit is
        pc. */
    void importSize();

    /** This function configures the snapshot to import a cuboid lined up with the coordinate axes,
        defined by six components (xmin,ymin,zmin,xmax,ymax,zmax). The default unit is pc. */
    void importBox();

    /** This function configures the snapshot to import a mass density per unit of volume (for use
        by media). The default unit is Msun/pc3. It is allowed to combine the importMassDensity()
        and importMass() options, supporting special use cases where the volume of the entity
        cannot be derived otherwise. However, combining the "mass" family functions with the
        "number" family functions is prohibited and leads to undefined behavior. */
    void importMassDensity();

    /** This function configures the snapshot to import a mass, i.e. mass density integrated over
        volume (for use by media). The default unit is Msun. It is allowed to combine the
        importMassDensity() and importMass() options, supporting special use cases where the volume
        of the entity cannot be derived otherwise. However, combining the "mass" family functions
        with the "number" family functions is prohibited and leads to undefined behavior. */
    void importMass();

    /** This function configures the snapshot to import a number density per unit of volume (for
        use by media). The default unit is 1/cm3. It is allowed to combine the
        importNumberDensity() and importNumber() options, supporting special use cases where the
        volume of the entity cannot be derived otherwise. However, combining the "mass" family
        functions with the "number" family functions is prohibited and leads to undefined behavior.
        */
    void importNumberDensity();

    /** This function configures the snapshot to import a number, i.e. number density integrated
        over volume (for use by media). The default unit is 1. It is allowed to combine the
        importNumberDensity() and importNumber() options, supporting special use cases where the
        volume of the entity cannot be derived otherwise. However, combining the "mass" family
        functions with the "number" family functions is prohibited and leads to undefined behavior.
        */
    void importNumber();

    /** This function configures the snapshot to import the current mass (for use by sources). Not
        to be confused with the initial mass, which is often requested by single stellar population
        %SED families. The default unit is Msun. */
    void importCurrentMass();

    /** This function configures the snapshot to import a (dimensionless) metallicity fraction. */
    void importMetallicity();

    /** This function configures the snapshot to import a temperature. The default unit is K. */
    void importTemperature();

    /** This function configures the snapshot to import a velocity with three components (x,y,z).
        The default unit is km/s. */
    void importVelocity();

    /** This function configures the snapshot to import a single velocity dispersion value,
        specifying a random offset to the bulk velocity with a spherically symmetric Gaussian
        distribution. The default unit is km/s. */
    void importVelocityDispersion();

    /** This function configures the snapshot to import a magnetic field vector with three
        components (x,y,z). The default unit is \f$\mu \mathrm{G}\f$. */
    void importMagneticField();

    /** This function configures the snapshot to import a sequence of parameters as described by
        the specified list of SnapshotParameter metadata objects. */
    void importParameters(const vector<SnapshotParameter>& parameters);

    /** This function sets the policy for calculating the mass or mass density at a given spatial
        position in the snapshot. If the policy is not set during configuration, or if neither of
        the importMass() and importDensity() functions was invoked, attempting to obtain mass or
        mass density values results in undefined behavior.

        If the policy is set during configuration, the argument values determine the heuristic for
        calculating the mass or density for each entity:

        - If the \em maxTemperature argument is positive, and the temperature field is being
        imported and its value for a given entity is positive as well, and the imported temperature
        is above the maximum temperature, then the mass or density corresponding to this entity is
        considered to be zero.

        - Otherwise, the intrinsic mass or mass density is obtained from the mass or density fields
        using a mechanism that depends on the snapshot type and possibly on the configuration. If
        the \em useMetallicity argument is true, and the metallicity field is being imported, this
        intrinsic value is multiplied by the metallicity fraction. Finally, it is multiplied by the
        value of the \em multiplier argument to obtain the final result. */
    void setMassDensityPolicy(double multiplier, double maxTemperature, bool useMetallicity);

    /** This function notifies the snapshot that one of the getEntities() functions may be called
        during the simulation, implying that the snapshot must prebuild the required search data
        structures. */
    void setNeedGetEntities();

protected:
    /** This function returns the column index of the first position field, or -1 if this is not
        being imported, for use by subclasses. */
    int positionIndex() const { return _positionIndex; }

    /** This function returns the column index of the size field, or -1 if this is not being
        imported, for use by subclasses. */
    int sizeIndex() const { return _sizeIndex; }

    /** This function returns the column index of the first box field, or -1 if this is not
        being imported, for use by subclasses. */
    int boxIndex() const { return _boxIndex; }

    /** This function returns the column index of the density field, or -1 if this is not being
        imported, for use by subclasses. */
    int densityIndex() const { return _densityIndex; }

    /** This function returns the column index of the mass field, or -1 if this is not being
        imported, for use by subclasses. */
    int massIndex() const { return _massIndex; }

    /** This function returns the column index of the initial mass field, or -1 if this is not
        being imported, for use by subclasses. */
    int initialMassIndex() const { return _initialMassIndex; }

    /** This function returns the column index of the current mass field, or -1 if this is not
        being imported, for use by subclasses. */
    int currentMassIndex() const { return _currentMassIndex; }

    /** This function returns the column index of the metallicity field, or -1 if this is not being
        imported, for use by subclasses. */
    int metallicityIndex() const { return _metallicityIndex; }

    /** This function returns the column index of the age field, or -1 if this is not being
        imported, for use by subclasses. */
    int ageIndex() const { return _ageIndex; }

    /** This function returns the column index of the temperature field, or -1 if this is not being
        imported, for use by subclasses. */
    int temperatureIndex() const { return _temperatureIndex; }

    /** This function returns the column index of the first velocity field, or -1 if this is not
        being imported, for use by subclasses. */
    int velocityIndex() const { return _velocityIndex; }

    /** This function returns the column index of the velocity dispersion field, or -1 if this is
        not being imported, for use by subclasses. */
    int velocityDispersionIndex() const { return _velocityDispersionIndex; }

    /** This function returns the column index of the first magnetic field field, or -1 if this is
        not being imported, for use by subclasses. */
    int magneticFieldIndex() const { return _magneticFieldIndex; }

    /** This function returns the column index of the first field in the parameter list, or -1 if
        this is not being imported, for use by subclasses. */
    int parametersIndex() const { return _parametersIndex; }

    /** This function returns the number of parameters being imported (which may be zero). */
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
        request to use the metallicity, and the metallicity field is being imported. Returns false
        otherwise. For use by subclasses. */
    bool useMetallicity() const { return _hasDensityPolicy && _useMetallicity && _metallicityIndex >= 0; }

    /** This function returns true if the user configured the mass or mass density policy with a
        nonzero maximum temperature and the temperature field is being imported. Returns false
        otherwise. For use by subclasses. */
    bool useTemperatureCutoff() const { return _hasDensityPolicy && _maxTemperature > 0 && _temperatureIndex >= 0; }

    /** This function returns true if one of the getEntities() functions may be called during the
        simulation, implying that the snapshot must prebuild the required search data structures.
        Returns false otherwise. For use by subclasses. */
    bool needGetEntities() const { return _needGetEntities; }

    /** This function issues log messages with statistics on the imported masses. It is implemented
        here for use by subclasses. */
    void logMassStatistics(int numIgnored, double totalOriginalMass, double totalMetallicMass,
                           double totalEffectiveMass);

    //============== Interrogation (to be implemented in subclass) =============

public:
    /** This function returns the extent of the complete spatial domain of the snapshot as a box
        lined up with the coordinate axes. */
    virtual Box extent() const = 0;

    /** This function returns the number of entities \f$N_\mathrm{ent}\f$ in the snapshot. */
    virtual int numEntities() const = 0;

    /** This function returns the volume of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the index is out of range, if no density policy has been set, or no
        mass/density information is being imported, the behavior is undefined. */
    virtual double volume(int m) const = 0;

    /** This function returns the mass or number density for the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the index is out of range, if no density policy has been set, or no
        mass/density information is being imported, the behavior is undefined. */
    virtual double density(int m) const = 0;

    /** This function returns the mass or number density represented by the snapshot at a given
        point \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero. If no
        density policy has been set or no mass/density information is being imported, the behavior
        is undefined. */
    virtual double density(Position bfr) const = 0;

    /** This function returns the total mass or number represented by the snapshot, which is
        equivalent to the mass or number density integrated over the complete spatial domain. If no
        density policy has been set or no mass/density information is being imported, the behavior
        is undefined. */
    virtual double mass() const = 0;

    /** This function returns a characteristic position for the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. Such a position is always available for all snapshot types, regardless
        of whether a position is explicitly being imported or not. If the index is out of range,
        the behavior is undefined. */
    virtual Position position(int m) const = 0;

    /** This function returns a random position within the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$, drawn from an appropriate probability distribution depending on the
        snapshot type (e.g. uniform for cells, and some smoothing kernel for particles). If the
        index is out of range, the behavior is undefined. */
    virtual Position generatePosition(int m) const = 0;

    /** This function returns a random position within the spatial domain of the snapshot, drawn
        from the mass density distribution represented by the snapshot. If no density policy has
        been set or no mass/density information is being imported, the behavior is undefined. */
    virtual Position generatePosition() const = 0;

protected:
    /** This function returns a reference to an array containing the imported properties (in column
        order) for the entity with index \f$0\le m \le N_\mathrm{ent}-1\f$. If the index is out of
        range, the behavior is undefined. */
    virtual const Array& properties(int m) const = 0;

    /** This function returns the index \f$0\le m \le N_\mathrm{ent}-1\f$ for the entity at or
        nearest to the specified point \f${\bf{r}}\f$, or -1 if the point is outside the domain or
        if there are no entities in the snapshot. For a cell-based snapshot, the function returns
        the index of the cell containing the given point. For a particle-based snapshot, the
        function returns the index of the particle whose center is nearest to the given point. */
    virtual int nearestEntity(Position bfr) const = 0;

public:
    /** This function replaces the contents of the specified entity collection by the set of
        entities that overlap the specified point \f${\bf{r}}\f$, with their corresponding weights.
        If the point is outside the domain or otherwise does not overlap any entity, the collection
        will be empty.

        For a cell-based snapshot, the function returns the cell containing the given point, if
        any. The weight is set to 1. For a particle-based snapshot, the function returns all
        particles with a smoothing kernel that overlaps the given point. The weight of a particle
        is given by the particle's smoothing kernel value at the given point. */
    virtual void getEntities(EntityCollection& entities, Position bfr) const = 0;

    /** This function replaces the contents of the specified entity collection by the set of
        entities that overlap the specified path with starting point \f${\bf{r}}\f$ and direction
        \f${\bf{k}}\f$, with their corresponding weights. If the path does not overlap any entity,
        the collection will be empty.

        For a cell-based snapshot, the weight of a cell is given by the length of the path segment
        inside the cell. For a particle-based snapshot, the weight of a particle is given by the
        effective length seen by the path as it crosses the particle's smoothing kernel. */
    virtual void getEntities(EntityCollection& entities, Position bfr, Direction bfk) const = 0;

    //============== Interrogation implemented here =============

public:
    /** This function returns true if the snapshot holds number (density) values, and false if it
        holds mass (density) values (or if no mass or density column is being imported). */
    bool holdsNumber() const { return _holdsNumber; }

    /** This function returns the volume of the complete domain of the snapshot, taken to be a box
        lined up with the coordinate axes. */
    double volume() const;

    /** This function returns true if the initial mass is being imported, and false otherwise. */
    bool hasInitialMass() const { return _initialMassIndex >= 0; }

    /** This function returns the initial mass of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the initial mass is not being imported, or the index is out of
        range, the behavior is undefined. */
    double initialMass(int m) const;

    /** This function returns the initial mass of the entity nearest to (or at) the specified point
        \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero. If the
        initial mass is not being imported, the behavior is undefined. */
    double initialMass(Position bfr) const;

    /** This function returns true if the current mass is being imported, and false otherwise. */
    bool hasCurrentMass() const { return _currentMassIndex >= 0; }

    /** This function returns the current mass of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the current mass is not being imported, or the index is out of
        range, the behavior is undefined. */
    double currentMass(int m) const;

    /** This function returns the current mass of the entity nearest to (or at) the specified point
        \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero. If the
        current mass is not being imported, the behavior is undefined. */
    double currentMass(Position bfr) const;

    /** This function returns true if the metallicity is being imported, and false otherwise. */
    bool hasMetallicity() const { return _metallicityIndex >= 0; }

    /** This function returns the metallicity of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the metallicity is not being imported, or the index is out of
        range, the behavior is undefined. */
    double metallicity(int m) const;

    /** This function returns the metallicity of the entity nearest to (or at) the specified point
        \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero. If the
        metallicity is not being imported, the behavior is undefined. */
    double metallicity(Position bfr) const;

    /** This function returns true if the age is being imported, and false otherwise. */
    bool hasAge() const { return _ageIndex >= 0; }

    /** This function returns the age of the entity with index \f$0\le m \le N_\mathrm{ent}-1\f$.
        If the age is not being imported, or the index is out of range, the behavior is undefined.
        */
    double age(int m) const;

    /** This function returns the age of the entity nearest to (or at) the specified point
        \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero. If the age
        is not being imported, the behavior is undefined. */
    double age(Position bfr) const;

    /** This function returns true if the temperature is being imported, and false otherwise. */
    bool hasTemperature() const { return _temperatureIndex >= 0; }

    /** This function returns the temperature of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the temperature is not being imported, or the index is out of
        range, the behavior is undefined. */
    double temperature(int m) const;

    /** This function returns the temperature of the entity nearest to (or at) the specified point
        \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero. If the
        temperature is not being imported, the behavior is undefined. */
    double temperature(Position bfr) const;

    /** This function returns true if the velocity is being imported, and false otherwise. */
    bool hasVelocity() const { return _velocityIndex >= 0; }

    /** This function returns the velocity of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the velocity is not being imported, or the index is out of range,
        the behavior is undefined. */
    Vec velocity(int m) const;

    /** This function returns the velocity of the entity nearest to (or at) the specified point
        \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero velocity. If
        the velocity is not being imported, the behavior is undefined. */
    Vec velocity(Position bfr) const;

    /** This function returns true if the velocity dispersion is being imported, and false
        otherwise. */
    bool hasVelocityDispersion() const { return _velocityDispersionIndex >= 0; }

    /** This function returns the velocity dispersion of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the velocity dispersion is not being imported, or the index is out
        of range, the behavior is undefined. */
    double velocityDispersion(int m) const;

    /** This function returns the velocity dispersion of the entity nearest to (or at) the
        specified point \f${\bf{r}}\f$. If the point is outside the domain, the function returns
        zero dispersion. If the velocity dispersion is not being imported, the behavior is
        undefined. */
    double velocityDispersion(Position bfr) const;

    /** This function returns true if the magnetic field is being imported, and false otherwise. */
    bool hasMagneticField() const { return _magneticFieldIndex >= 0; }

    /** This function returns the magnetic field vector of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the magnetic field is not being imported, or the index is out of
        range, the behavior is undefined. */
    Vec magneticField(int m) const;

    /** This function returns the magnetic field vector of the entity nearest to (or at) the
        specified point \f${\bf{r}}\f$. If the point is outside the domain, the function returns a
        zero magnetic field. If the magnetic field is not being imported, the behavior is
        undefined. */
    Vec magneticField(Position bfr) const;

    /** This function returns true if parameters are being imported (i.e. if the number of imported
         parameters is nonzero), and false otherwise. */
    bool hasParameters() const { return _numParameters > 0; }

    /** This function stores the parameters of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$ into the given array. If parameters are not being imported, or the
        index is out of range, the behavior is undefined. */
    void parameters(int m, Array& params) const;

    /** This function stores the parameters of the entity nearest to (or at) the specified point
        \f${\bf{r}}\f$ into the given array. If the point is outside the domain, the function
        returns the appropriate number of zero parameter values. If parameters are not being
        imported, the behavior is undefined. */
    void parameters(Position bfr, Array& params) const;

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
    int _boxIndex{-1};
    int _densityIndex{-1};
    int _massIndex{-1};
    int _initialMassIndex{-1};
    int _currentMassIndex{-1};
    int _metallicityIndex{-1};
    int _ageIndex{-1};
    int _temperatureIndex{-1};
    int _velocityIndex{-1};
    int _velocityDispersionIndex{-1};
    int _magneticFieldIndex{-1};
    int _parametersIndex{-1};
    int _numParameters{0};

    // mass and mass density policy
    double _multiplier{0.};
    double _maxTemperature{0.};
    bool _useMetallicity{false};
    bool _hasDensityPolicy{false};
    bool _holdsNumber{false};  // true if snapshot holds number (density); false if it holds mass (density)
    bool _needGetEntities{false};
};

////////////////////////////////////////////////////////////////////

#endif
