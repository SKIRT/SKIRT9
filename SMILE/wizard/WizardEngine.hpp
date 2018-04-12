/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WIZARDENGINE_HPP
#define WIZARDENGINE_HPP

#include "Basics.hpp"
#include <QObject>
class Item;
class PropertyHandler;
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** The WizardEngine class implements the back end of the wizard through which users can create and
    edit SMILE datasets. The MainWindow creates a single WizardEngine instance. This WizardEngine
    object manages the dataset under construction, keeps track of the wizard's state (i.e. which
    question should currently be displayed), allowing it to advance and retreat, and generates the
    user interface pane corresponding to the current state. */
class WizardEngine : public QObject
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor places the WizardEngine in its initial state,
        corresponding to the question "what would you like to do". */
    explicit WizardEngine(QObject* parent);

    /** The destructor deallocates the dataset under construction and releases any other resources
        being held. */
    ~WizardEngine();

    // ================== State Handling ====================

public:
    /** This function returns true if the wizard can currently advance; false otherwise. */
    bool canAdvance();

    /** This function returns true if the wizard can currently retreat; false otherwise. */
    bool canRetreat();

    /** This function returns true if the wizard currently holds unsaved information; false
        otherwise. */
    bool isDirty();

    /** This function returns the path of the file to which the current dataset has already been
        saved (although it may have changed since then), or the empty string if it has never been
        saved. */
    QString filepath();

    /** This function returns a pointer to the schema definition for the current dataset (without
        transfer of ownership), or the null pointer if no dataset has yet been opened or created.
        */
    const SchemaDef* schema();

private:
    /** This function returns the property index for the specified child item in the specified
        parent item, or zero if the specified child item is not a child of the specified parent
        item. */
    int propertyIndexForChild(Item* parent, Item* child);

    /** If no argument is given, this function creates and returns a unique pointer to a property
        handler for the current property of the current item. If an argument is given, it serves as
        an offset to the property index, i.e. 1 means the next property and -1 means the previous
        property of the current item. The function does \em not check whether the adjusted index is
        within range. */
    std::unique_ptr<PropertyHandler> createCurrentPropertyHandler(int offset=0);

    /** This function descends the hierarchy as deeply as possible, starting from the current item
        and property index, updating the current item and property index as a result. */
    void descendToDeepest();

public slots:
    /** This function advances the wizard to the next state. It should only be called if
        canAdvance() returns true. */
    void advance();

    /** This function retreats the wizard to the previous state. It should only be called if
        canRetreat() returns true. */
    void retreat();

    /** This function advances the wizard to a state that starts editing the specified item in the
        current item list property. */
    void advanceToEditSubItem(int subItemIndex);

signals:
    /** This signal is emitted when the return value of the canAdvance() function may have changed.
        The argument specifies the new value. */
    void canAdvanceChangedTo(bool canAdvance);

    /** This signal is emitted when the return value of the canRetreat() function may have changed.
        The argument specifies the new value. */
    void canRetreatChangedTo(bool canRetreat);

    /** This signal is emitted when the state of the wizard has changed. */
    void stateChanged();

    /** This signal is emitted when the filename in which the current dataset was saved has
        changed. */
    void titleChanged();

    /** This signal is emitted when the dirty state of the current dataset has changed. */
    void dirtyChanged();

public:
    /** This function emits the stateChanged(), canAdvanceChangedTo() and canRetreatChangedTo()
        signals. */
    void emitStateChanged();

    // ================== State Updating ====================

public slots:
    /** This function updates the basic choice to the specified values. */
    void setBasicChoice(bool openExisting, string libraryPath, string schemaName);

    /** If the current root does not have the specified type (or if there is no current root), this
        function deletes the current dataset (if present), and replaces it by a newly created root
        item of the specified type. If the current root already has the specified type, this
        function does nothing. */
    void setRootType(string newRootType);

    /** This function deletes the current dataset (if present), and replaces it by the new dataset
        specified through it root item. The function adopts ownership for specified dataset. This
        function further clears the dirty flag and remembers the filepath from which the dataset
        was loaded. */
    void hierarchyWasLoaded(Item* root, QString filepath);

    /** This function updates the flag that indicates whether the value of the property currently
        being handled is valid. */
    void setPropertyValid(bool valid);

    /** This function sets the dirty flag. */
    void hierarchyWasChanged();

    /** This function clears the dirty flag and remembers the filepath in which the dataset was
        saved. */
    void hierarchyWasSaved(QString filepath);

    // ================== GUI Generation ====================

public:
    /** This function creates a fresh QWidget object corresponding to the current wizard state,
        returns a pointer to it, and transfers ownership to the caller. The QWidget has no parent,
        but is otherwise fully equipped to handle the keyboard and mouse events for any UI elements
        it contains. For example, the UI elements are equipped so that they can properly update the
        corresponding portion of the dataset under construction. */
    QWidget* createPane();

    /** This function returns a human-readable description of the path to the property currently
        being handled in the dataset. */
    QString hierarchyPath();

    // ================== Data Members ====================

private:
    // ---- data members representing the state of the wizard ----

    // the top-level state; always valid
    enum class State { BasicChoice, CreateRoot, OpenHierarchy, ConstructHierarchy, SaveHierarchy };
    State _state = State::BasicChoice;

    // part of the basic choice: true = open existing dataset, false = create new dataset; always valid
    bool _openExisting = false;

    // part of the basic choice: the (file)name of the schema for the dataset under construction;
    // always valid but remains empty until BasicChoice has been completed at least once
    string _schemaName;

    // the schema for the dataset under construction; pointer has ownership;
    // always valid but remains null until BasicChoice has been completed at least once
    std::unique_ptr<SchemaDef> _schema;

    // the root item of the dataset under construction; pointer has ownership;
    // always valid but remains null until CreateRoot has been completed at least once
    std::unique_ptr<Item> _root;

    // the item in the dataset currently being handled; pointer is a reference without ownership;
    // valid only during ConstructHierarchy
    Item* _current = nullptr;

    // the zero-based index of the property currently being handled (in the current item);
    // valid only during ConstructHierarchy
    int _propertyIndex = -1;

    // the zero-based index of the currently selected sub-item of the current item list property,
    // or -1 when editing the item list property itself;
    // valid only if the current property is an item list
    int _subItemIndex = -1;

    // true if the value of the property being handled is valid, false otherwise
    // valid only during ConstructHierarchy
    bool _propertyValid = false;

    // ---- data members related to the state of the wizard ----

    // true if the current dataset holds unsaved information, false otherwise; always valid
    bool _dirty = false;

    // the path of the file to which the current dataset has already been saved (although it may have changed),
    // or the empty string if it has never been saved; always valid
    QString _filepath;
};

////////////////////////////////////////////////////////////////////

#endif
