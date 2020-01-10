/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WIZARDENGINE_HPP
#define WIZARDENGINE_HPP

#include "Basics.hpp"
#include "NameManager.hpp"
#include <QObject>
#include <stack>
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
    /** This function creates and returns a unique pointer to a property handler for the property
        of the current item with the specified property index. The function does \em not check
        whether the property index is within range. */
    std::unique_ptr<PropertyHandler> createPropertyHandler(int propertyIndex);

    /** This function returns the property index for the specified child item in its parent item,
        or -1 if the specified child item does not have a parent or if it is not a child of its
        parent (which should never happen). The result does not depend on the evaluation of
        conditional expressions against the current name sets. */
    int propertyIndexForChild(Item* child);

public slots:
    /** This function advances the wizard to the next state. It should only be called if
        canAdvance() returns true. Th optional flags are used for recursive calls. If \em state is
        set to false, the current advance state is not saved. If \em descend is set to false, the
        function will not descend into the current compound property. */
    void advance(bool state = true, bool descend = true);

    /** This function advances the wizard to a state that starts editing the specified item in the
        current item list property. */
    void advanceToEditSubItem(int subItemIndex);

    /** This function retreats the wizard to the previous state. It should only be called if
        canRetreat() returns true. */
    void retreat();

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

    /** This function intializes the state of the wizard as if it has just been started, and causes
        the basic choice pane to be shown. */
    void restartWizard();

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
    // ---- data members reflecting the basic choice and dataset root creation ----

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

    // ---- data members determining the information currently being handled ----

    // the wizard stage; always valid
    enum class Stage { BasicChoice, CreateRoot, OpenHierarchy, ConstructHierarchy, SaveHierarchy };
    Stage _stage = Stage::BasicChoice;

    // the item in the dataset currently being handled; pointer is a reference without ownership;
    // valid only during ConstructHierarchy
    Item* _current = nullptr;

    // the zero-based index of the first or only property currently being handled (in the current item);
    // valid only during ConstructHierarchy
    int _firstPropertyIndex = -1;

    // the zero-based index of the last property currently being handled (in the current item);
    // if a single property is being handled, then _lastPropertyIndex ==  _firstPropertyIndex;
    // valid only during ConstructHierarchy
    int _lastPropertyIndex = -1;

    // the zero-based index of the currently selected sub-item of the current item list property,
    // or -1 when editing the item list property itself;
    // valid only during ConstructHierarchy and if the current property is an item list
    int _subItemIndex = -1;

    // ---- data members implementing the advance/retreat state stack ----

    class State
    {
    public:
        State(Stage stage, Item* current, int firstIndex, int lastIndex, int subIndex)
            : _stage(stage), _current(current), _firstIndex(firstIndex), _lastIndex(lastIndex), _subIndex(subIndex)
        {}
        void getState(Stage& stage, Item*& current, int& firstIndex, int& lastIndex, int& subIndex)
        {
            stage = _stage;
            current = _current;
            firstIndex = _firstIndex;
            lastIndex = _lastIndex;
            subIndex = _subIndex;
        }

    private:
        Stage _stage;
        Item* _current;
        int _firstIndex, _lastIndex, _subIndex;
    };

    // stack on which the state is pushed before each advance, and popped after each retreat (with some complications)
    std::stack<State> _stateStack;

    // stack on which indices into the state stack are pushed to prevent retreating into subitem editing
    std::stack<size_t> _stateIndexStack;

    // ---- other data members related to the state of the wizard ----

    // true if the value of the property being handled is valid, false otherwise
    // valid only during ConstructHierarchy
    bool _propertyValid = false;

    // true if the current dataset holds unsaved information, false otherwise; always valid
    bool _dirty = false;

    // the path of the file to which the current dataset has already been saved (although it may have changed),
    // or the empty string if it has never been saved; always valid
    QString _filepath;

    // the name manager used to evaluate the conditional attributes of a property
    NameManager _nameMgr;
};

////////////////////////////////////////////////////////////////////

#endif
