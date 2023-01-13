/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BASICCHOICEWIZARDPANE_HPP
#define BASICCHOICEWIZARDPANE_HPP

#include "Basics.hpp"
#include <QWidget>
class QAbstractButton;
class QButtonGroup;
class QLabel;
class QPushButton;
class QRadioButton;
class QLayout;

////////////////////////////////////////////////////////////////////

/** A BasicChoiceWizardPane instance displays the user interface for the question "what do you want
    to do" at the start of the wizard guidance process, and it passes the selection on to the
    target object. The pane also offers a button to select the path to the SMILE schema library,
    and a button to revert to the default built-in schema library path. The currently selected
    schema library path is persistent across invocations of the application. */
class BasicChoiceWizardPane : public QWidget
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. The first
        two arguments specify the choice that will be selected when the pane is initially
        displayed. If the schema name does not match any of the possible choices, no choices will
        be selected. The third argument indicates the dirty state of the current simulation item
        hierarchy. The last argument specifies the object that will be notified of changes in the
        selection through invocation of the object's setBasicChoice() slot. */
    explicit BasicChoiceWizardPane(bool initialOpenExisting, string initialSchemaName, bool dirty, QObject* target);

    // ==================== Event Handling ====================

public slots:
    /** This function updates the basic choice in reaction to a user click on one of the radio
        buttons. If the current hierarchy is dirty, and the new choice differs from the previous
        one, the function asks for confirmation from the user before actually updating the choice.
        If the new choice is the same as the previous one, the function does nothing. */
    void setBasicChoice(QAbstractButton* button);

    /** This function allows the user to select a new directory containing the SMILE schema library
        in reaction to a user click on the corresponding push button. The function stores the new
        library path in the persistent application settings, and then then calls the
        loadSchemaInfo() and updateChoiceInterface() functions to update the internal schema
        information and the user interface accordingly. */
    void selectLibrary();

    /** This function sets the SMILE schema library path to the built-in path, assuming it exists,
        in reaction to a user click on the corresponding push button. The function stores the new
        library path in the persistent application settings, and then calls the loadSchemaInfo()
        and updateChoiceInterface() functions to update the internal schema information and the
        user interface accordingly. */
    void setDefaultLibrary();

private:
    /** This function determines the absolute path to the built-in schema library, if it is found, or
        the empty string otherwise. The function searches relative to the application's executable
        for the directories used by the build system on various host systems. */
    void findDefaultLibraryPath();

    /** This function loads the current schema library path from persistent application
        preferences. If the application preference does not yet exist (perhaps because the
        application is being run for the first time), current schema library path is set to the
        empty string. */
    void loadPersistentLibraryPath();

    /** This function stores the current schema library path to persistent application preferences.
        */
    void storePersistentLibraryPath();

    /** This function retrieves information about the current schema library, including the names
        and user-friendly descriptions of all schemas in the library (or an error message if
        something went wrong). */
    void loadSchemaInfo();

    /** This function updates the user interface after a new schema library info was loaded. It
        updates the libary path text field, and creates new basic choice radio buttons, inserting
        them into the basic choice layout and button group. */
    void updateChoiceInterface();

signals:
    /** This signal is emitted after the basic choice was changed. */
    void basicChoiceWasChanged(bool openExisting, string libraryPath, string schemaName);

    // ==================== Data members ======================

private:
    // copy of constructor arguments
    bool _openExisting;  // open existing dataset or create new one
    string _schemaName;  // name of schema file
    bool _dirty;         // current dataset needs saving

    // schema library info
    string _defaultLibraryPath;    // directory path to the default schema library
    string _libraryPath;           // directory path to the schema library
    string _error;                 // error message: non-empty if an error occurred while loading schema library info
    vector<string> _schemaNames;   // names of the schema files in the library
    vector<string> _schemaTitles;  // corresponding user-friendly descriptions

    // GUI widgets
    QLayout* _choiceLayout;         // the layout containing the basic choice radio buttons
    QButtonGroup* _buttonGroup;     // the button group containing the basic choice radio buttons
    QList<QLabel*> _labels;         // the basic choice error labels currently in the layout, if any
    QList<QRadioButton*> _buttons;  // the basic choice radio buttons currently in the layout, if any
    QLabel* _libraryLabel;          // the text field to display the current schema library path
};

////////////////////////////////////////////////////////////////////

#endif
