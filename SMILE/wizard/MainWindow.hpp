/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include "Basics.hpp"
#include <QMainWindow>
class WizardEngine;
class QLabel;
class QLayout;

////////////////////////////////////////////////////////////////////

/** A single instance of the MainWindow class is created and shown during application startup in
    the main() function. As the name implies, the class represents the application's main (and
    only) window. It provides the application menu and a status bar, and it hosts the wizard that
    guides the user through the creation and editing of SMILE datasets. */
class MainWindow : public QMainWindow
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates the GUI for the main window, including the
        application menu, a status bar, and the areas hosting the wizard. It also creates
        appropriate actions and connections to support the top-level events. */
    MainWindow();

private:
    /** This function positions and resizes the main window according to the settings saved during
        the previous session. */
    void readSettings();

    /** This function saves the current position and size of the main window so that the
        information can be used to reproduce the same configuration in a future session. */
    void writeSettings();

    // ==================== Event Handling ====================

public slots:
    /** This function replaced the wizard pane displayed in the central area by a pane newly
        retrieved from the wizard engine. */
    void replaceWizardPane();

    /** This function updates the window title to reflect the current filename, which can be
        retrieved from the wizard engine. */
    void updateTitle();

    /** This function updates the window title bar to reflect the current dirty state, which can be
        retrieved from the wizard engine. */
    void updateDirtyState();

    /** This function assumes that it is called from a signal, and that the sender QObject has a
        property named "URL" of type QString. The function retrieves the value of that property and
        then opens the specified URL in the user's default browser. */
    void browseUrl();

protected:
    /** This function is invoked for key presses that aren't handled in a sub-pane. It handles the
        keyboard shortcuts for the advance and retreat actions. */
    void keyPressEvent(QKeyEvent* event) override;

    /** This function is invoked when the user ctrl-clicks or right-clicks in the window. It
        creates and shows a context menu with help options opening a relevant URL in the default
        browser. */
    void contextMenuEvent(QContextMenuEvent* event) override;

    /** This function is invoked when the user attempts to close the main window or to quit the
        application. If there are any unsaved changes, the function offers the user a chance to
        cancel the close or quit operation. Otherwise, or if the user decides to quit anyway, the
        function stores the main window's position and size for future reference and allows the
        application to exit. */
    void closeEvent(QCloseEvent* event) override;

    // ==================== Data Members ====================

private:
    WizardEngine* _wizard;
    QWidget* _wizardPane;
    QLayout* _wizardLayout;
    QLabel* _pathLabel;
    bool _acceptedCloseEvent = false;  // becomes true after the user accepted a window close event
                                       // to avoid the question being asked twice when quitting
};

////////////////////////////////////////////////////////////////////

#endif
