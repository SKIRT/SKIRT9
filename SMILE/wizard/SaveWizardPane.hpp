/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SAVEWIZARDPANE_HPP
#define SAVEWIZARDPANE_HPP

#include <Basics.hpp>
#include <QWidget>
class Item;
class SchemaDef;
class QLabel;
class QPushButton;

////////////////////////////////////////////////////////////////////

/** A SaveWizardPane instance displays the user interface for saving the current SMILE dataset
     into a SMILE XML file. */
class SaveWizardPane : public QWidget
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. The first
        argument specifies the schema definition for the dataset to be saved. The second argument
        provides a pointer to the root of the dataset item hierarchy. The subsequent arguments
        provide the filepath in which the simulation item hierarchy has been previously saved, if
        any, and the current dirty state of the hierarchy. The last argument specifies the object
        that will be notified of a successful save operation through invocation of the object's
        relevant slot. */
    explicit SaveWizardPane(const SchemaDef* schema, Item* root, QString filepath, bool dirty, QObject* target);

    // ==================== Event Handling ====================

public slots:
    /** If the dataset was previously saved to a known filepath, this function
        saves the dataset again to the same file path, replacing the previous
        file, and notifies the target object by emitting a hierarchyWasSaved() signal. If no
        previous filepath is known, the function does nothing. */
    void save();

    /** This function displays the appropriate dialog to obtain a file path from the user, saves
        the dataset to the selected file path, and notifies the target object by
        emitting a hierarchyWasSaved() signal. */
    void saveAs();

    /** This function attempts to restart the wizard. */
    void restart();

    /** This function attempts to quit the application. */
    void quit();

private:
    /** This private function saves the dataset to the specfied file path, and
        notifies the target object by emitting a hierarchyWasSaved() signal. */
    void saveToFile(QString filepath);

    /** This function enables or disables the Save push button depending on the filename and dirty
        state, and puts the current filepath into the corresponding label. */
    void updateSaveInfo();

signals:
    /** This signal is emitted after the simulation item hierarchy has been successfully saved. */
    void hierarchyWasSaved(QString filepath);

    /** This signal is emitted after the user pushes the restart button (and declines to save any
        changes, if applicable). */
    void restartWizard();

    // ==================== Data members ======================

private:
    const SchemaDef* _schema;
    Item* _root;
    QString _filepath;
    bool _dirty;
    QLabel* _filepathLabel;
    QPushButton* _saveButton;
    QPushButton* _saveAsButton;
    QPushButton* _restartButton;
    QPushButton* _quitButton;
};

////////////////////////////////////////////////////////////////////

#endif
