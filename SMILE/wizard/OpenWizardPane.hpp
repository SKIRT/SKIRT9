/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OPENWIZARDPANE_HPP
#define OPENWIZARDPANE_HPP

#include <Basics.hpp>
#include <QWidget>
class Item;
class SchemaDef;
class QLabel;
class QPushButton;

////////////////////////////////////////////////////////////////////

/** An OpenWizardPane instance displays the user interface for loading a new SMILE dataset from a
    SMILE XML file. */
class OpenWizardPane : public QWidget
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor creates and initializes the GUI for this pane. The first
        argument specifies the schema definition for the dataset under construction. The second
        argument provides a filepath that has been previously used during this session, if any.
        This file path is used to position the open dialog in the file system. The third argument
        provides the dirty state of the currently open dataset, which will be overridden. This
        information is used to propertly warn the user. The last argument specifies the object that
        will be notified of a successful open and load operation through invocation of the object's
        relevant slot. */
    explicit OpenWizardPane(const SchemaDef* schema, QString filepath, bool dirty, QObject* target);

    // ==================== Event Handling ====================

public slots:
    /** This function obtains an open file path from the user, and then loads the SMILE dataset
        from that file. After a successful load operation, the function notifies the target object
        by emitting a hierarchyWasLoaded() signal. */
    void open();

signals:
    /** This signal is emitted after the dataset has been successfully loaded. */
    void hierarchyWasLoaded(Item* root, QString filepath);

    // ==================== Data members ======================

private:
    const SchemaDef* _schema;
    QString _filepath;
    bool _dirty;
    QLabel* _filepathLabel;
    QPushButton* _openButton;
};

////////////////////////////////////////////////////////////////////

#endif
