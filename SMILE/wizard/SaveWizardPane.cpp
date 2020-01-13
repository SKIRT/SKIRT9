/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SaveWizardPane.hpp"
#include "Item.hpp"
#include "SchemaDef.hpp"
#include "XmlHierarchyWriter.hpp"
#include <QApplication>
#include <QFileDialog>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <QStandardPaths>
#include <QVBoxLayout>

////////////////////////////////////////////////////////////////////

SaveWizardPane::SaveWizardPane(const SchemaDef* schema, Item* root, QString filepath, bool dirty, QObject* target)
    : _schema(schema), _root(root), _filepath(filepath), _dirty(dirty)

{
    // connect ourselves to the target
    connect(this, SIGNAL(hierarchyWasSaved(QString)), target, SLOT(hierarchyWasSaved(QString)));
    connect(this, SIGNAL(restartWizard()), target, SLOT(restartWizard()));

    // get a description for the type of dataset, removing the first letter (which should be "a")
    QString filetype = QString::fromStdString(schema->schemaTitle().substr(1));

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // ---- save ----
    {
        // add the caption
        layout->addWidget(new QLabel("Press this button to save the" + filetype + " into the same file:"));

        // add the button
        _saveButton = new QPushButton("Save");
        auto buttonLayout = new QHBoxLayout;
        layout->addLayout(buttonLayout);
        buttonLayout->addWidget(_saveButton, 1);

        // add the filepath label
        _filepathLabel = new QLabel();
        _filepathLabel->setWordWrap(true);
        buttonLayout->addWidget(_filepathLabel, 4);

        // connect the button
        connect(_saveButton, SIGNAL(clicked()), this, SLOT(save()));
    }

    // ---- save as ----
    {
        // add the caption
        layout->addWidget(new QLabel("Press this button to save the " + filetype + " into a new file:"));

        // add the button
        _saveAsButton = new QPushButton("Save As...");
        auto buttonLayout = new QHBoxLayout;
        layout->addLayout(buttonLayout);
        buttonLayout->addWidget(_saveAsButton, 1);
        buttonLayout->addStretch(4);

        // connect the button
        connect(_saveAsButton, SIGNAL(clicked()), this, SLOT(saveAs()));
    }

    // ---- restart ----
    {
        // add the caption
        layout->addWidget(new QLabel("Press this button to restart the wizard:"));

        // add the button
        _restartButton = new QPushButton("Restart");
        auto buttonLayout = new QHBoxLayout;
        layout->addLayout(buttonLayout);
        buttonLayout->addWidget(_restartButton, 1);
        buttonLayout->addStretch(4);

        // connect the button
        connect(_restartButton, SIGNAL(clicked()), this, SLOT(restart()));
    }

    // ---- quit ----
    {
        // add the caption
        layout->addWidget(new QLabel("Press this button or close the window to exit the wizard:"));

        // add the button
        _quitButton = new QPushButton("Quit");
        auto buttonLayout = new QHBoxLayout;
        layout->addLayout(buttonLayout);
        buttonLayout->addWidget(_quitButton, 1);
        buttonLayout->addStretch(4);

        // connect the button
        connect(_quitButton, SIGNAL(clicked()), this, SLOT(quit()));
    }

    // --------

    // finalize the layout and assign it to ourselves
    layout->addStretch();
    setLayout(layout);

    // enable/disable save button and fill the filepath label
    updateSaveInfo();
}

////////////////////////////////////////////////////////////////////

void SaveWizardPane::save()
{
    // if the previous path is known, save again; otherwise ask anyway
    if (!_filepath.isEmpty())
        saveToFile(_filepath);
    else
        saveAs();
}

////////////////////////////////////////////////////////////////////

void SaveWizardPane::saveAs()
{
    // get a file path from the user
    QString directory =
        !_filepath.isEmpty() ? _filepath : QStandardPaths::writableLocation(QStandardPaths::DesktopLocation);
    QString caption = qApp->applicationName() + " - Save " + QString::fromStdString(_schema->schemaTitle());
    QString extension = QString::fromStdString(_schema->schemaExtension());
    QString filter = extension + " files (*." + extension + ")";
    QString filepath = QFileDialog::getSaveFileName(this, caption, directory, filter);

    // if the user did not cancel, save the file
    if (!filepath.isEmpty())
    {
        // add filename extension if needed
        if (!filepath.toLower().endsWith("." + extension.toLower())) filepath += "." + extension;
        saveToFile(filepath);
    }
}

////////////////////////////////////////////////////////////////////

void SaveWizardPane::restart()
{
    if (_dirty)
    {
        auto ret = QMessageBox::warning(this, qApp->applicationName(), "Do you want to discard your unsaved changes?",
                                        QMessageBox::Discard | QMessageBox::Cancel, QMessageBox::Cancel);
        if (ret == QMessageBox::Cancel) return;
    }

    // tell the wizard to restart
    emit restartWizard();
}

////////////////////////////////////////////////////////////////////

void SaveWizardPane::quit()
{
    qApp->closeAllWindows();
}

////////////////////////////////////////////////////////////////////

void SaveWizardPane::saveToFile(QString filepath)
{
    // save the hierarchy in the specified file
    XmlHierarchyWriter::write(_root, _schema, filepath.toStdString(),
                              (qApp->applicationName() + " " + qApp->applicationVersion()).toStdString());
    _filepath = filepath;
    _dirty = false;

    // notify the target
    emit hierarchyWasSaved(filepath);

    // update our UI
    updateSaveInfo();
}

////////////////////////////////////////////////////////////////////

void SaveWizardPane::updateSaveInfo()
{
    _saveButton->setEnabled(_dirty);
    _filepathLabel->setText(_filepath);
}

////////////////////////////////////////////////////////////////////
