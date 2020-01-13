/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpenWizardPane.hpp"
#include "FatalError.hpp"
#include "Item.hpp"
#include "ItemUtils.hpp"
#include "SchemaDef.hpp"
#include "XmlHierarchyCreator.hpp"
#include <QApplication>
#include <QFileDialog>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <QStandardPaths>
#include <QVBoxLayout>

////////////////////////////////////////////////////////////////////

OpenWizardPane::OpenWizardPane(const SchemaDef* schema, QString filepath, bool dirty, QObject* target)
    : _schema(schema), _filepath(filepath), _dirty(dirty)
{
    // connect ourselves to the target
    connect(this, SIGNAL(hierarchyWasLoaded(Item*, QString)), target, SLOT(hierarchyWasLoaded(Item*, QString)));

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the caption
    layout->addWidget(new QLabel("Press this button to load " + QString::fromStdString(schema->schemaTitle())));

    // add the button
    _openButton = new QPushButton("Open...");
    auto buttonLayout = new QHBoxLayout;
    layout->addLayout(buttonLayout);
    buttonLayout->addWidget(_openButton, 1);

    // add the filepath label
    _filepathLabel = new QLabel();
    _filepathLabel->setWordWrap(true);
    _filepathLabel->setText(_filepath);
    buttonLayout->addWidget(_filepathLabel, 4);

    // connect the button
    connect(_openButton, SIGNAL(clicked()), this, SLOT(open()));

    // finalize the layout and assign it to ourselves
    layout->addStretch();
    setLayout(layout);
}

////////////////////////////////////////////////////////////////////

void OpenWizardPane::open()
{
    // if the current hierarchy is dirty, give the user a chance to opt out
    if (_dirty)
    {
        auto ret = QMessageBox::warning(this, qApp->applicationName(), "Do you want to discard your unsaved changes?",
                                        QMessageBox::Discard | QMessageBox::Cancel, QMessageBox::Cancel);
        if (ret == QMessageBox::Cancel) return;
    }

    // get a file path from the user
    QString directory =
        !_filepath.isEmpty() ? _filepath : QStandardPaths::writableLocation(QStandardPaths::DesktopLocation);
    QString caption = qApp->applicationName() + " - Open " + QString::fromStdString(_schema->schemaTitle());
    QString extension = QString::fromStdString(_schema->schemaExtension());
    QString filter = extension + " files (*." + extension + ")";
    QString filepath = QFileDialog::getOpenFileName(this, caption, directory, filter);

    // if the user did not cancel, load the file
    if (!filepath.isEmpty())
    {
        // attempt to load the dataset from the specified file
        std::unique_ptr<Item> root;
        vector<string> messageLines;
        try
        {
            root = XmlHierarchyCreator::readFile(_schema, filepath.toStdString());
        }
        catch (FatalError error)
        {
            messageLines = error.message();
        }

        // if successful, process the result
        if (root)
        {
            _filepath = filepath;
            _dirty = false;

            // set all properties to the "user-configured-this-property" state
            // so that property wizard panes won't replace the value by a fresh default
            ItemUtils::setHierarchyConfigured(_schema, root.get());

            // set all items to the "complete" state
            // so that the wizard doesn't force users to descend into each subitem in item list
            ItemUtils::setHierarchyComplete(root.get());

            // notify the target, handing over ownership for the new hierarchy
            emit hierarchyWasLoaded(root.release(), _filepath);

            // update our UI
            _filepathLabel->setText(_filepath);
        }

        // otherwise alert the user
        else
        {
            QString message = "An error occurred while opening or loading the file:";
            if (messageLines.size() > 0) message += "\n" + QString::fromStdString(messageLines[0]);
            if (messageLines.size() > 1) message += "\n" + QString::fromStdString(messageLines[1]);
            QMessageBox::critical(this, qApp->applicationName(), message, QMessageBox::Ok);
        }
    }
}

////////////////////////////////////////////////////////////////////
