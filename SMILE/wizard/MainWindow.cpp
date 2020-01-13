/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MainWindow.hpp"
#include "SchemaDef.hpp"
#include "WizardEngine.hpp"
#include <QApplication>
#include <QDesktopServices>
#include <QFileInfo>
#include <QFileOpenEvent>
#include <QHBoxLayout>
#include <QLabel>
#include <QMenu>
#include <QMessageBox>
#include <QPushButton>
#include <QScrollArea>
#include <QSettings>
#include <QStatusBar>
#include <QUrl>
#include <QVBoxLayout>

////////////////////////////////////////////////////////////////////

MainWindow::MainWindow()
{
    // setup the window, restoring previous position and size
    setMinimumSize(600, 400);
    readSettings();

    // create the status bar
    statusBar()->addPermanentWidget(new QLabel(qApp->applicationName() + " " + qApp->applicationVersion()));

    // create the label that will display the current path in the SMILE dataset
    _pathLabel = new QLabel;
    _pathLabel->setWordWrap(true);
    auto pathLayout = new QHBoxLayout;
    pathLayout->addWidget(_pathLabel);
    auto pathArea = new QFrame;
    pathArea->setFrameStyle(QFrame::StyledPanel);
    pathArea->setLayout(pathLayout);

    // create the pane that holds the buttons to drive the wizard
    auto advanceButton = new QPushButton("Continue");
    auto retreatButton = new QPushButton("Back");
    auto buttonGroupLayout = new QHBoxLayout;
    buttonGroupLayout->addWidget(retreatButton);
    buttonGroupLayout->addWidget(advanceButton);
    auto buttonLayout = new QHBoxLayout;
    buttonLayout->addWidget(pathArea, 2);
    buttonLayout->addLayout(buttonGroupLayout, 1);

    // create the pane that will hold the wizard UI
    _wizardPane = new QWidget;
    _wizardLayout = new QHBoxLayout;
    _wizardLayout->setSizeConstraint(QLayout::SetMinAndMaxSize);
    _wizardLayout->addWidget(_wizardPane);
    auto wizardFrame = new QFrame;
    wizardFrame->setLayout(_wizardLayout);
    wizardFrame->setToolTip("Right-click for help");
    auto wizardArea = new QScrollArea;
    wizardArea->setFrameStyle(QFrame::StyledPanel);
    wizardArea->setWidgetResizable(true);
    wizardArea->setWidget(wizardFrame);

    // create the central area
    auto centralLayout = new QVBoxLayout;
    auto centralArea = new QWidget;
    centralLayout->addWidget(wizardArea);
    centralLayout->addLayout(buttonLayout);
    centralArea->setLayout(centralLayout);
    setCentralWidget(centralArea);

    // create the wizard engine and connect it into our UI
    _wizard = new WizardEngine(this);
    connect(advanceButton, SIGNAL(clicked()), _wizard, SLOT(advance()));
    connect(retreatButton, SIGNAL(clicked()), _wizard, SLOT(retreat()));
    connect(_wizard, SIGNAL(canAdvanceChangedTo(bool)), advanceButton, SLOT(setEnabled(bool)));
    connect(_wizard, SIGNAL(canRetreatChangedTo(bool)), retreatButton, SLOT(setEnabled(bool)));
    connect(_wizard, SIGNAL(stateChanged()), this, SLOT(replaceWizardPane()));
    connect(_wizard, SIGNAL(titleChanged()), this, SLOT(updateTitle()));
    connect(_wizard, SIGNAL(dirtyChanged()), this, SLOT(updateDirtyState()));
    _wizard->emitStateChanged();

    // set the window title (assumes the wizard has been setup)
    updateTitle();
    updateDirtyState();
}

////////////////////////////////////////////////////////////////////

void MainWindow::readSettings()
{
    QSettings settings;
    QPoint pos = settings.value("mainpos", QPoint(200, 200)).toPoint();
    QSize size = settings.value("mainsize", QSize(400, 400)).toSize();
    resize(size);
    move(pos);
}

////////////////////////////////////////////////////////////////////

void MainWindow::writeSettings()
{
    QSettings settings;
    settings.setValue("mainpos", pos());
    settings.setValue("mainsize", size());
}

////////////////////////////////////////////////////////////////////

void MainWindow::replaceWizardPane()
{
    _wizardLayout->removeWidget(_wizardPane);
    delete _wizardPane;
    _wizardPane = _wizard->createPane();
    _wizardLayout->addWidget(_wizardPane);

    _pathLabel->setText(_wizard->hierarchyPath());
}

////////////////////////////////////////////////////////////////////

void MainWindow::updateTitle()
{
    QString filepath = _wizard->filepath();
    QString file = filepath.isEmpty() ? QString("Untitled") : QFileInfo(filepath).fileName();
    QString appName = qApp->applicationName();
    setWindowTitle(file + "[*] - " + appName);
}

////////////////////////////////////////////////////////////////////

void MainWindow::updateDirtyState()
{
    setWindowModified(_wizard->isDirty());
}

////////////////////////////////////////////////////////////////////

void MainWindow::browseUrl()
{
    if (sender())
    {
        QString url = sender()->property("URL").toString();
        if (!url.isEmpty()) QDesktopServices::openUrl(QUrl(url));
    }
}

////////////////////////////////////////////////////////////////////

void MainWindow::keyPressEvent(QKeyEvent* event)
{
    switch (event->key())
    {
        case Qt::Key_Left:
        case Qt::Key_PageUp:
            if (_wizard->canRetreat()) _wizard->retreat();
            break;
        case Qt::Key_Right:
        case Qt::Key_PageDown:
        case Qt::Key_Enter:
        case Qt::Key_Return:
            if (_wizard->canAdvance()) _wizard->advance();
            break;
        default: QMainWindow::keyPressEvent(event);
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    QString uncapitalize(QString name)
    {
        QString result;
        for (auto c : name)
        {
            if (c.isUpper())
            {
                result += '_';
                result += c.toLower();
            }
            else
            {
                result += c;
            }
        }
        return result;
    }
}

////////////////////////////////////////////////////////////////////

void MainWindow::contextMenuEvent(QContextMenuEvent* event)
{
    QMenu menu(this);
    menu.addAction("Online help")->setEnabled(false);
    menu.addSeparator();

    // add help item for current data set, if it provides a URL;
    // add help item for the widget under the mouse click, if it provides a status tip
    //     in the form "type" or "type:property", or if it is the path label
    auto schema = _wizard->schema();
    if (schema && !schema->schemaUrl().empty())
    {
        QString schemaUrl = QString::fromStdString(schema->schemaUrl());

        auto widget = qApp->widgetAt(event->globalPos());
        if (widget == _pathLabel)
        {
            QStringList segments = _wizard->hierarchyPath().split(" ");
            int n = segments.size();
            if (n > 2 && n % 2 == 1)
            {
                // add help item for the property represented by the last segments
                string property = segments[n - 1].toStdString();
                string type = segments[n - 3].toStdString();
                QString propertyTitle = QString::fromStdString(schema->propertyTitle(type, property));
                QString definingType = QString::fromStdString(schema->definingType(type, property));
                QString propertyUrl = schemaUrl + "/" + "class" + uncapitalize(definingType) + ".html";
                menu.addAction("About " + propertyTitle + "...", this, SLOT(browseUrl()))
                    ->setProperty("URL", propertyUrl);

                // add help items for all of the types in the hierarchy
                for (int i = n - 3; i >= 0; i -= 2)
                {
                    string type = segments[i].toStdString();
                    QString typeTitle = QString::fromStdString(schema->title(type));
                    QString typeUrl = schemaUrl + "/" + "class" + uncapitalize(QString::fromStdString(type)) + ".html";
                    menu.addAction("About " + typeTitle + "...", this, SLOT(browseUrl()))->setProperty("URL", typeUrl);
                }
            }
        }
        else if (widget && !widget->statusTip().isEmpty())
        {
            // extract the item type and optional property name from the status tip
            QStringList segments = widget->statusTip().split(" : ");
            string type = segments[0].toStdString();
            string property = segments.size() > 1 ? segments[1].toStdString() : "";

            if (!property.empty())
            {
                // add help item for the property represented by the current widget
                QString propertyTitle = QString::fromStdString(schema->propertyTitle(type, property));
                QString definingType = QString::fromStdString(schema->definingType(type, property));
                QString propertyUrl = schemaUrl + "/" + "class" + uncapitalize(definingType) + ".html";
                menu.addAction("About " + propertyTitle + "...", this, SLOT(browseUrl()))
                    ->setProperty("URL", propertyUrl);
            }

            // add help item for the type represented by the current widget
            QString typeTitle = QString::fromStdString(schema->title(type));
            QString typeUrl = schemaUrl + "/" + "class" + uncapitalize(QString::fromStdString(type)) + ".html";
            menu.addAction("About " + typeTitle + "...", this, SLOT(browseUrl()))->setProperty("URL", typeUrl);
        }

        // add help item for the current data set
        menu.addAction("About " + QString::fromStdString(schema->schemaName()) + "...", this, SLOT(browseUrl()))
            ->setProperty("URL", schemaUrl);
    }

    // in any case, add help item for MakeUp itself
    menu.addAction("About MakeUp...", this, SLOT(browseUrl()))->setProperty("URL", "http://www.skirt.ugent.be/makeup9");
    menu.exec(event->globalPos());
}

////////////////////////////////////////////////////////////////////

void MainWindow::closeEvent(QCloseEvent* event)
{
    if (_wizard->isDirty() && !_acceptedCloseEvent)
    {
        auto ret = QMessageBox::warning(this, qApp->applicationName(), "Do you want to discard your unsaved changes?",
                                        QMessageBox::Discard | QMessageBox::Cancel, QMessageBox::Cancel);
        if (ret == QMessageBox::Cancel)
        {
            event->ignore();
            return;
        }
        _acceptedCloseEvent = true;
    }

    writeSettings();
    event->accept();
}

////////////////////////////////////////////////////////////////////
