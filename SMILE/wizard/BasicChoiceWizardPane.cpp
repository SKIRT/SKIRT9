/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BasicChoiceWizardPane.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <QApplication>
#include <QButtonGroup>
#include <QFileDialog>
#include <QFileInfo>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <QRadioButton>
#include <QSettings>
#include <QStandardPaths>
#include <QVBoxLayout>

////////////////////////////////////////////////////////////////////

BasicChoiceWizardPane::BasicChoiceWizardPane(bool initialOpenExisting, string initialSchemaName, bool dirty,
                                             QObject* target)
    : _openExisting(initialOpenExisting), _schemaName(initialSchemaName), _dirty(dirty)
{
    // create the overall layout and add a title
    auto layout = new QVBoxLayout;
    layout->addWidget(new QLabel("Welcome to the " + qApp->applicationName() + " wizard!"));

    // create the choice layout and the corresponding button group
    _choiceLayout = new QVBoxLayout;
    layout->addLayout(_choiceLayout);
    _buttonGroup = new QButtonGroup;

    // connect the button group to ourselves, and ourselves to the target
    connect(_buttonGroup, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(setBasicChoice(QAbstractButton*)));
    connect(this, SIGNAL(basicChoiceWasChanged(bool, string, string)), target,
            SLOT(setBasicChoice(bool, string, string)));

    // add white space and a horizontal line
    layout->addStretch();
    auto line = new QFrame;
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);
    layout->addWidget(line);

    // add the caption
    layout->addWidget(new QLabel("Use these buttons to select a new SMILE schema library:"));

    // add the buttons
    auto selectLibraryButton = new QPushButton("Browse for library...");
    auto defaultLibraryButton = new QPushButton("Set built-in library");
    auto buttonLayout = new QVBoxLayout;
    buttonLayout->addWidget(selectLibraryButton);
    buttonLayout->addWidget(defaultLibraryButton);

    // add the library path label
    _libraryLabel = new QLabel();
    _libraryLabel->setWordWrap(true);
    auto pathLayout = new QHBoxLayout;
    pathLayout->addLayout(buttonLayout, 1);
    pathLayout->addWidget(_libraryLabel, 4);
    layout->addLayout(pathLayout);

    // connect the buttons
    connect(selectLibraryButton, SIGNAL(clicked()), this, SLOT(selectLibrary()));
    connect(defaultLibraryButton, SIGNAL(clicked()), this, SLOT(setDefaultLibrary()));

    // assign the layout to ourselves
    setLayout(layout);

    // get the default and current library paths
    findDefaultLibraryPath();
    loadPersistentLibraryPath();
    if (_libraryPath.empty() && !_defaultLibraryPath.empty())
    {
        _libraryPath = _defaultLibraryPath;
        storePersistentLibraryPath();
    }

    // obtain schema information and further complete the user interface
    loadSchemaInfo();
    updateChoiceInterface();
}

////////////////////////////////////////////////////////////////////

void BasicChoiceWizardPane::setBasicChoice(QAbstractButton* /*button*/)
{
    int buttonIndex = _buttonGroup->checkedId();
    bool openExisting = (buttonIndex % 2) != 0;
    string schemaName = _schemaNames[buttonIndex >> 1];

    if (_openExisting != openExisting || _schemaName != schemaName)
    {
        bool update = true;

        // if the current hierarchy is dirty, give the user a chance to cancel the change
        if (_dirty)
        {
            auto ret =
                QMessageBox::warning(this, qApp->applicationName(), "Do you want to discard your unsaved changes?",
                                     QMessageBox::Discard | QMessageBox::Cancel, QMessageBox::Cancel);
            update = (ret == QMessageBox::Discard);
        }

        // if allowed, update the choice
        if (update)
        {
            _openExisting = openExisting;
            _schemaName = schemaName;
            _dirty = false;
            emit basicChoiceWasChanged(_openExisting, _libraryPath, _schemaName);
        }
        // otherwise, revert the selected radio button
        else
        {
            int buttonIndex = (StringUtils::indexOf(_schemaNames, _schemaName) << 1) + (_openExisting ? 1 : 0);
            auto selected = _buttonGroup->button(buttonIndex);
            if (selected) selected->setChecked(true);
        }
    }
}

////////////////////////////////////////////////////////////////////

void BasicChoiceWizardPane::selectLibrary()
{
    // get a new library directory path from the user
    QString directory = !_libraryPath.empty() ? QString::fromStdString(_libraryPath)
                                              : QStandardPaths::writableLocation(QStandardPaths::DesktopLocation);
    QString caption = qApp->applicationName() + " - Select SMILE schema library";
    string newPath = QFileDialog::getExistingDirectory(this, caption, directory).toStdString();
    if (!newPath.empty() && newPath != _libraryPath)
    {
        _libraryPath = newPath;
        storePersistentLibraryPath();
    }

    // obtain schema information and update the user interface
    loadSchemaInfo();
    updateChoiceInterface();
}

////////////////////////////////////////////////////////////////////

void BasicChoiceWizardPane::setDefaultLibrary()
{
    if (_defaultLibraryPath != _libraryPath)
    {
        _libraryPath = _defaultLibraryPath;
        storePersistentLibraryPath();
    }

    // obtain schema information and update the user interface
    loadSchemaInfo();
    updateChoiceInterface();
}

////////////////////////////////////////////////////////////////////

void BasicChoiceWizardPane::findDefaultLibraryPath()
{
    _defaultLibraryPath.clear();

    // get the location of the executable
    QString appPath = qApp->applicationDirPath();

    // iterate over the relative paths
    for (auto path : {"schemas", "../schemas", "../../schemas", "../../../schemas", "../../../../schemas"})
    {
        QFileInfo test(appPath + "/" + path);
        if (test.isDir())
        {
            _defaultLibraryPath = test.canonicalFilePath().toStdString();
            break;
        }
    }
}

////////////////////////////////////////////////////////////////////

void BasicChoiceWizardPane::loadPersistentLibraryPath()
{
    _libraryPath = QSettings().value("SchemaLibraryPath").toString().toStdString();
}

////////////////////////////////////////////////////////////////////

void BasicChoiceWizardPane::storePersistentLibraryPath()
{
    QSettings().setValue("SchemaLibraryPath", QString::fromStdString(_libraryPath));
}

////////////////////////////////////////////////////////////////////

void BasicChoiceWizardPane::loadSchemaInfo()
{
    _error.clear();
    _schemaNames.clear();
    _schemaTitles.clear();

    // list schema files in library directory
    for (string name : System::filesInDirectory(_libraryPath))
    {
        if (StringUtils::endsWith(name, ".smile"))
        {
            string title = SchemaDef::getSchemaTitle(StringUtils::joinPaths(_libraryPath, name));
            if (!title.empty())
            {
                _schemaNames.push_back(name);
                _schemaTitles.push_back(title);
            }
        }
    }
    if (_schemaNames.empty())
    {
        _error = "There are no SMILE schema files in the schema library";
    }
}

////////////////////////////////////////////////////////////////////

void BasicChoiceWizardPane::updateChoiceInterface()
{
    // show the schema library path
    _libraryLabel->setText(QFileInfo(QString::fromStdString(_libraryPath)).absoluteFilePath());

    // remove any basic choice labels and buttons currently in the layout
    for (auto label : _labels) delete label;
    for (auto button : _buttons)
    {
        _buttonGroup->removeButton(button);
        delete button;
    }
    _labels.clear();
    _buttons.clear();

    // if there is a problem with the schema library, show error information
    if (!_error.empty())
    {
        auto label1 = new QLabel("<<< " + QString::fromStdString(_error) + " >>>");
        auto label2 = new QLabel("Please select a SMILE schema library below.");
        _choiceLayout->addWidget(label1);
        _choiceLayout->addWidget(label2);
        _labels << label1 << label2;
    }
    else
    {
        auto label = new QLabel("What would you like to do? Select one of the following options:");
        _choiceLayout->addWidget(label);
        _labels << label;

        // create the radio buttons that allow the user to make the basic choice
        int buttonID = 0;
        for (const auto& title : _schemaTitles)
        {
            auto choice1 = new QRadioButton("Create and configure " + QString::fromStdString(title));
            auto choice2 = new QRadioButton("Open and edit " + QString::fromStdString(title));
            choice1->setFocusPolicy(Qt::NoFocus);
            choice2->setFocusPolicy(Qt::NoFocus);
            _buttonGroup->addButton(choice1, buttonID++);
            _buttonGroup->addButton(choice2, buttonID++);
            _choiceLayout->addWidget(choice1);
            _choiceLayout->addWidget(choice2);
            _buttons << choice1 << choice2;
        }

        // select the initial choice
        int buttonIndex = (StringUtils::indexOf(_schemaNames, _schemaName) << 1) + (_openExisting ? 1 : 0);
        auto selected = _buttonGroup->button(buttonIndex);
        if (selected) selected->setChecked(true);
    }
}

////////////////////////////////////////////////////////////////////
