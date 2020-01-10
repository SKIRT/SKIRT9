/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CreateRootWizardPane.hpp"
#include "SchemaDef.hpp"
#include <QButtonGroup>
#include <QLabel>
#include <QRadioButton>
#include <QVBoxLayout>
#include <QVariant>

////////////////////////////////////////////////////////////////////

CreateRootWizardPane::CreateRootWizardPane(const SchemaDef* schema, string initialType, QObject* target)
{
    // connect ourselves to the target
    connect(this, SIGNAL(rootTypeChanged(string)), target, SLOT(setRootType(string)));

    // if there is currently no root, or it is of the wrong type, create one using the first option in the list
    if (initialType.empty() || !schema->inherits(initialType, schema->schemaType()))
    {
        initialType = schema->descendants(schema->schemaType())[0];
        emit rootTypeChanged(initialType);
    }

    // create the layout so that we can add stuff one by one
    auto layout = new QVBoxLayout;

    // add the question
    QString baseTitle = QString::fromStdString(schema->title(schema->schemaType()));
    layout->addWidget(new QLabel("Select one of the following options for the type of " + baseTitle + ":"));

    // add the radio buttons reflecting the possible choices, putting them into a button group as well
    auto buttonGroup = new QButtonGroup;
    for (auto choiceType : schema->descendants(schema->schemaType()))
    {
        QString choiceTitle = QString::fromStdString(schema->title(choiceType));
        if (!choiceTitle.isEmpty()) choiceTitle.replace(0, 1, choiceTitle[0].toUpper());
        auto choiceButton = new QRadioButton(choiceTitle);
        choiceButton->setFocusPolicy(Qt::NoFocus);
        buttonGroup->addButton(choiceButton);
        layout->addWidget(choiceButton);

        // associate the item type corresponding to this button with the button object
        choiceButton->setProperty("choiceType", QString::fromStdString(choiceType));
        choiceButton->setStatusTip(QString::fromStdString(choiceType));

        // select the button corresponding to the initial choice
        if (choiceType == initialType) choiceButton->setChecked(true);
    }

    // connect the button group to ourselves
    connect(buttonGroup, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(selectTypeFor(QAbstractButton*)));

    // finalize the layout and assign it to ourselves
    layout->addStretch();
    setLayout(layout);
}

////////////////////////////////////////////////////////////////////

void CreateRootWizardPane::selectTypeFor(QAbstractButton* button)
{
    emit rootTypeChanged(button->property("choiceType").toString().toStdString());
}

////////////////////////////////////////////////////////////////////
