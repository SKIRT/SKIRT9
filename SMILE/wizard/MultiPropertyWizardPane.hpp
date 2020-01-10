/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIPROPERTYWIZARDPANE_HPP
#define MULTIPROPERTYWIZARDPANE_HPP

#include "Basics.hpp"
#include <QHash>
#include <QWidget>
class PropertyWizardPane;
class QVBoxLayout;

////////////////////////////////////////////////////////////////////

/** A MultiPropertyWizardPane instance combines one or more PropertyWizardPane instances into a
    single pane. It is used to show the user interface for multiple non-compound properties on a
    single pane rather than on consecutive wizard panes. */
class MultiPropertyWizardPane : public QWidget
{
    Q_OBJECT

    // ============= Construction and Destruction =============

public:
    /** The default (and only) constructor connects the multiPropertyValidChanged() and
        multiPropertyValueChanged() signals to the specified target object. These signals are
        emitted by the multi-pane in response to the corresponding signals emitted by the subpanes
        being managed. */
    explicit MultiPropertyWizardPane(QObject* target);

    /** This function adds the specified wizard pane to the multi-pane's layout and includes it in
        the list of subpanes being managed by the multi-pane. The multi-wizard pane assumes
        ownership of the provided PropertyWizardPane instance. */
    void addPane(PropertyWizardPane* pane);

    // ==================== Event Handling ====================

private:
    /** This function hides subpanes that are currently silent for the purpose of this wizard and
        shows the other subpanes. */
    void updateVisibility();

    /** This function emits the multiPropertyValidChanged signal with the aggregate validity of the
        visible subpanes. */
    void emitValidity();

protected:
    /** This function ensures that the first focus-enabled widget in the pane receives the focus
        when the pane is shown. */
    void showEvent(QShowEvent* event) override;

signals:
    /** This signal is emitted when the valid state of the multi-pane may have changed, i.e. when
        the valid state of a managed subpane changes. */
    void multiPropertyValidChanged(bool valid);

    /** This signal is emitted when the value of one of the properties being handled by multi-pane
        pane has changed. */
    void multiPropertyValueChanged();

public slots:
    /** This function emits a multiPropertyValidChanged() signal in response to a signal received
        from one of the subpanes being managed. It tracks the validity state for each of the
        subpanes so that it can pass on an aggregate value, i.e. the multi-pane is valid only when
        all subpanes are valid. */
    void setPropertyValid(bool valid);

    /** This function emits a multiPropertyValueChanged() signal in response to a signal received
        from one of the subpanes being managed. */
    void hierarchyWasChanged();

    // ================== Data Members ====================

private:
    QVBoxLayout* _multiLayout;
    QList<PropertyWizardPane*> _subPanes;
    QHash<QObject*, bool> _subPaneState;
};

////////////////////////////////////////////////////////////////////

#endif
