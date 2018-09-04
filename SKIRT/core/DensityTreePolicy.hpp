/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DENSITYTREEPOLICY_HPP
#define DENSITYTREEPOLICY_HPP

#include "TreePolicy.hpp"

//////////////////////////////////////////////////////////////////////

/** DensityTreePolicy represents the configurable options and the corresponding implementation
    mechanisms for constructing spatial tree grids based on the density distribution of the media
    in the medium system.

    This policy offers several options for configuring the recursive subdivision of the
    hierarchical tree. First of all, the minimum and maximum tree subdvision levels (actually
    offered by the base class) override the other subdvision criteria described below. Tree nodes
    are always subdivided up to the minimum level, and nodes are never subdivided beyond the
    maximum level, regardless of the outcome of the other criteria.

    The remaing subdvision criteria consist of the maximum mass fraction \f$\delta_\text{max}\f$
    for each material type (dust, electrons, gas), the maximum diagonal dust optical depth
    \f$\tau_{\lambda,\text{max}}\f$ at wavelength \f$\lambda\f$, and the maximum dust density
    dispersion \f$q_\text{max}\f$. A node is subdivided as long as the value calculated for the
    node for one or more of these five criteria exceeds the corresponding configured maximum value.
    A criterion is automatically disabled if the corresponding material type is not present in the
    model, and it can be explicitly disabled by configuring a zero maximum value. Configuring an
    impossibly high maximum value has the same effect, but may require substantial calculation to
    verify the criterion for each node.

    We first discuss the three criteria related to dust. For this material type, we use mass and
    mass density (as opposed to number and number density) because it is the appropriate quantity
    for dust in case multiple dust medium components have a different mass per hydrogen atom value.
    The total dust mass in the model, \f$M_\text{model}\f$, and the dust density at a given
    position, \f$\rho(\bf{r})\f$, are obtained by summing the corresponding quantity for each dust
    medium. The average dust density \f$\rho\f$ inside a given node is estimated from density
    samples in \f$N\f$ random positions \f$\bf{r}_i\f$ distributed uniformly across the volume
    \f$V\f$ of the node: \f[\rho = \frac{1}{N}\,\sum_{i=1}^{N}\rho(\bf{r}_i) .\f] The fraction of
    the mass \f$\delta\f$ within the node is then easily found as \f[\delta = \frac{\rho V}
    {M_\text{model}}.\f]

    The estimated optical depth \f$\tau_\lambda\f$ at wavelength \f$\lambda\f$ across the diagonal
    \f$\Delta s\f$ of a node can now be expressed as \f[\tau_\lambda = \kappa_\lambda
    \,\rho\,\Delta s\f] where \f$\kappa_\lambda\f$ is a representative extinction mass coefficient
    for the dust in the medium. For the sake of performance this value is assumed to be constant
    across the spatial domain, and it is determined as the plain average of the
    \f$\kappa_\lambda\f$ values of the dust media components, taken at the origin of the model
    coordinate system.

    Finally, a measure for the dust density dispersion \f$q\f$ within the node is determined as \f[
    q = \begin{cases} \;\dfrac{\rho_{\text{max}}-\rho_{\text{min}}}{\rho_{\text{max}}} &
    \quad\text{if $\rho_{\text{max}}>0$,} \\ \;0 & \quad\text{if $\rho_{\text{max}}=0$.}
    \end{cases} \f] where \f$\rho_{\text{min}}\f$ and \f$\rho_{\text{max}}\f$ are the smallest and
    largest sampled density values from the list of \f$N\f$ sampled positions in the node. The
    quantity \f$q\f$ is a simple measure for the uniformity of the density within the node: for a
    constant density, \f$q=0\f$, whereas \f$q\f$ approaches 1 if a steep gradient is present. The
    special case \f$\rho_{\text{max}}=0\f$ is included because it is possible that a node is empty,
    in which case the uniform value \f$q=0\f$ should be returned. With a configured value of
    \f$0<q_\text{max}<1\f$, nodes that contain a sharp edge with empty space one side will continue
    to be subdivided for ever, because such cells have a density dispersion measure of \f$q=1\f$.
    It is thus important to always specify a reasonable maximum subdivision level when using this
    subdivision criterion.

    For electrons and for gas, only the maximum mass fraction criterion is offered. For these
    material types, the number \f$\mathcal{N}\f$ and number density \f$n\f$ are used instead of the
    mass \f$M\f$ and mass density \f$\rho\f$. Other than this, the procedure is the same as the one
    described for dust. */
class DensityTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(DensityTreePolicy, TreePolicy,
                  "a tree grid construction policy using the medium density distribution")

    // TO DO.

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function constructs the hierarchical tree and all (interconnected) nodes forming the
        tree as described for the corresponding pure virtual function in the base class. The
        implementation for this class loops over the tree subdivision levels (up to the maximum
        level configured in the TreePolicy base class). For each level, the function alternates
        between evaluating all of the nodes (i.e. determining which nodes need subdivision) and
        actually subdividing the nodes that need it.

        These operations are split over two phases because the first one can be parallelized (the
        only output is a Boolean flag), while the second one cannot (the tree structure is updated
        in various ways). Parallelizing the first operation is often meaningful, because
        determining whether a node needs subdivision can be resource-intensive (for example, it may
        require sampling densities in the source distribution). */
    vector<TreeNode*> constructTree(TreeNode* root) override;

private:
    /** This function returns true if the given node needs to be subdivided according to the
        criteria configured for this policy, and false otherwise. The minimum and maximum level are
        not checked, because this function is never called for nodes that don't conform to the
        level criteria. */
    bool needsSubdivide(TreeNode* node);
};

//////////////////////////////////////////////////////////////////////

#endif
