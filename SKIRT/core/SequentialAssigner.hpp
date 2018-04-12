/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SEQUENTIALASSIGNER_HPP
#define SEQUENTIALASSIGNER_HPP

#include "ProcessAssigner.hpp"

//////////////////////////////////////////////////////////////////////

/** The SequentialAssigner class is a subclass of the ProcessAssigner class, representing objects
    that assign work to different processes. The SequentialAssigner does this by dividing the work
    (consisting of many parts) in sequential sections, where each section contains more or less the
    same number of parts. Then, each process is assigned to a different section, according to their
    rank. If a certain method in another class incorporates an object of the SequentialAssigner
    class for performing a set of tasks (or parts of work), each process will execute a different
    subset of these tasks. After performing this work in parallel, communication is typically
    needed to accumulate the results stored at different processes. */
class SequentialAssigner : public ProcessAssigner
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by SKIRT classes that wish to hard-code the creation of a
        new ProcessAssigner object of this type. Before the constructor returns, the newly created
        object is hooked up as a child to the specified parent in the simulation hierarchy (so it
        will automatically be deleted) and the setup of the ProcessAssigner base class is invoked,
        which sets the _comm attribute that points to the object of type PeerToPeerCommunicator
        that is found in the simulation hierarchy. As a second argument, it takes the number of
        parts of work \f$n\f$ that need to be performed. Based on this number, the number of
        processes and the rank of this process, the number of tasks that are assigned to this
        process are determined and stored in the \c _assigned member. The algorithm of this function
        goes as follows. First, from the PeerToPeerCommunicator object, the rank \f$i\f$ and size
        \f$N_{P}\f$ are acquired and stored in temporary variables. Then, the quotient and the
        remainder of the integer division of \f$n\f$ and \f$N_{P}\f$ are calculated: \f[ q = \left
        \lfloor{\frac{n}{N_{P}}}\right \rfloor \f] \f[ r = n \bmod{N_{P}} \f] These two values are
        respectively stored in the private \c _quotient and \c _remainder members, for use of the
        rankForIndex() function. Next, based on \f$q\f$ and \f$r\f$, the number of values assigned
        to the process is determined. This is done based on the following simple principle: <ol>
        <li> First, hand out \f$q\f$ values to each process. <li> Then, give the first \f$r\f$
        process one value extra. </ol> With the above method, all \f$n\f$ values get
        assigned to a process, and the work load is maximally balanced (the difference in number of
        tasks between two arbitrary processes is no more than one). The number of values (or tasks)
        as calculated based on the method described above, is stored in the _assigned member. The
        last thing the assign function calculates is the starting value for the particular process.
        This is done by differentiating between two seperate cases: <ul> <li> Either the rank
        \f$i\f$ of the process is smaller than \f$r\f$, <li> or \f$i\f$ is greater than or equal to
        \f$r\f$. </ul> In the former case, the starting index is given by: \f[ t_0 = i \cdot ( q +
        1 ) \f] In the latter case, this index is given by: \f[ t_0 = r \cdot ( q + 1) + (i - r)
        \cdot q \f] The resulting value is stored in the private \c _start member of this class,
        for use in the absoluteIndex() function. */
    explicit SequentialAssigner(size_t size, SimulationItem* parent);

    //======================== Other Functions =======================

public:
    /** This function takes the relative index of a certain part of the work assigned to this
        process as an argument and returns the absolute index of that part, a value from zero to
        the total amount of parts that need to be executed in the simulation. This is done by
        adding the relative index to the index of the starting value for this process, stored in
        the \c _start member. */
    size_t absoluteIndex(size_t relativeIndex) const override;

    /** This function takes the absolute index of a certain part of the work as an argument and
        returns the relative index of that part, a value from zero to the number of parts that were
        assigned to this process, \c _assigned. This is done by subtracting the index of the
        starting value for this process, _start, from the absolute index. */
    size_t relativeIndex(size_t absoluteIndex) const override;

    /** This function returns the rank of the process that is assigned to a certain part of the
        work. This part is identified by its absolute index, passed to this function as an
        argument. The algorithm of this function differentiates between two cases: <ol> <li> The
        index \f$t\f$ is smaller than \f$r \cdot (q + 1) \f$, <li> or \f$i\f$ is greater than or
        equal to \f$r \cdot (q + 1) \f$. </ol> In the former case, the rank \f$j\f$ is found by the
        following formula: \f[ j = \frac{t}{q + 1} \f] In the latter case, the rank is determined
        by: \f[ j = r + \frac{t^{*}}{q} \f] where \f$ t^{*} = t - r \cdot (q + 1) \f$ */
    int rankForIndex(size_t index) const override;

    //======================== Data Members ========================

private:
    size_t _start{0};          // the index of the first value assigned to this process
    size_t _quotient{0};       // the quotient
    size_t _remainder{0};      // the remainder
};

////////////////////////////////////////////////////////////////////

#endif
