/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef THREADLOCALMEMBER_HPP
#define THREADLOCALMEMBER_HPP

#include "Basics.hpp"
#include <mutex>
#include <unordered_map>
#include <unordered_set>

////////////////////////////////////////////////////////////////////

/** An instance of the ThreadLocalMember<T> class template provides a thread-local data member of
    type T. In other words, there is a seperate copy of the data member for each instance of the
    embedding object, and for each thread referring to the data member. A reference to the
    thread-local copy of the data member can be obtained at the cost of a single std::unordered_map
    look-up, which has constant-time performance on average.

    The type T must be default-constructable, because each copy is constructed automatically in the
    background when a reference to it is requested for the first time. A copy is destructed when
    the thread for which it was created (the \em associated thread) is destructed.

    The class also offers a function to obtain a list of pointers to the copies of the data member
    for all threads that ever referred to it (and that are not yet destructed). For example, when a
    data member is employed as a thread-local buffer, this capability allows a client to "flush"
    the contents of the buffer after multi-threaded operation has ceased (but before the threads
    are actually destroyed).

    The mechanisms and functions of this class are fully thread-safe. However, the all() function
    breaks thread-locality; any access to data member copies associated with other threads should
    be properly synchronized. Also, the pointers returned by the all() function become invalid when
    the threads associated with the referenced data members are destructed.
*/
template<class T> class ThreadLocalMember
{
public:
    /** This function returns a pointer to the thread-local copy of the data member. */
    T* local() { return &repository.map[this]; }

    /** This function returns a list of pointers to the copies of the data member for all threads
        that ever referred to it and that are not yet destructed. The list is in arbitrary order.
        Access to data member copies associated with other threads should be properly synchronized.
        Also, the pointers returned by this function become invalid when the threads associated
        with the referenced data members are destructed */
    vector<T*> all()
    {
        vector<T*> list;
        std::unique_lock<std::mutex> lock(mutex);
        for (auto repo : allRepositories)
            if (repo->map.count(this)) list.push_back(&repo->map.at(this));
        return list;
    }

private:
    // repository that contains all copies of the data member in the local thread
    // and automatically registers itself to the global repository registry
    struct Repository
    {
        Repository()
        {
            std::unique_lock<std::mutex> lock(mutex);
            allRepositories.insert(this);
        }
        ~Repository()
        {
            std::unique_lock<std::mutex> lock(mutex);
            allRepositories.erase(this);
        }
        std::unordered_map<ThreadLocalMember<T>*, T> map;
    };
    static thread_local Repository repository;

    // set containing a pointer to all repositories, across the threads
    using RepositoryRegistry = std::unordered_set<Repository*>;
    static RepositoryRegistry allRepositories;
    static std::mutex mutex;
};

template<class T> thread_local typename ThreadLocalMember<T>::Repository ThreadLocalMember<T>::repository;
template<class T> typename ThreadLocalMember<T>::RepositoryRegistry ThreadLocalMember<T>::allRepositories;
template<class T> std::mutex ThreadLocalMember<T>::mutex;

////////////////////////////////////////////////////////////////////

#endif
