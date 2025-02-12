#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# For use on any Unix system, including Mac OS X and Linux
#
# Execute this script with "git" as default directory to
# automatically reformat all C++ code in the SKIRT project
#
# Requires clang-format version 18.1 to be installed
# in the default path.
#

# --------------------------------------------------------------------

# Look for clang-format-18 or clang-format
CLANGFORMATPATH="$(which clang-format-18)"
if [ "$CLANGFORMATPATH" == "" ]
then
    CLANGFORMATPATH="$(which clang-format)"
fi

# Exit with an error if we don't find it
if [ "$CLANGFORMATPATH" == "" ]
then
    echo
    echo Fatal error: there is no clang-format in the default path
    echo
    exit
fi

# Verify the clang-format version
CLANGFORMATVERSION=$($CLANGFORMATPATH -version)
if ! [[ $CLANGFORMATVERSION == *"18.1."* ]]
then
    echo
    echo Fatal error: $CLANGFORMATPATH is not version 18.1 but
    echo $CLANGFORMATVERSION
    echo
    exit
fi

# Verify that we are probably in the SKIRT git directory
if ! [[ -f "makeSKIRT.sh" ]]
then
    echo
    echo Fatal error: the current directory does not seem to be the SKIRT git directory
    echo
    exit
fi

# Format all .hpp and .cpp files, which skips third-party source code because it has a different filename extension
# (we use xargs so that multiple invocations of clang-format can run in parallel)
echo Using $CLANGFORMATPATH -- $CLANGFORMATVERSION...
find . \( -name '*.hpp' -or -name '*.cpp' \) -print0 | xargs -0L1 -P0 $CLANGFORMATPATH -style=file -i 2> formaterror.txt

# If the error file is nonempty, show its contents before removing it
if [ -s formaterror.txt ]
then
    cat formaterror.txt
    rm -f formaterror.txt
    echo Failed!
else
    rm -f formaterror.txt
    echo Done
fi
