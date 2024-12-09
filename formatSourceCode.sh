#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# For use on any Unix system, including Mac OS X and Linux
#
# Execute this script with "git" as default directory to
# automatically reformat all C++ code in the SKIRT project
#
# Requires clang-format version 16.0 to be installed
# in the default path, in the home directory, or inside Xcode.
#

# --------------------------------------------------------------------

# Look for clang-format or clang-format-16 in various places;
# exit with an error if we don't find it
CLANGFORMATPATH="$(which clang-format-16)"
if [ "$CLANGFORMATPATH" == "" ]
then
    CLANGFORMATPATH="$(which clang-format)"
fi
if [ "$CLANGFORMATPATH" == "" ]
then
    CANDIDATEPATH="$HOME/clang/bin/clang-format-16"
    if [[ -x $CANDIDATEPATH ]]
    then
        CLANGFORMATPATH=$CANDIDATEPATH
    fi
fi
if [ "$CLANGFORMATPATH" == "" ]
then
    CANDIDATEPATH="$HOME/clang/bin/clang-format"
    if [[ -x $CANDIDATEPATH ]]
    then
        CLANGFORMATPATH=$CANDIDATEPATH
    fi
fi
if [ "$CLANGFORMATPATH" == "" ]
then
    CANDIDATEPATH="/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang-format"
    if [[ -x $CANDIDATEPATH ]]
    then
        CLANGFORMATPATH=$CANDIDATEPATH
    fi
fi
if [ "$CLANGFORMATPATH" == "" ]
then
    echo
    echo Fatal error: there is no clang-format in the default path or in ~/clang/bin/
    echo
    exit
fi

# Verify the clang-format version
CLANGFORMATVERSION=$($CLANGFORMATPATH -version)
if ! [[ $CLANGFORMATVERSION == *"16.0."* ]]
then
    echo
    echo Fatal error: $CLANGFORMATPATH is not version 16.0 but
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
echo Using $CLANGFORMATPATH...
find . \( -name '*.hpp' -or -name '*.cpp' \) -exec $CLANGFORMATPATH -style=file -i {} \;
echo Done
