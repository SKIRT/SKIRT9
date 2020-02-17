#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# For use on any Unix system, including Mac OS X and Linux
#
# Execute this script with "git" as default directory to
# automatically reformat all C++ code in the SKIRT project
#
# Requires clang-format version 9.0.0 to be installed
# in the default path or in the home directory
#
# See https://releases.llvm.org/9.0.0/tools/clang/docs/ClangFormat.html
#

# --------------------------------------------------------------------

# Look for clang-format or clang-format-9 in the default path or in the home directory;
# exit with an error if we don't find it
CLANGFORMATPATH="$(which clang-format-9)"
if [ "$CLANGFORMATPATH" == "" ]
then
    CLANGFORMATPATH="$(which clang-format)"
fi
if [ "$CLANGFORMATPATH" == "" ]
then
    CANDIDATEPATH="$HOME/clang/bin/clang-format-9"
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
    echo
    echo Fatal error: there is no clang-format in the default path or in ~/clang/bin/
    echo
    exit
fi

# Verify the clang-format version
CLANGFORMATVERSION=$($CLANGFORMATPATH -version)
if ! [[ $CLANGFORMATVERSION == *"9.0.0"* ]]
then
    echo
    echo Fatal error: $CLANGFORMATPATH is not version 9.0.0 but
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
