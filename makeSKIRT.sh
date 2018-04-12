#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# For use on any Unix system, including Mac OS X and Linux
#
# Execute this script with "git" as default directory to
# build a release version of skirt in the "release" directory
# using your current local copy of the code
#
# By default the build uses a single thread; you can specify the
# number of parallel threads as the first command line argument
#

# --------------------------------------------------------------------

# Look for cmake in the default path; exit with an error if we don't find it
CMAKEPATH="$(which cmake)"
if [ "$CMAKEPATH" == "" ]
then
echo
echo Fatal error: there is no cmake in the default path
echo
exit
else
echo
echo Using $CMAKEPATH to generate build files
echo
fi

# Generate the build files
$CMAKEPATH -E make_directory ../release
$CMAKEPATH -E chdir ../release $CMAKEPATH -DCMAKE_BUILD_TYPE:STRING=Release -L ../git
echo

# Perform the build
make -j ${1:-1} -C ../release
echo
