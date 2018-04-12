#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# For use on any Unix system, including Mac OS X and Linux
#
# Execute this script with "git" as default directory to configure
# the release build options for your current local copy of the code
#
# To adjust a build option, enter:    ./configSKIRT.sh <option-name>=<new-value>
# For example, to enable MPI, enter:  ./configSKIRT.sh BUILD_WITH_MPI=ON
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

# Assemble the user options from the script arguments
USEROPTIONS=""
for ARGUMENT in "$@"
do
    USEROPTIONS="$USEROPTIONS -D$ARGUMENT"
done

# Generate the build files
$CMAKEPATH -E make_directory ../release
$CMAKEPATH -E chdir ../release $CMAKEPATH $USEROPTIONS -DCMAKE_BUILD_TYPE:STRING=Release -L ../git

# Provide instructions to the user
echo
echo "To adjust build options, enter:     ./configSKIRT.sh <option1>=<value1> <option2>=<value2> ..."
echo "For example, to enable MPI, enter:  ./configSKIRT.sh BUILD_WITH_MPI=ON"
echo
