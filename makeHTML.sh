#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# Execute this script with "git" as default directory; use on Mac OS X only
#

# generate the html documentation in a folder next to the git folder
/Applications/Doxygen.app/Contents/Resources/doxygen Docs/doxygen/html.doxygen
