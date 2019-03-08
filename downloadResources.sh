#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# For use on any Unix system, including Mac OS X and Linux
#
# Execute this script with "git" as default directory to download
# the SKIRT 9 resource files provided on the public SKIRT server, and
# place them in the 'resources' directory next to the git directory.
#

# list of archives to be downloaded and extracted (separate filenames with a space)
FILELIST=( StoredTables.tar.gz )

# download using wget (Linux) or curl (Mac OS X)
if which wget >/dev/null
then
    DOWNLOAD="wget --no-check-certificate"
elif which curl >/dev/null
then
    DOWNLOAD="curl --insecure -O"
else
    echo error: no wget or curl available to download files
    exit
fi

# download and extract each of the files in the list
filecount=$((${#FILELIST[@]} - 1))
for i in $(eval echo {0..$filecount})
do
    FILENAME=${FILELIST[$i]}
    echo "----------------------------------------"
    echo downloading $FILENAME ...
    mkdir -p ../resources
    cd ../resources
    $DOWNLOAD https://sciences.ugent.be/skirtextdat/SKIRT9/Resources/$FILENAME
    tar -xvf $FILENAME
    rm $FILENAME
    cd ../git
done

echo "----------------------------------------"
echo Done.
