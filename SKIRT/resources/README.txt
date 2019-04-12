The SKIRT project -- advanced radiative transfer
© Astronomical Observatory, Ghent University

SKIRT first looks for resource files in this directory and then looks for resource files installed in the "resources"
directory next to the SKIRT git directory, i.e. outside of the build tree. This mechanism allows to provide small
and frequently-used resource files as part of the SKIRT build tree in the source code repository, while requiring
larger resource files to be downloaded seperately from the SKIRT web site.

In the current version, however, the SKIRT source code repository contains no resources. Instead, the file
ExpectedResources.txt placed in this directory (next to this README) lists the names and version numbers of the
resource packs that should be downloaded seperately from the SKIRT web site. The "Core" resource pack is required
for the basic operation of SKIRT; other resource packs are optional and must be installed only if the corresponding
SKIRT functionality is actually being used.

The convenience shell script downloadResources.sh provided in the SKIRT git directory reads the ExpectedResources.txt
file and, after confirmation by the user, downloads and installs each of the expected resource packs.

Refer to the SKIRT installation guide on the SKIRT web site (www.skirt.ugent.be) for more information.
