The SKIRT project -- advanced radiative transfer
© Astronomical Observatory, Ghent University

SKIRT first looks for resource files in this directory and then looks for resource files installed in the "resources"
directory next to the SKIRT \c git directory, i.e. outside of the build tree. This mechanism allows to provide small
and frequently-used resource files as part of the SKIRT build tree in the source code repository, while requiring
larger resource files to be downloaded seperately from the SKIRT web site.

In the current version, however, the SKIRT source code repository contains no resources. All resources must be
downloaded seperately from the SKIRT web site using a shell script provided for this purpose. Refer to the SKIRT
installation guide on the SKIRT web site (www.skirt.ugent.be)
