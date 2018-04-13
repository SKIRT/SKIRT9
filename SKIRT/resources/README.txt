The SKIRT project -- advanced radiative transfer
© Astronomical Observatory, Ghent University

SKIRT first looks for built-in resource files (in this directory and its subdirectories) and then looks for externally
provided resource files (in a "resources" directory next to the SKIRT \c git directory, i.e. outside of the build
tree). This mechanism allows to provide small and frequently-used resource files as part of the SKIRT build tree in the
source code repository, while requiring larger resource files to be downloaded seperately from the SKIRT web site using
a shell script provided for this purpose.
