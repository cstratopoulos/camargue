# Uninstallation-type script for Camargue
# This script cleans all files, removes the config header, and replaces the
# Makefile in the project directory with a template

make clean
rm includes/config.hpp
rm Makefile
cp scripts/Makefile.template Makefile
