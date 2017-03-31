#!/bin/sh
# Uninstallation-type script for Camargue
# This script cleans all files, removes the config header, and replaces the
# Makefile in the project directory with a template

if [ -f Makefile ]
then
    make clean
    rm Makefile
fi

if [ -f includes/config.hpp ]
then
    rm includes/config.hpp
fi

cp scripts/Makefile.template Makefile
