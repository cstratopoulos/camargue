# Installation script for Camargue
# This script is meant to be used to automatically grab external dependencies
# for a Camargue install. As per the README, you should first put a Makefile
# in the camargue/ directory with a C++ compiler and paths to cplex.h and
# cplex.a indicated. After that, this script can be used to try and grab
# desired external dependencies (assuming you are connected to the internet).

usage='
Usage:   ./cmr_install.sh [One or more options, see below]
The first two options are catch-all flags. Zero or one may be specified
-B Bare install: no external dependencies
-F Full install: all external dependencies
If you specified neither -B nor -F, use one or more of the options below to
select externals individually.
-s Safe MIR: the Cook, Dash, Fukasawa, Goycoolea safe Gomory code.
-c catch.hpp: the header for Catch unit tests
-t timsort.hpp: the header for sorting candidate teeth with Timsort
'

if [ "$#" -eq 0 ]; then
    (>&2 echo "Error: no arguments specified")
    (>&2 echo "$usage")
fi

bare=0
full=0
safegmi=0
catch=0
tim=0

while getopts BFsct option
do
    case "$option"
    in
	B) bare=1;;
	F) full=1;;
	s) safegmi=1;;
	c) catch=1;;
	t) tim=1;;
	*) (>&2 echo "$usage")
	    exit 1;;
    esac
done

if [ \( "$bare" -eq 1 \) -a \( "$full" -eq 1 \) ]; then
    (>&2 echo "Error: cannot specify Bare AND Full")
    (>&2 echo "$usage")
    exit 1
fi

if [ \( "$bare" -eq 1 \) -o \( "$full" -eq 1 \) ]; then
    if [ \( "$safegmi" -eq 1 \) -o \( "$catch" -eq 1 \) -o \( "$tim" -eq 1 \) ]
    then
	(>&2 echo "Error: cannot specify catch-all flag with individual flags")
	(>&2 echo "$usage")
	exit 1
    fi
fi

tim_url=https://raw.githubusercontent.com/gfx/cpp-TimSort/master/timsort.hpp
catch_url=https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp
gmi_url=https://www.informs.org/content/download/257165/2428239/file/safemir090309.tar.gz

if [ "$full" -eq 1 ]; then
    safegmi=1
    catch=1
    tim=1
elif [ "$bare" -eq 1 ]; then
    echo "Doing bare install with no externals."
fi

if [ "$catch" -eq 1 ]; then
    echo "Installing with Catch...."
    if [ -e externals/catch.hpp ]; then
	echo "Looks like Catch is already installed."
    else
	curl -# -o "externals/catch.hpp" "$catch_url"
	if [ "$?" -eq 0 ]; then
	    echo "Downloaded catch.hpp"
	else
	    (>&2 echo "Error downloading catch.hpp")
	    exit 1
	fi
    fi
fi

if [ "$tim" -eq 1 ]; then
    echo "Installing with Timsort...."
    if [ -e externals/timsort.hpp ]; then
	echo "Looks like Timsort is already installed."
    else
	curl -# -o "externals/timsort.hpp" "$tim_url"
	if [ "$?" -eq 0 ]; then
	    echo "Downloaded timsort.hpp"
	else
	    (>&2 echo "Error downloading timsort.hpp")
	    exit 1
	fi
    fi
fi



if [ "$safegmi" -eq 1 ]; then
    echo "Installing with Safe GMI code...."
    if [ -e externals/safemir ]; then
	echo "Looks like Safe GMI is already installed."
    else
	curl -# -o "externals/safemir.tar.gz" "$gmi_url"
	if [ "$?" -eq 1 ]; then
	    (>&2 echo "Error downloading safemir.tar.gz")
	fi
	echo "Downloaded compressed Safe GMI code, extracting..."
	tar -xzf externals/safemir.tar.gz  && \
	    rm externals/safemir.tar.gz && mv safemir externals/safemir
	if [ "$?" -eq 1 ]; then
	    (>&2 echo "Error extracting safemir.tar.gz")
	fi
	echo "Extracted, now modifying...."
	scripts/edit_safemir.sh
	if [ "$?" -eq 1 ]; then
	    (>&2 echo "Error modifying safemir, project probably won't compile")
	    exit 1
	fi
    fi
fi

echo "Done installing externals, configuring...."
./configure.sh
