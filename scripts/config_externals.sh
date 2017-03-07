# Configuration subroutines for Camargue
# The script in this file will test for the presence of the following external
# includes:
# - Catch
# - safeGMI
# - OpenMP
# For each of these, it will modify camargue/includes/config.hpp to add or
# remove the appropriate preprocessor compilation directives. See config.hpp
# for info on what each of these are, and for how to download them if you don't
# have them already! In the case of OpenMP, the Makefile will be modified
# as well.
#
# Finally, this script will attempt to detect whether Concorde has been
# configured with the QSOPT LP Solver, in which case the QSOPT library qsopt.a
# will be added to the Makefile.
#
#
# This script is meant to be invoked by the main config.sh, it will produce
# incorrect results if run on its own.
#
# For compatibility with the Mac version of sed, this script makes the changes
# one at a time by making a duplicate and renaming, rather than in place. This
# is slightly inefficient but not a huge deal realistically.

target=includes/config.hpp


######################### Catch Unit Testing ##################################

printf 'Checking for Catch....'

if [ -f externals/catch.hpp ]
then
    printf 'yes\n'
    sed 's/#undef CMR_HAVE_CATCH/#define CMR_HAVE_CATCH 1/' "$target" \
	> includes/cfg_catch
else
    printf 'no\n'
    sed 's/#define CMR_HAVE_CATCH 1/#undef CMR_HAVE_CATCH/' "$target" \
	> includes/cfg_catch
fi

if [ "$?" -eq 0 ]
then
    mv includes/cfg_catch "$target"
else
    test -e includes/cfg_catch && rm includes/cfg_catch
    (>&2 echo "Error modifying for Catch")
    exit 1
fi

############################## Safe MIR ########################################

printf 'Checking for Safe GMI....'

if [ -e externals/safemir ]
then
    printf 'yes\n'
    sed 's/#undef CMR_HAVE_SAFEGMI/#define CMR_HAVE_SAFEGMI 1/' "$target" \
	> includes/cfg_gmi
else
    printf 'no\n'
    sed 's/#define CMR_HAVE_SAFEGMI 1/#undef CMR_HAVE_SAFEGMI/' "$target" \
	> includes/cfg_gmi
fi

if [ "$?" -eq 0 ]
then
    mv includes/cfg_gmi "$target"
else
    test -e includes/cfg_gmi && rm includes/cfg_gmi
    (>&2 echo "Error modifying for Safe MIR")
    exit 1
fi

############################## OpenMP ########################################

scripts/check_omp.sh
omp_exit="$?"
cfg_exit=1
make_exit=1

if [ "$omp_exit" -eq 1 ]
then
    echo "Fatal error configuring compiler"
    exit 1
fi

if [ "$omp_exit" -eq 0 ]
then
    sed 's/#undef CMR_HAVE_OPENMP/#define CMR_HAVE_OPENMP 1/' "$target" \
	> includes/cfg_omp
    cfg_exit="$?"
    sed 's/\(FOMP *:=\).*/\1 -fopenmp/' Makefile > tmp_make
    make_exit="$?"
else
    echo "Compiler does not appear to support OMP"
    sed 's/#define CMR_HAVE_OPENMP 1/#undef CMR_HAVE_OPENMP/' "$target" \
	> includes/cfg_omp
    cfg_exit="$?"
    sed 's/\(FOMP *:=\).*/\1/' Makefile > tmp_make
    make_exit="$?"
fi

if [ "$make_exit" -eq 0 ]
then
    mv tmp_make Makefile
else
    test -e tmp_make && rm tmp_make
    test -e includes/cfg_omp && rm includes/cfg_omp
    (>&2 echo "Error modifying Makefile for OpenMP")
    exit 1
fi

if [ "$cfg_exit" -eq 0 ]
then
    mv includes/cfg_omp "$target"
else
    test -e includes/cfg_omp && rm includes/cfg_omp
    (>&2 echo "Error modifying config.hpp for OpenMP")
    exit 1
fi


####################### Concorde LP Solver Detection ###########################

echo 'Checking for Concorde QSOPT install....'

CC_LP=$(grep 'LPSOLVER_LIB =' externals/concorde/TSP/Makefile | \
	    sed 's/LPSOLVER_LIB = //')

qsopt=$(cat "$CC_LP" | grep 'qsopt' | wc -l)

if [ "$qsopt" -eq 0 ]; then
    echo "Concorde doesn't seem to be using QSOPT, no changes needed"
    exit 0;
fi

printf 'QSOPT install detected, editing Makefile....'

sed "s@\(QSOPT *:= \)@\1$qsopt@" > tmp_make_qsopt && mv tmp_make_qsopt Makefile

if [ "$?" -ne 0 ]; then
    (>&2 echo "Error modifying Makefile for QSOPT Concorde")
    test -e tmp_make_qsopt && rm tmp_make_qsopt
    exit 1
fi

printf 'done\n'


exit 0
