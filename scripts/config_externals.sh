# Configuration subroutines for Camargue
# The script in this file will test for the presence of the following external
# includes:
# - Catch
# - Timsort
# - safeGMI
# - OpenMP
# For each of these, it will modify camargue/includes/config.hpp to add or
# remove the appropriate preprocessor compilation directives. See config.hpp
# for info on what each of these are, and for how to download them if you don't
# have them already! In the case of OpenMP, the Makefile will be modified
# as well. 
#
# This script is meant to be invoked by the main config.sh, it will produce
# incorrect results if run on its own.
#
# For compatibility with the Mac version of sed, this script makes the changes
# one at a time by making a duplicate and renaming, rather than in place. This
# is slightly inefficient but not a huge deal realistically. 

target=includes/config.hpp


######################### Catch Unit Testing ##################################

have_catch=$(ls externals | grep 'catch.hpp' | wc -l)
echo Checking for Catch: "$have_catch"

if [ "$have_catch" -eq 1 ]
then
    sed 's/#undef CMR_HAVE_CATCH/#define CMR_HAVE_CATCH 1/' "$target" \
	> includes/cfg_catch
else
    sed 's/#define CMR_HAVE_CATCH 1/#undef CMR_HAVE_CATCH/' "$target" \
	> includes/cfg_catch
fi

if [ "$?" -eq 0 ]
then
    mv includes/cfg_catch "$target"
else
    rm includes/cfg_catch
    echo "Error modifying for Catch"
    exit 1
fi


########################## Timsort Implementation #############################

have_tim=$(ls externals | grep 'timsort.hpp' | wc -l)
echo Checking for Timsort: "$have_tim"

if [ "$have_tim" -eq 1 ]
then
    sed 's/#undef CMR_HAVE_TIMSORT/#define CMR_HAVE_TIMSORT 1/' "$target" \
	> includes/cfg_tim
else
    sed 's/#define CMR_HAVE_TIMSORT 1/#undef CMR_HAVE_TIMSORT/' "$target" \
	> includes/cfg_tim
fi

if [ "$?" -eq 0 ]
then
    mv includes/cfg_tim "$target"
else
    rm includes/cfg_tim
    echo "Error modifying for Timsort"
    exit 1
fi


############################## Safe MIR ########################################

have_gmi=$(ls externals | grep 'safemir' | wc -l)
echo Checking for Safe MIR: "$have_gmi"

if [ "$have_gmi" -eq 1 ]
then
    sed 's/#undef CMR_HAVE_SAFEGMI/#define CMR_HAVE_SAFEGMI 1/' "$target" \
	> includes/cfg_gmi
else
    sed 's/#define CMR_HAVE_SAFEGMI 1/#undef CMR_HAVE_SAFEGMI/' "$target" \
	> includes/cfg_gmi
fi

if [ "$?" -eq 0 ]
then
    mv includes/cfg_gmi "$target"
else
    rm includes/cfg_gmi
    echo "Error modifying for Safe MIR"
    exit 1
fi

############################## OpenMP ########################################

scripts/check_omp.sh
omp_exit="$?"
cfg_exit=1
make_exit=1

if [ "$omp_exit" -eq 0 ]
then
    sed 's/#undef CMR_HAVE_OPENMP/#define CMR_HAVE_OPENMP 1/' "$target" \
	> includes/cfg_omp
    cfg_exit="$?"
    sed 's/\(FOMP *:=\).*/\1 -fopenmp/' Makefile > tmp_make
    make_exit="$?"
else
    echo Compiler does not appear to support OMP
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
    rm tmp_make
    rm includes/cfg_omp
    echo "Error modifying Makefile for OpenMP"
    exit 1
fi

if [ "$cfg_exit" -eq 0 ]
then
    mv includes/cfg_omp "$target"
else
    rm includes/cfg_omp
    echo "Error modifying config.hpp for OpenMP"
    exit 1
fi



exit 0

